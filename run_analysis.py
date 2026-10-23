#!/usr/bin/env python3
"""Reproducible descriptive NSCLC Visium technical-robustness pilot.
No claim of validated artifact detection or inferentially calibrated biological boundaries.
Usage: python run_analysis.py --input-dir /mnt/data --output-dir /mnt/data/DKFZ_NSCLC_spatial_pilot
"""
from __future__ import annotations
import argparse, io, json, tarfile, warnings, time
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

SEED=26023
RNG=np.random.default_rng(SEED)
SAMPLES={
 'P5': ('GSM8855707_NSCLC_P5_filtered_feature_bc_matrix.h5', 'GSM8855707_NSCLC_P5_spatial.tar.gz'),
 'P12': ('GSM8855714_NSCLC_P12_filtered_feature_bc_matrix.h5', 'GSM8855714_NSCLC_P12_spatial.tar.gz'),
}

def read_data(h5path,tarpath):
 with h5py.File(h5path,'r') as f:
  m=f['matrix']; shape=tuple(m['shape'][:]); barcodes=np.char.decode(m['barcodes'][:]); genes=np.char.decode(m['features/name'][:]); feature_types=np.char.decode(m['features/feature_type'][:]);
  X=sp.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=shape).T.tocsr().astype(np.int32)
 # Keep only RNA Gene Expression features; exclude 35 Antibody Capture features.
 feature_mask=(feature_types=='Gene Expression')
 X=X[:,feature_mask].tocsr()
 genes=genes[feature_mask]
 feature_types=feature_types[feature_mask]
 with tarfile.open(tarpath,'r:gz') as t:
  def load_suffix(s):
   f=next((f for f in t if f.name.endswith('/'+s)),None)
   if f is None: raise ValueError('Spatial archive missing '+s)
   return t.extractfile(f).read()
  coords=pd.read_csv(io.BytesIO(load_suffix('tissue_positions.csv'))).set_index('barcode')
  scales=json.loads(load_suffix('scalefactors_json.json'))
  thumb=Image.open(io.BytesIO(load_suffix('tissue_lowres_image.png'))).convert('RGB')
 assert coords.index.is_unique and len(set(barcodes))==len(barcodes)
 assert set(barcodes)==set(coords.index[coords['in_tissue']==1]), 'Expression-spatial barcode mismatch'
 coords=coords.loc[barcodes].copy()
 assert bool(coords['in_tissue'].eq(1).all())
 return X,barcodes,genes,feature_types,coords,scales,thumb

def qc(X,genes):
 depth=np.asarray(X.sum(axis=1)).ravel().astype('float64')
 detected=np.diff(X.indptr)
 mito= np.char.startswith(np.char.upper(genes.astype(str)), 'MT-')
 pct_mt=(np.asarray(X[:,mito].sum(axis=1)).ravel() / np.maximum(depth,1))*100 if mito.any() else np.full(len(depth),np.nan)
 return depth,detected,pct_mt

def normalize(X,depth,genes_idx=None):
 if genes_idx is not None: X=X[:,genes_idx]
 # library size is based on ALL genes, not just the selected HVGs
 X=X.tocsr().astype('float32',copy=True)
 sf=(1e4/np.maximum(depth,1)).astype('float32')
 X=sp.diags(sf,format='csr') @ X
 X.data=np.log1p(X.data)
 return X

def select_hvg(X,depth,n=1200):
 Y=normalize(X,depth)
 mean=np.asarray(Y.mean(axis=0)).ravel()
 squaremean=np.asarray(Y.power(2).mean(axis=0)).ravel()
 var=np.maximum(0,squaremean-mean**2)
 # select by log-normalized variance after minimum detection; fixed for perturbation tests
 detected=np.asarray(X.getnnz(axis=0)).ravel()
 keep=(detected>=max(15,round(X.shape[0]*0.02))) & (mean>0.01)
 indices=np.where(keep)[0]
 return indices[np.argsort(var[indices])[-min(n,len(indices)):]]

def cluster(X,depth,hvg,k,seed=SEED):
 Z=normalize(X,depth,hvg).toarray()
 Z-=Z.mean(axis=0,keepdims=True)
 pca=PCA(n_components=20,svd_solver='randomized',random_state=seed,iterated_power=2)
 pcs=pca.fit_transform(Z)
 labels=KMeans(n_clusters=k,n_init=8,random_state=seed).fit_predict(pcs[:, :15])
 return labels,pcs,pca.explained_variance_ratio_

def match_labels(reference, candidate,k):
 confusion=np.zeros((k,k),dtype='int64')
 np.add.at(confusion,(reference,candidate),1)
 rr,cc=linear_sum_assignment(-confusion)
 mapper=dict(zip(cc,rr))
 return np.array([mapper[v] for v in candidate],dtype='int32')

def graph_edges(coords):
 # Use full-resolution tissue pixel coordinates; true Visium spot centers.
 pts=coords[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy().astype(float)
 tree=cKDTree(pts)
 distances,indices=tree.query(pts,k=7)
 nearest=np.median(distances[:,1])
 # local graph, cap long distances across tissue holes
 rows=np.repeat(np.arange(len(pts)),6); cols=indices[:,1:].ravel(); dist=distances[:,1:].ravel()
 good=(dist<nearest*1.30)
 a=np.minimum(rows[good],cols[good]); b=np.maximum(rows[good],cols[good]);
 edges=np.unique(np.stack((a,b),axis=1),axis=0)
 return edges,nearest

def moran(values,edges):
 centered=values-values.mean()
 denom=np.dot(centered,centered)
 if denom==0:return np.nan
 return len(values)/len(edges) * np.sum(centered[edges[:,0]]*centered[edges[:,1]]) / denom

def block_labels(coords,nblocks=4):
 x=coords['pxl_col_in_fullres'].to_numpy().astype(float)
 y=coords['pxl_row_in_fullres'].to_numpy().astype(float)
 bx=np.minimum(nblocks-1,np.floor((x-x.min())/(np.ptp(x)+1e-12)*nblocks).astype(int))
 by=np.minimum(nblocks-1,np.floor((y-y.min())/(np.ptp(y)+1e-12)*nblocks).astype(int))
 return bx+nblocks*by

def permutation_stats(labels,depth,edges,blocks,nperm=199):
 observed=np.mean(labels[edges[:,0]]==labels[edges[:,1]])
 rng=np.random.default_rng(SEED+122)
 global_null=[]; block_null=[]
 unique_blocks=[np.where(blocks==b)[0] for b in np.unique(blocks)]
 for i in range(nperm):
  lab=rng.permutation(labels)
  global_null.append(np.mean(lab[edges[:,0]]==lab[edges[:,1]]))
  lab2=labels.copy()
  for ix in unique_blocks:
   if len(ix)>1:lab2[ix]=lab2[rng.permutation(ix)]
  block_null.append(np.mean(lab2[edges[:,0]]==lab2[edges[:,1]]))
 m=moran(np.log1p(depth),edges)
 depth_null=[]
 for i in range(nperm): depth_null.append(moran(rng.permutation(np.log1p(depth)),edges))
 return dict(edge_coherence=float(observed), global_null_mean=float(np.mean(global_null)),
  block_null_mean=float(np.mean(block_null)),block_null_q95=float(np.quantile(block_null,.95)),
  global_shuffle_p=(1+np.sum(np.array(global_null)>=observed))/(nperm+1),
  block_shuffle_p=(1+np.sum(np.array(block_null)>=observed))/(nperm+1),
  depth_moran_I=float(m), depth_naive_permutation_p=(1+np.sum(np.array(depth_null)>=m))/(nperm+1))

def depth_thinning(X,depth,mode,seed):
 rng=np.random.default_rng(seed)
 if mode=='half': q=np.full(len(depth),0.5,dtype='float32')
 elif mode=='quarter':q=np.full(len(depth),0.25,dtype='float32')
 elif mode=='equalise':
  target=np.quantile(depth,0.25)
  q=np.minimum(1,target/np.maximum(depth,1)).astype('float32')
 else:raise ValueError(mode)
 Y=X.copy()
 rep=np.repeat(np.arange(X.shape[0]),np.diff(Y.indptr))
 Y.data=rng.binomial(Y.data,q[rep]).astype(np.int32)
 Y.eliminate_zeros()
 d=np.asarray(Y.sum(axis=1)).ravel()
 return Y,d

def visualise(sample,out,xy,thumb,scale,depth,detected,base_labels,edge_support,edges,pcs,k):
 x=xy['pxl_col_in_fullres'].to_numpy()*scale['tissue_lowres_scalef'];y=xy['pxl_row_in_fullres'].to_numpy()*scale['tissue_lowres_scalef']
 fig,axes=plt.subplots(2,3,figsize=(16,10),constrained_layout=True)
 for ax in axes.flat:
  ax.imshow(thumb);ax.set_xlim(0,thumb.width);ax.set_ylim(thumb.height,0);ax.axis('off')
 im=axes[0,0].scatter(x,y,c=np.log10(depth+1),s=5,cmap='viridis',edgecolors='none');axes[0,0].set_title('Log10 total UMI counts');fig.colorbar(im,ax=axes[0,0],shrink=.65)
 im=axes[0,1].scatter(x,y,c=detected,s=5,cmap='magma',edgecolors='none');axes[0,1].set_title('Detected genes');fig.colorbar(im,ax=axes[0,1],shrink=.65)
 cmap=plt.get_cmap('tab10',k)
 axes[0,2].scatter(x,y,c=base_labels,s=5,cmap=cmap,vmin=-.5,vmax=k-.5,edgecolors='none');axes[0,2].set_title(f'Expression clustering (k={k}; NOT histology)')
 # Edge consensus: mean support that corresponding pair is assigned different label under perturbations
 node_support=np.bincount(edges.ravel(),weights=np.repeat(edge_support,2),minlength=len(depth))/np.maximum(np.bincount(edges.ravel(),minlength=len(depth)),1)
 im=axes[1,0].scatter(x,y,c=node_support,s=6,cmap='coolwarm',vmin=0,vmax=1,edgecolors='none');axes[1,0].set_title('Fraction of perturbed runs with neighbor disagreement');fig.colorbar(im,ax=axes[1,0],shrink=.65)
 axes[1,1].axis('on');axes[1,1].cla();axes[1,1].scatter(pcs[:,0],pcs[:,1],c=base_labels,cmap=cmap,s=4,alpha=.7,edgecolors='none');axes[1,1].set_title('PCA of expression');axes[1,1].set_aspect('auto');axes[1,1].set_xlabel('PC1');axes[1,1].set_ylabel('PC2')
 axes[1,2].axis('on');axes[1,2].cla();axes[1,2].scatter(np.log10(depth+1),pcs[:,0],s=4,alpha=.35);axes[1,2].set_aspect('auto');axes[1,2].set_xlabel('Log10 total UMI counts');axes[1,2].set_ylabel('PC1');axes[1,2].set_title('Depth–PC1 association (descriptive)')
 fig.suptitle(f'{sample} | Lung NSCLC Visium | Single specimen, observational pilot',fontsize=17)
 fig.savefig(out / f'{sample}_spatial_analysis.png',dpi=175)
 plt.close(fig)
 return node_support

def run_sample(sample,paths,out):
 tick=time.time(); print(f'[{sample}] loading',flush=True)
 X,barcodes,genes,feature_types,coords,scale,thumb=read_data(*paths)
 # matrix contains Gene Expression and possibly controls; exclude non-gene expression if present
 # Keep exact matrix for traceability; use all genes in sums but HVGs filter low detection.
 counts,detected,mito=qc(X,genes)
 hvg=select_hvg(X,counts)
 edges,nn=graph_edges(coords)
 blocks=block_labels(coords)
 print(f'[{sample}] spots={X.shape[0]} genes={X.shape[1]} edges={len(edges)} mean_depth={counts.mean():.0f} median_depth={np.median(counts):.0f}; PCA/clustering',flush=True)
 k=4
 base,pcs,evr=cluster(X,counts,hvg,k)
 aris=[];runrecords=[];all_labels=[]
 modes=['half','quarter','equalise']
 for mode in modes:
  for rep in range(3):
   pert,nd=depth_thinning(X,counts,mode,SEED+rep+{'half':100,'quarter':300,'equalise':500}[mode])
   labels,_,_=cluster(pert,nd,hvg,k,seed=SEED+rep+{'half':100,'quarter':300,'equalise':500}[mode])
   ar=adjusted_rand_score(base,labels)
   aligned=match_labels(base,labels,k)
   all_labels.append(aligned)
   runrecords.append(dict(sample=sample,scenario=mode,replicate=rep+1,ARI_vs_baseline=float(ar),median_UMI_after=float(np.median(nd)),fraction_counts_retained=float(nd.sum()/counts.sum())))
   print(f'[{sample}] mode={mode} rep={rep+1} ARI={ar:.3f}',flush=True)
 all_labels=np.asarray(all_labels)
 edge_support=np.mean(all_labels[:,edges[:,0]]!=all_labels[:,edges[:,1]],axis=0)
 spot_concord=np.mean(all_labels==base[None,:],axis=0)
 sstats=permutation_stats(base,counts,edges,blocks,nperm=199)
 sstats['spearman_pc1_log_counts']=float(spearmanr(pcs[:,0],np.log1p(counts)).statistic)
 sstats['spot_concordance_median']=float(np.median(spot_concord))
 sstats['spot_concordance_q10']=float(np.quantile(spot_concord,.1))
 sstats['edge_disagreement_mean_perturbed']=float(np.mean(edge_support))
 sstats['baseline_boundary_fraction']=float(np.mean(base[edges[:,0]]!=base[edges[:,1]]))
 sstats['spots']=int(len(counts)); sstats['genes']=int(X.shape[1]); sstats['n_edges']=int(len(edges));sstats['n_hvg']=int(len(hvg));sstats['sample']=sample
 sstats['median_UMI']=float(np.median(counts));sstats['median_detected_genes']=float(np.median(detected));sstats['mitochondrial_genes_count']=int(np.char.startswith(np.char.upper(genes.astype(str)),'MT-').sum());sstats['pc1_variance_explained']=float(evr[0]);sstats['mutual_replication_note']='two independent specimens, separate clustering; no biological group contrast estimable'
 print(f'[{sample}] plotting',flush=True)
 spot_support=visualise(sample,out,coords,thumb,scale,counts,detected,base,edge_support,edges,pcs,k)
 obs=pd.DataFrame(dict(barcode=barcodes,patient=sample,x_fullres=coords['pxl_col_in_fullres'].to_numpy(),y_fullres=coords['pxl_row_in_fullres'].to_numpy(),UMI=counts,detected_genes=detected,mitochondrial_pct=mito,cluster_k4=base,PC1=pcs[:,0],PC2=pcs[:,1],spot_cluster_concordance=spot_concord,neighbor_disagreement_perturbed=spot_support))
 obs.to_csv(out / f'{sample}_spot_level_metrics.csv.gz',index=False,compression='gzip')
 pd.DataFrame(dict(gene=genes[hvg])).to_csv(out / f'{sample}_selected_HVGs.csv',index=False)
 print(f'[{sample}] done {time.time()-tick:.1f}s',flush=True)
 return sstats,runrecords

def main():
 parser=argparse.ArgumentParser()
 parser.add_argument('--input-dir',default='/mnt/data')
 parser.add_argument('--output-dir',default='/mnt/data/DKFZ_NSCLC_spatial_pilot')
 args=parser.parse_args();inp=Path(args.input_dir); out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
 rows=[];reps=[]
 for sample,(h,t) in SAMPLES.items():
  s,r=run_sample(sample,(inp/h,inp/t),out);rows.append(s);reps+=r
 pd.DataFrame(rows).to_csv(out/'sample_results.csv',index=False)
 pd.DataFrame(reps).to_csv(out/'perturbation_runs.csv',index=False)
 print('ALL DONE',flush=True)
if __name__=='__main__':main()
