"""Additional exploratory checks (not independent biological validation)."""
from pathlib import Path
import argparse
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score,roc_auc_score
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_analysis import read_data,qc,normalize,select_hvg,cluster,graph_edges,SAMPLES,SEED
from PIL import Image

parser=argparse.ArgumentParser()
parser.add_argument('--input-dir',required=True)
parser.add_argument('--output-dir',required=True)
args=parser.parse_args()
ROOT=Path(args.input_dir);OUT=Path(args.output_dir)
OUT.mkdir(parents=True,exist_ok=True)
rows=[]; marker_rows=[]
for sample,(h,t) in SAMPLES.items():
 print('Extended',sample,flush=True)
 X,barcodes,genes,types,coords,scales,thumb=read_data(ROOT/h,ROOT/t)
 counts,detected,mito=qc(X,genes)
 hvg=select_hvg(X,counts)
 baseline,pcs,evr=cluster(X,counts,hvg,4)
 # Same normalized expression, but regress depth nuisance from every gene (stress test, not ground truth)
 Y=normalize(X,counts,hvg).toarray().astype('float32')
 x=np.log1p(counts).astype('float32'); x-=x.mean(); Y-=Y.mean(axis=0,keepdims=True)
 betas=x@Y/(x@x)
 residual=Y-x[:,None]*betas[None,:]
 pca=PCA(n_components=20,svd_solver='randomized',random_state=SEED,iterated_power=2)
 respcs=pca.fit_transform(residual)
 reslab=KMeans(n_clusters=4,n_init=8,random_state=SEED).fit_predict(respcs[:,:15])
 ar_depth=adjusted_rand_score(baseline,reslab)
 k3=KMeans(n_clusters=3,n_init=8,random_state=SEED).fit_predict(pcs[:,:15])
 k5=KMeans(n_clusters=5,n_init=8,random_state=SEED).fit_predict(pcs[:,:15])
 edges,nn=graph_edges(coords)
 a,b=edges.T
 crossing=(baseline[a]!=baseline[b]).astype(int)
 local_depth=np.abs(np.log1p(counts[a])-np.log1p(counts[b]))
 auroc=float(roc_auc_score(crossing,local_depth)) if len(np.unique(crossing))>1 else np.nan
 # eta squared log-depth by cluster, descriptive only
 logd=np.log1p(counts)
 sst=np.sum((logd-logd.mean())**2)
 ssb=sum(np.sum(baseline==k)*(logd[baseline==k].mean()-logd.mean())**2 for k in range(4))
 eta=float(ssb/sst)
 rows.append(dict(sample=sample, depth_regressed_ARI_vs_baseline=ar_depth,k3_vs_k4_ARI=adjusted_rand_score(k3,baseline),k5_vs_k4_ARI=adjusted_rand_score(k5,baseline),depth_difference_boundary_AUROC=auroc,log_depth_eta_squared_by_clusters=eta,mean_boundary_depth_jump=float(np.mean(local_depth[crossing==1])),mean_nonboundary_depth_jump=float(np.mean(local_depth[crossing==0]))))
 # Normalized exploratory gene overlays, not histopathology labels or markers validated as biomarkers
 candidates=['EPCAM','KRT19','PTPRC','COL1A1','SFTPC','SCGB1A1','MKI67']
 available=[(g,int(np.where(genes==g)[0][0])) for g in candidates if g in genes]
 xcoords=coords['pxl_col_in_fullres'].to_numpy()*scales['tissue_lowres_scalef']
 ycoords=coords['pxl_row_in_fullres'].to_numpy()*scales['tissue_lowres_scalef']
 fig,axes=plt.subplots(1,3,figsize=(13,5),constrained_layout=True)
 for ax,(gene,ix) in zip(axes,available[:3]):
  vec=np.asarray(X[:,ix].toarray()).ravel().astype('float64')
  vals=np.log1p(1e4*vec/np.maximum(counts,1))
  ax.imshow(thumb)
  im=ax.scatter(xcoords,ycoords,c=vals,s=5,cmap='viridis',edgecolors='none')
  ax.set_xlim(0,thumb.width);ax.set_ylim(thumb.height,0);ax.set_title(f'{gene} (log1p CPM)');ax.axis('off')
  fig.colorbar(im,ax=ax,shrink=.65)
  marker_rows.append(dict(sample=sample,gene=gene,spots_detected=int(np.sum(vec>0)),median_logCPM=float(np.median(vals))))
 fig.suptitle(f'{sample} | Exploratory gene expression spatial overlays (not validated biomarkers)',fontsize=13)
 fig.savefig(OUT/f'{sample}_gene_overlay.png',dpi=165);plt.close(fig)

df=pd.DataFrame(rows);df.to_csv(OUT/'additional_sensitivity_results.csv',index=False)
pd.DataFrame(marker_rows).to_csv(OUT/'gene_overlay_descriptives.csv',index=False)
base=pd.read_csv(OUT/'perturbation_runs.csv')
agg=base.groupby(['sample','scenario'],as_index=False)['ARI_vs_baseline'].agg(['mean','min','max']).reset_index()
fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
for ax,s in zip(axes,['P5','P12']):
 d=base[base['sample']==s]
 means=d.groupby('scenario').ARI_vs_baseline.mean().reindex(['half','quarter','equalise'])
 mins=d.groupby('scenario').ARI_vs_baseline.min().reindex(['half','quarter','equalise'])
 maxs=d.groupby('scenario').ARI_vs_baseline.max().reindex(['half','quarter','equalise'])
 positions=np.arange(3)
 ax.bar(positions,means.values,alpha=.75)
 ax.errorbar(positions,means.values,yerr=np.array([means.values-mins.values,maxs.values-means.values]),fmt='none',capsize=4,color='black')
 ax.set_xticks(positions,['50% counts','25% counts','Equal-depth\nthinning']);ax.set_ylim(0,1.05);ax.set_ylabel('ARI vs original clustering');ax.set_title(s+' — count-depth sensitivity')
 ax.grid(axis='y',alpha=.2)
fig.suptitle('Cluster stability decreases under count thinning | Three independent simulations per setting')
fig.savefig(OUT/'cross_sample_thinning_stability.png',dpi=200);plt.close(fig)
print(df.to_string(index=False),flush=True)
