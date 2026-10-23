#!/usr/bin/env python3
"""Plot descriptive local cluster-boundary perturbation maps from previous analysis output.
This script does NOT estimate posterior P(biological boundary | data)."""
from pathlib import Path
import argparse
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

def run(sample,results,out):
 df=pd.read_csv(results/f'{sample}_spot_level_metrics.csv.gz')
 xy=df[['x_fullres','y_fullres']].to_numpy(dtype=float)
 labels=df.cluster_k4.to_numpy().astype(int)
 score=np.clip(df.neighbor_disagreement_perturbed.to_numpy(dtype=float),0,1)
 dist,idx=cKDTree(xy).query(xy,k=7)
 near=np.median(dist[:,1]); rows=np.repeat(np.arange(len(xy)),6)
 cols=idx[:,1:].ravel(); distances=dist[:,1:].ravel(); mask=distances<1.3*near
 edges=np.unique(np.stack([np.minimum(rows[mask],cols[mask]),np.maximum(rows[mask],cols[mask])],axis=1),axis=0)
 crossing=edges[labels[edges[:,0]]!=labels[edges[:,1]]]
 proxy=(score[crossing[:,0]]+score[crossing[:,1]])/2
 selected=crossing[proxy>=np.quantile(proxy,.75)]
 fig,axs=plt.subplots(1,2,figsize=(12.3,5.5),constrained_layout=True)
 for ax in axs:
  ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('Full-resolution tissue x (pixels)')
 axs[0].scatter(xy[:,0],xy[:,1],c=labels,cmap='tab10',s=8,alpha=.93,linewidths=0,vmin=-.5,vmax=3.5)
 axs[0].set_title('A. Expression-derived clusters (k = 4)')
 sc=axs[1].scatter(xy[:,0],xy[:,1],c=score,cmap='viridis',vmin=0,vmax=1,s=9,linewidths=0)
 for i,j in selected:axs[1].plot(xy[[i,j],0],xy[[i,j],1],color='black',linewidth=.48,alpha=.47)
 axs[1].set_title('B. Local perturbation disagreement + boundary edges')
 fig.colorbar(sc,ax=axs[1],shrink=.72,label='Descriptive disagreement frequency (0–1)')
 fig.suptitle(f'{sample} | Public NSCLC Visium tissue: exploratory boundary sensitivity',fontsize=14)
 fig.text(.5,-.01,f'{len(df):,} spots; {len(crossing):,} baseline crossing edges. Lines = top quartile of local boundary-adjacent disagreement. Not calibrated biological probabilities.',ha='center',fontsize=8)
 fig.savefig(out/f'{sample}_boundary_perturbation_map.png',dpi=190,bbox_inches='tight');plt.close(fig)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--results-dir',default='results');p.add_argument('--output-dir',default='figures');a=p.parse_args()
 out=Path(a.output_dir);out.mkdir(exist_ok=True)
 for sample in ['P5','P12']:run(sample,Path(a.results_dir),out)
