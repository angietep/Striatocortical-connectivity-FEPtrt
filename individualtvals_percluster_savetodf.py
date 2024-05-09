#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 12:32:02 2024

@author: brainsur
#%%
# This script reads txt files with t-vals for each subject and 
# for each voxel inside a mask 
# of significant pvals (cluster corrected)
# and generates a column in a df to later plot interaction plots
# (in another script)
  
"""

#%%

import numpy as np
import os, sys
import nibabel as nib
import argparse
import pandas as pd
from scipy.ndimage import label

#%%
def parse():

        options = argparse.ArgumentParser(description="Run 1st level analysis. Created by ...")
        options.add_argument('-v', '--voxel',dest="voxel", action='store', type=str, required=False,
                            help='voxel file')
        options.add_argument('-c', '--covariates',dest="covariates", action='store', type=str, required=False,
                            help='covariates file')
        options.set_defaults(covariates=None)

        options.add_argument('-t', '--tvalpath',dest="tvalpath", action='store', type=str, required=False,
                             help='the work directory for the tvals')
        
     
        options.add_argument('-o', '--output',dest="output", action='store', type=str, required=True,
                            help='path to output directory')
        options.add_argument('-r', '--betas_path',dest="betas_path", action='store', type=str, required=True,
                            help='path to clustercorrected results')

        #print(options.parse_args())

        return options.parse_args()


#%%
"""
pseudocódigo

Leer clustercorrected output (CCoutput)
Detectar seedname from filename
Detectar variable of interest from filename (VOI)
Encontrar pval folder for that seed
Usar CCoutput como mask
Leer beta of VOI from pval files in voxels inside mask
Armar nueva nifti

"""

#%%
def main ():
    os.environ["ROOTDIR"] = '/Users/brainsur/Desktop'  # seth path
    #os.environ["ROOTDIR"] = '/Volumes/TOSHIBA/'  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir,"striatconnTRT")
        masks = os.path.join(workdir,"masks")
        secondlevel = os.path.join(workdir, "secondlevel")
        output = os.path.join(secondlevel, "results")
        betas_path = os.path.join(output, "betas")
        
    else :
        options = parse()
        workdir = options.workdir
        masks = options.mask
        output  = options.outputs
        betas_path = options.betas_path
      

    seednames = ['DCPutamen',
                  'DorsalCaudate',
                  'DRPutamen',
                  'InfVentralCaudate',
                  'SupVentralCaudate',
                  'VRPutamen' #_space-MNI152NLin2009cAsym.nii.gz
                  ]
    
    VOIs = ['interactionDITxgroup',
            'DIT',
            'Group',
            #'APdose',
            #'PANSS_TP'            
            ]
    sample = pd.read_csv(os.path.join(workdir,'cleansample_covars.csv'))
    
    #################
    Vgm_nii = nib.load(os.path.join(masks,'GrayMattermask_thalamus_space-MNI152_dim-9110991.nii.gz'))
    #read vol
    Vgm_vol = Vgm_nii.get_fdata()
    #save origianl dimensions (voxels_x, voxels_y, voxels_z)
    dim3d = Vgm_vol.shape
    #reshape to 2D
    Vgm_2d = Vgm_vol.reshape(-1, np.prod(dim3d)).T  # -1 means auto-calculate size of dimension
    #save indexes in which Vgm == 1 (indexes for gray matter location)
    idx_GM = np.where(Vgm_2d)[0]

    # Get the list of files in the directory
    results_files = os.listdir(betas_path)
       
    for res_filename in results_files:
        result_nii = nib.load(os.path.join(betas_path,res_filename))
        result_vol = result_nii.get_fdata()
        
        # Perform clustering using second-nearest neighbor criterion
        clusters, num_clusters = label(result_vol, structure=np.ones((3,3,3)))
        clusters_2d= clusters.reshape(-1,np.prod(dim3d)).T        
        cluster_size = np.bincount(clusters.flatten())[1]
    
        for seedname in seednames:
            if seedname in res_filename:
                print(f"Seed name '{seedname}' found in file: {res_filename}")
                for index, VOI in enumerate(VOIs, start=1):
                    if VOI in res_filename.split('_'):
                        print(f"VOI {VOI} found in file: {res_filename}")
                
                        tvals_filepath = os.path.join(secondlevel,'tvals', f'{seedname}')
                        
                        idx_cluster = np.where(clusters_2d)[0]
                        all_tvals = []
                        for i in idx_cluster:
                            voxel = np.where(idx_GM == i)[0][0]
                            vox_filename = f'voxel_{voxel}_{seedname}.txt'
                            df = pd.read_csv(os.path.join(tvals_filepath, vox_filename), sep='\t', header = None)
                            
                            all_tvals.append(df[0].values)
                        
                        all_tvals = np.array(all_tvals).T
                        tval_mean = all_tvals.mean(axis=1)
                        
                        colname = f'tvals_{seedname}_{VOI}_{res_filename[-16:-7]}_size_{cluster_size}'
                        sample[colname] = tval_mean
                    
    sample.to_csv(os.path.join(workdir,'cleansample_covars.csv'), index=False)

#%%

if __name__ == '__main__':
    main()


