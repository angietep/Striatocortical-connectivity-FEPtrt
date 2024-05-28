#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:53:22 2024

@author: brainsur
"""

# This script reads txt-files with p-values of each voxel,
# and generates nifti file with 1-pval on each voxel of the gm mask.

#%%

import numpy as np
import os, sys
import nibabel as nib
import argparse
import pandas as pd
import tarfile
import shutil

#%%
#%%
def parse():

        options = argparse.ArgumentParser(description="Run 2nd level analysis Assemble pvalues into nifti . Created by ...")
        options.add_argument('-p', '--pvalpath',dest="pvalpath", action='store', type=str, required=True,
                             help='path to pvalues')
        options.add_argument('-m', '--mask',dest="masks", action='store', type=str, required=True,
                             help='path to mask')
        options.add_argument('-o', '--output',dest="outputs", action='store', type=str, required=True,
                             help='path to save nifti')
        #print(options.parse_args())

        return options.parse_args()

#%%
def main ():
    os.environ["ROOTDIR"] = '/Users/brainsur/'  # seth path
    #os.environ["ROOTDIR"] = '/Users/angeles/'  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir,"Desktop/striatconnTRT")
        masks = os.path.join(workdir,"masks")
        tmp = os.path.join(workdir, "secondlevel")
        output = os.path.join(tmp, "results")
        
        
    else :

        options = parse()
        #participants = options.participants
        tmp = options.workdir
        masks = options.mask
        #rawdata = options.rawdata
        #derivat = options.derivatives
        output  = options.outputs


    seednames = ['DCPutamen',
                 'DorsalCaudate',
                 'DRPutamen',
                 'InfVentralCaudate',
                 'SupVentralCaudate',
                 'VRPutamen' #_space-MNI152NLin2009cAsym.nii.gz
                 ]


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



    for seed in range(len(seednames)):
        
        print(f"\t ----- Working on seed {seed}: {seednames[seed]} ----")
        #print(f"Extracting folder {seednames[seed]}.tar.gz") I've extracted outside for performance
        # Specify the path to the compressed archive (tar.gz file)
        #compressed_folder = os.path.join(tmp,f'{seednames[seed]}.tar.gz')
        #outputseed = os.path.join(tmp,f'{seednames[seed]}')
        # Extract the contents of the compressed archive
        #with tarfile.open(compressed_folder, 'r:gz') as archive:
        #    archive.extractall(outputseed)
        #os.remove(compressed_folder)
        pval_filepath = os.path.join(tmp,'pvals_suppFig1', f'{seednames[seed]}')

        pval_list=[]
        beta_list=[]
        
        for voxel in range(len(idx_GM)):
            #print(f"Voxel {voxel} out of {len(idx_GM)} for seed {seed}: {seednames[seed]}")

            filename = f'voxel_{voxel}_{seednames[seed]}.txt'
            #If file does not exist. (i.e.: model did not converge) - replace with nans
            if not os.path.isfile(os.path.join(pval_filepath, filename)):
                pval = np.nan
                beta = np.nan
                
            # Read the text file into a DataFrame
            else:
                df = pd.read_csv(os.path.join(pval_filepath, filename), sep='\t')
                
                # Access the p-values from the DataFrame
                pval = df["pval_intercept"].values[0]
                if pval < 0.01:
                    beta = df["beta_intercept"].values[0]
                else: beta = np.nan
                
            #pval_list.append(pval)
            beta_list.append(beta)

        # Convert p-values to a NumPy array and invert (1-pval)
        
        #pvalues = np.array(pval_list)
        #inv_pvals = 1 - pvalues
        
        betas = np.array(beta_list)
              

        #Generate nifti files for SEED p-values using brainmask affine info
        #HC PVAL
        seedmap_vol1 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol1 = seedmap_vol1.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol1[:] = 0 # clean img
        seedmap_vol1[0,idx_GM] = betas[:]
        seedmap_vol1 = seedmap_vol1.reshape(dim3d) #reshape to 3D


        #Generate nifti objects
        seedmap_nii1 = nib.Nifti1Image(seedmap_vol1, Vgm_nii.affine)
   

        # Save as nifti files

        filename1 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_suppFig1_betas-for-pvalsunder001-uncorrected.nii.gz"
   

        if not os.path.exists(output):
            os.makedirs(output)

        seedmap_nii1.to_filename(os.path.join(output, filename1))
      

        print(f"\t ----- Finished seed {seed}: {seednames[seed]} ----")

        # # Create the compressed archive
        # print(f'Compressing folder {outputseed}')
        # shutil.make_archive(outputseed, 'gztar',outputseed)
        # # Remove the original folder
        # print(f'Folder compressed.\nRemoving original folder {outputseed}')
        # shutil.rmtree(outputseed)


#%%

if __name__ == '__main__':
    main()








