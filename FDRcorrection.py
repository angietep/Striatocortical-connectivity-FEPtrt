#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 12:00:05 2024

@author: angeles
"""

# FDR CORRECTION
# 1) Read uncorrected results 
# 2) Remove nans
# 3) Apply FDR

#%%

import os, sys
import argparse

import numpy as np
import nibabel as nib
from statsmodels.stats.multitest import multipletests


#%%
def parse():

        options = argparse.ArgumentParser(description=" Created by ...")
        options.add_argument('-ud', '--uncorrected_data',dest="uncorrected_data", action='store', type=str, required=True,
                             help='path to saved niftis from second level')
        #print(options.parse_args())

        return options.parse_args()

#%%

def main ():
    #os.environ["ROOTDIR"] = '/Users/brainsur/Desktop'  # seth path
    os.environ["ROOTDIR"] = '/Volumes/TOSHIBA'  # seth path  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir,"striatconnTRT")
        tmp = os.path.join(workdir, "secondlevel")
        uncorrected_data = os.path.join(tmp, "results","uncorrected")
        output = os.path.join(uncorrected_data, "../FDRcorrected")
        masks = os.path.join(workdir,"masks")
        
    else :

        options = parse()
        uncorrected_data  = options.uncorrected_data
       

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
        
        for _ in range(7):
            # Dictionary mapping filenames to numbers
            filenames = {
                  1: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_HC_1-pvals-uncorrected.nii.gz",
                  2: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_TRS_1-pvals-uncorrected.nii.gz",
                  3: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_time_1-pvals-uncorrected.nii.gz",
                  4: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_timexHC_1-pvals-uncorrected.nii.gz",
                  5: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_timexTRS_1-pvals-uncorrected.nii.gz",
                  6: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_APdose_1-pvals-uncorrected.nii.gz",
                  7: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_PANSSTP_1-pvals-uncorrected.nii.gz"
            }
        
            # Choose filename based on loop iteration
            filename = filenames[_ + 1]
            print(f'Working on {filename}')
              
            image_path = os.path.join(uncorrected_data,filename)
          
            # Load the NIfTI image
            nii_img = nib.load(image_path)
            data = nii_img.get_fdata()
            dim3d = data.shape
            
            data_2d = data.reshape(-1, np.prod(dim3d)).T 
            p_values = 1 - data_2d[idx_GM] #p-values are inverted in input img
            
            # Remove NaN values
            p_values_no_nan = p_values[~np.isnan(p_values)]
            
            # Perform FDR correction using the Benjamini-Hochberg method
            reject_no_nan, pvals_corrected_no_nan, _, _ = multipletests(p_values_no_nan, alpha=0.05, method='fdr_bh')
            
            # Reinsert NaN values (to keep GM mask size)
            reject = np.full(p_values.shape, np.nan)
            pvals_corrected = np.full(p_values.shape, np.nan)
            reject[~np.isnan(p_values)] = reject_no_nan
            pvals_corrected[~np.isnan(p_values)] = pvals_corrected_no_nan
            
            # Inverse p-values
            inverse_p_values = 1 - pvals_corrected
            #Replace nan
            inverse_p_values[np.isnan(inverse_p_values)] = 0
            reject[np.isnan(reject)] = 0
            
            corrected = data_2d
            corrected[:] = 0 # clean img
            corrected[idx_GM] = inverse_p_values
            corrected = corrected.reshape(dim3d)
            
            if any(reject==True):
           
                # Create a new NIfTI image with the modified data
                new_nii_img = nib.Nifti1Image(corrected, affine=nii_img.affine, header=nii_img.header)     
           
                if not os.path.exists(output):
                    os.makedirs(output)
    
                # Save the new NIfTI image
                nib.save(new_nii_img, os.path.join(output,filename))
                          

#%%

if __name__ == '__main__':
    main()

