#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 15:27:41 2024

@author: angeles
"""

#%%
#%%

import os, sys
import argparse

import numpy as np
import nibabel as nib
from scipy.ndimage import label

#%%
def parse():

        options = argparse.ArgumentParser(description=" Created by ...")
        options.add_argument('-c', '--minclustersize',dest="cluster_thr", action='store', type=float, required=True,
                             help='min cluster size for cluster correction')
        options.add_argument('-o', '--output',dest="outputs", action='store', type=str, required=True,
                             help='path to saved niftis from second level')
        #print(options.parse_args())

        return options.parse_args()

#%%

def main ():
    os.environ["ROOTDIR"] = '/Users/brainsur/Desktop'  # seth path
    #os.environ["ROOTDIR"] = '/Volumes/TOSHIBA'  # seth path  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir,"striatconnTRT")
        tmp = os.path.join(workdir, "secondlevel")
        output = os.path.join(tmp, "results")
        cluster_thr = 363 #FROM 3DCLUSTSIZE (OUTPUT OF CLUSTERCORRECTION_PART1.PY)
        
        
    else :

        options = parse()
        output  = options.outputs
        cluster_thr = options.clustersize


    seednames = ['DCPutamen',
                 'DorsalCaudate',
                 'DRPutamen',
                 'InfVentralCaudate',
                 'SupVentralCaudate',
                 'VRPutamen' #_space-MNI152NLin2009cAsym.nii.gz
                 ]


    #################
  
    for seed in range(len(seednames)):
        
        for _ in range(5):
            # Dictionary mapping filenames to numbers
            filenames = {
                  1: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_interactionDITxgroup_1-pvals-uncorrected.nii.gz",
                  2: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_DIT_1-pvals-uncorrected.nii.gz",
                  3: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_Group_1-pvals-uncorrected.nii.gz",
                  4: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_APdose_1-pvals-uncorrected.nii.gz",
                  5: "longitudinalTRT_seed-" + seednames[seed] + "_space-MNI152_dim-9110991_PANSS_TP_1-pvals-uncorrected.nii.gz"
            }
        
            # Choose filename based on loop iteration
            filename = filenames[_ + 1]
            print(f'Working on {filename}')
              
            image_path = os.path.join(output,filename)
          
            # Load the NIfTI image
            nii_img = nib.load(image_path)
            data = nii_img.get_fdata()
    
            # Inverse p-values
            inverse_p_values = 1 - data

            # Threshold for p < 0.01
            thresholded_data = inverse_p_values < 0.01

            # NN3 Perform clustering using third-nearest neighbor criterion NN=3
            #clusters, num_clusters = label(thresholded_data, structure=np.ones((3, 3, 3)))
           
            
           # NN2 Perform clustering using third-nearest neighbor criterion NN=2
            neighborhood_structure = np.array([[[0, 1, 0], [1, 0, 1], [0, 1, 0]],
                                   [[1, 1, 1], [0, 1, 0], [1, 1, 1]],
                                   [[0, 1, 0], [1, 0, 1], [0, 1, 0]]])

            # Perform clustering using second-nearest neighbor criterion
            clusters, num_clusters = label(thresholded_data, structure=neighborhood_structure)
            
            # Find clusters with 364 or more voxels
            cluster_sizes = np.bincount(clusters.flatten())
            large_clusters_indices = np.where(cluster_sizes >= cluster_thr)[0]

            # Create a new image with only clusters > 364 and their corresponding inverse p-values
            new_data = np.zeros_like(data)
            
            if not os.path.exists(os.path.join(output,'clustercorrectedNN2')):
                os.makedirs(os.path.join(output,'clustercorrectedNN2'))

            for idx in large_clusters_indices:
                if idx == 0:
                    continue  # Skip background cluster (label 0)
                cluster_mask = clusters == idx
                new_data[cluster_mask] = inverse_p_values[cluster_mask]

                # Create a new NIfTI image with the modified data
                new_nii_img = nib.Nifti1Image(new_data, affine=nii_img.affine, header=nii_img.header)     
           
                # Save the new NIfTI image
                nib.save(new_nii_img, os.path.join(output,'clustercorrectedNN2',filename))
                
            
#%%

if __name__ == '__main__':
    main()

