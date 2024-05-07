#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 16:55:13 2024

@author: angeles

"""

#%%

import numpy as np
import os, sys
import nibabel as nib
import argparse
import pandas as pd
import tarfile
import shutil
#%%
# This script reads txt files with beta vals for each voxel inside a mask 
# of significant pvals (cluster corrected)
# and generates nifti files 
  

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
        options.add_argument('-r', '--CCresults',dest="output", action='store', type=str, required=True,
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
        output = os.path.join(tmp, "results")
        betas = os.path.join(output, "betas")
        CCresults = os.path.join(output,'clustercorrectedNN2')
        
    else :

        options = parse()
        tmp = options.workdir
        masks = options.mask
        output  = options.outputs
        CCresults = options.ccresults


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
            'APdose',
            'PANSS_TP'            
            ]


    # Get the list of files in the directory
    results_files = os.listdir(CCresults)
    if not os.path.exists(betas):
        os.makedirs(betas)
    
    for filename in results_files:
        result_nii = nib.load(os.path.join(CCresults,filename))
        result_vol = result_nii.get_fdata()
        dim3d = result_vol.shape
        result_2d = result_vol.reshape(-1,np.prod(dim3d)).T
        idx_pvals = np.where(result_2d)[0]
        
        for seedname in seednames:
            if seedname in filename:
                print(f"Seed name '{seedname}' found in file: {filename}")
                for index, VOI in enumerate(VOIs, start=1):
                   if VOI in filename:
                        print(f"VOI {VOI} found in file: {filename}")
                
                        beta_filepath = os.path.join(secondlevel,'pvals', f'{seedname}')
                        targetbeta = {
                            1: "beta_interaction",
                            2: "beta_DIT",
                            3: "beta_Group",
                            4: "beta_APdose",
                            5: "beta_PANSSTP"
                            }
                        beta_list = []
                        for voxel in range(len(idx_pvals)):
                            filename = f'voxel_{voxel}_{seedname}.txt'
                            df = pd.read_csv(os.path.join(beta_filepath, filename), sep='\t')
                            
                            beta = df[targetbeta[index]].values[0]
                            
                            beta_list.append(beta)

                        betas = np.array(beta_list)
                        
                        beta_vol = result_vol
                        beta_vol = beta_vol.reshape(-1, np.prod(dim3d))
                        beta_vol[:] = 0
                        beta_vol[0,idx_pvals] = betas[:]
                        beta_vol = beta_vol.reshape(dim3d)
                        
                        beta_nii = nib.Nifti1Image(beta_vol, result_nii.affine)

                        # Save as nifti file

                        beta_filename = "longitudinalTRT_seed-" + \
                                    seedname + \
                                    "_space-MNI152_dim-9110991_" + \
                                    VOI + \
                                    "_betas.nii.gz"                            

                        beta_nii.to_filename(os.path.join(betas, beta_filename))


#%%

if __name__ == '__main__':
    main()







def main ():
    #os.environ["ROOTDIR"] = '/Users/brainsur/Desktop'  # seth path
    os.environ["ROOTDIR"] = '/Volumes/TOSHIBA'  # seth path  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir,"striatconnTRT")
        tmp = os.path.join(workdir, "secondlevel")
        output = os.path.join(tmp, "results")
                
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

