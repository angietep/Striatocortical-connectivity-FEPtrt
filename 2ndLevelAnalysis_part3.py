# SECOND LEVEL PART 3: This script reads txt-files with p-values of each voxel,
# and generates nifti files with 1-pval on each voxel of the gm mask.

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
    if hasattr(sys, "ps1"):
        options = {}
        tmp = "/data/zugmana2/INPD/derivatives/angeles/INPD_secondlevel-2/tmp/results/" #os.environ["ROOTDIR"]
        masks = os.path.join("/data/zugmana2/INPD","code","Striatocortical-connectivity","masks")
        output = "/data/zugmana2/INPD/derivatives/angeles/INPD_secondlevel-2/results"

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
    Vgm_nii = nib.load(os.path.join(masks,'GrayMattermask_thalamus_space-MNI152NLin2009cAsym.nii.gz'))
    #read vol
    Vgm_vol = Vgm_nii.get_fdata()
    #save origianl dimensions (voxels_x, voxels_y, voxels_z)
    dim3d = Vgm_vol.shape
    #reshape to 2D
    Vgm_2d = Vgm_vol.reshape(-1, np.prod(dim3d)).T  # -1 means auto-calculate size of dimension
    #save indexes in which Vgm == 1 (indexes for gray matter location)
    idx_GM = np.where(Vgm_2d)[0]



    for seed in range(len(seednames)):
        #print(f"Extracting folder {seednames[seed]}.tar.gz") I've extracted outside for performance
        # Specify the path to the compressed archive (tar.gz file)
        #compressed_folder = os.path.join(tmp,f'{seednames[seed]}.tar.gz')
        #outputseed = os.path.join(tmp,f'{seednames[seed]}')
        # Extract the contents of the compressed archive
        #with tarfile.open(compressed_folder, 'r:gz') as archive:
        #    archive.extractall(outputseed)
        #os.remove(compressed_folder)
        pval_model_list=[]
        pval_age_list=[]
        pval_age2_list=[]
        pval_filepath = os.path.join(tmp, f'{seednames[seed]}','pvals')

        for voxel in range(len(idx_GM)):
            #print(f"Voxel {voxel} out of {len(idx_GM)} for seed {seed}: {seednames[seed]}")

            filename = f'voxel_{voxel}_{seednames[seed]}.txt'
            #If file does not exist. (i.e.: model did not converge) - replace with nans
            if not os.path.isfile(os.path.join(pval_filepath, filename)):
                pval_modelvsnull = np.nan
                pval_age = np.nan
                pval_age2 = np.nan
                pval_model_list.append(pval_modelvsnull)
                pval_age_list.append(pval_age)
                pval_age2_list.append(pval_age2)
            # Read the text file into a DataFrame
            else:
                df = pd.read_csv(os.path.join(pval_filepath, filename), sep='\t')

                # Access the p-values from the DataFrame
                pval_modelvsnull = df["pval_modelvsnull"].values[0]
                pval_age = df["pval_age"].values[0]
                pval_age2 = df["pval_age2"].values[0]


                pval_model_list.append(pval_modelvsnull)
                pval_age_list.append(pval_age)
                pval_age2_list.append(pval_age2)

        # Convert p-values to a NumPy array and invert (1-pval)
        pvalues_model = np.array(pval_model_list)
        inv_pvals_model = 1 - pvalues_model

        pvalues_age = np.array(pval_age_list)
        inv_pvals_age = 1 - pvalues_age

        pvalues_age2 = np.array(pval_age2_list)
        inv_pvals_age2 = 1 - pvalues_age2

        #Generate nifti files for SEED p-values using brainmask affine info
        #MODEL VS NULL PVAL
        seedmap_vol1 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol1 = seedmap_vol1.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol1[:] = 0 # clean img
        seedmap_vol1[0,idx_GM] = inv_pvals_model[:]
        seedmap_vol1 = seedmap_vol1.reshape(dim3d) #reshape to 3D

        # AGE PVAL
        seedmap_vol2 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol2 = seedmap_vol2.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol2[:] = 0 # clean img
        seedmap_vol2[0,idx_GM] = inv_pvals_age[:]
        seedmap_vol2 = seedmap_vol2.reshape(dim3d) #reshape to 3D

        # AGE2 PVAL
        seedmap_vol3 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol3 = seedmap_vol3.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol3[:] = 0 # clean img
        seedmap_vol3[0,idx_GM] = inv_pvals_age2[:]
        seedmap_vol3 = seedmap_vol3.reshape(dim3d) #reshape to 3D

        #Generate nifti objects
        seedmap_nii1 = nib.Nifti1Image(seedmap_vol1, Vgm_nii.affine)
        seedmap_nii2 = nib.Nifti1Image(seedmap_vol2, Vgm_nii.affine)
        seedmap_nii3 = nib.Nifti1Image(seedmap_vol3, Vgm_nii.affine)

        # Save as nifti files
        filepath = os.path.join(output,'results')

        filename1 = "developmentalchanges_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152NLin2009cAsym" + \
                    "_modelvsnull_1-pvals-uncorrected.nii.gz"
        filename2 = "developmentalchanges_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152NLin2009cAsym" + \
                    "_age_1-pvals-uncorrected.nii.gz"

        filename3 = "developmentalchanges_seed-" + \
                  seednames[seed] + \
                  "_space-MNI152NLin2009cAsym" + \
                  "_age2_1-pvals-uncorrected.nii.gz"

        if not os.path.exists(filepath):
            os.makedirs(filepath)

        seedmap_nii1.to_filename(os.path.join(filepath, filename1))
        seedmap_nii2.to_filename(os.path.join(filepath, filename2))
        seedmap_nii3.to_filename(os.path.join(filepath, filename3))

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








