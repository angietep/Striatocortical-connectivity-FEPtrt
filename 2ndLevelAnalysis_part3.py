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
        pval_filepath = os.path.join(tmp,'pvals', f'{seednames[seed]}')

        pval_interaction_list=[]
        pval_DIT_list=[]
        pval_Group_list=[]
        pval_APdose_list = []
        pval_PANSS_TP_list = []
        
        for voxel in range(len(idx_GM)):
            #print(f"Voxel {voxel} out of {len(idx_GM)} for seed {seed}: {seednames[seed]}")

            filename = f'voxel_{voxel}_{seednames[seed]}.txt'
            #If file does not exist. (i.e.: model did not converge) - replace with nans
            if not os.path.isfile(os.path.join(pval_filepath, filename)):
                pval_interaction = np.nan
                pval_DIT = np.nan
                pval_Group = np.nan
                pval_APdose = np.nan
                pval_PANSS_TP = np.nan
                
            # Read the text file into a DataFrame
            else:
                df = pd.read_csv(os.path.join(pval_filepath, filename), sep='\t')

                # Access the p-values from the DataFrame
                pval_interaction = df["pval_interaction"].values[0]
                pval_DIT = df["pval_DIT"].values[0]
                pval_Group = df["pval_Group"].values[0]
                pval_APdose = df["pval_APdose"].values[0]
                pval_PANSS_TP = df["pval_PANSSTP"].values[0]

            pval_interaction_list.append(pval_interaction)
            pval_DIT_list.append(pval_DIT)
            pval_Group_list.append(pval_Group)
            pval_APdose_list.append(pval_APdose)
            pval_PANSS_TP_list.append(pval_PANSS_TP)

        # Convert p-values to a NumPy array and invert (1-pval)
        pvalues_interaction = np.array(pval_interaction_list)
        inv_pvals_interaction = 1 - pvalues_interaction

        pvalues_DIT = np.array(pval_DIT_list)
        inv_pvals_DIT = 1 - pvalues_DIT

        pvalues_Group = np.array(pval_Group_list)
        inv_pvals_Group = 1 - pvalues_Group
        
        pvalues_APdose = np.array(pval_APdose_list) 
        inv_pvals_APdose = 1 - pvalues_APdose
        
        pvalues_PANSS_TP = np.array(pval_PANSS_TP_list)
        inv_pvals_PANSS_TP = 1 - pvalues_PANSS_TP

        #Generate nifti files for SEED p-values using brainmask affine info
        #INTERACTION PVAL
        seedmap_vol1 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol1 = seedmap_vol1.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol1[:] = 0 # clean img
        seedmap_vol1[0,idx_GM] = inv_pvals_interaction[:]
        seedmap_vol1 = seedmap_vol1.reshape(dim3d) #reshape to 3D

        # DIT PVAL
        seedmap_vol2 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol2 = seedmap_vol2.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol2[:] = 0 # clean img
        seedmap_vol2[0,idx_GM] = inv_pvals_DIT[:]
        seedmap_vol2 = seedmap_vol2.reshape(dim3d) #reshape to 3D

        # Group PVAL
        seedmap_vol3 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol3 = seedmap_vol3.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol3[:] = 0 # clean img
        seedmap_vol3[0,idx_GM] = inv_pvals_Group[:]
        seedmap_vol3 = seedmap_vol3.reshape(dim3d) #reshape to 3D

        # Group APdose
        seedmap_vol4 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol4 = seedmap_vol4.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol4[:] = 0 # clean img
        seedmap_vol4[0,idx_GM] = inv_pvals_APdose[:]
        seedmap_vol4 = seedmap_vol4.reshape(dim3d) #reshape to 3D

        # Group PANSS_TP
        seedmap_vol5 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol5 = seedmap_vol5.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol5[:] = 0 # clean img
        seedmap_vol5[0,idx_GM] = inv_pvals_PANSS_TP[:]
        seedmap_vol5 = seedmap_vol5.reshape(dim3d) #reshape to 3D


        #Generate nifti objects
        seedmap_nii1 = nib.Nifti1Image(seedmap_vol1, Vgm_nii.affine)
        seedmap_nii2 = nib.Nifti1Image(seedmap_vol2, Vgm_nii.affine)
        seedmap_nii3 = nib.Nifti1Image(seedmap_vol3, Vgm_nii.affine)
        seedmap_nii4 = nib.Nifti1Image(seedmap_vol4, Vgm_nii.affine)
        seedmap_nii5 = nib.Nifti1Image(seedmap_vol5, Vgm_nii.affine)

        # Save as nifti files

        filename1 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_interactionDITxgroup_1-pvals-uncorrected.nii.gz"
        
        filename2 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_DIT_1-pvals-uncorrected.nii.gz"

        filename3 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_Group_1-pvals-uncorrected.nii.gz"

        filename4 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_APdose_1-pvals-uncorrected.nii.gz"

        filename5 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_PANSS_TP_1-pvals-uncorrected.nii.gz"
                      

        if not os.path.exists(output):
            os.makedirs(output)

        seedmap_nii1.to_filename(os.path.join(output, filename1))
        seedmap_nii2.to_filename(os.path.join(output, filename2))
        seedmap_nii3.to_filename(os.path.join(output, filename3))
        seedmap_nii4.to_filename(os.path.join(output, filename4))
        seedmap_nii5.to_filename(os.path.join(output, filename5))

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








