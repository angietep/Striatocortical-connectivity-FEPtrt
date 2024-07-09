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

#%% MAIN old (no for loop )
def main_no ():
    os.environ["ROOTDIR"] = '/Users/brainsur/'  # seth path
    #os.environ["ROOTDIR"] = '/Users/angeles/'  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir,"Desktop/striatconnTRT")
        masks = os.path.join(workdir,"masks")
        tmp = os.path.join(workdir, "secondlevel")
        output = os.path.join(tmp, "results","uncorrected")
        
        
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

        pval_HC_list=[]
        pval_TRS_list=[]
        pval_time_list=[]
        pval_timexHC_list = []
        pval_timexTRS_list = []
        pval_APdose_list = []
        pval_PANSSTP_list = []
        
        for voxel in range(len(idx_GM)):
            #print(f"Voxel {voxel} out of {len(idx_GM)} for seed {seed}: {seednames[seed]}")

            filename = f'voxel_{voxel}_{seednames[seed]}.txt'
            #If file does not exist. (i.e.: model did not converge) - replace with nans
            if not os.path.isfile(os.path.join(pval_filepath, filename)):
                pval_HC = np.nan
                pval_TRS = np.nan
                pval_time = np.nan
                pval_timexHC = np.nan
                pval_timexTRS = np.nan
                pval_APdose = np.nan
                pval_PANSSTP = np.nan
                
            # Read the text file into a DataFrame
            else:
                df = pd.read_csv(os.path.join(pval_filepath, filename), sep='\t')
                
                # Access the p-values from the DataFrame
                pval_HC = df["pval_HC"].values[0]
                pval_TRS = df["pval_TRS"].values[0]
                pval_time = df["pval_time"].values[0]
                pval_timexHC = df["pval_timexHC"].values[0]
                pval_timexTRS = df["pval_timexTRS"].values[0]
                pval_APdose = df["pval_APdose"].values[0]
                pval_PANSSTP = df["pval_PANSSTP"].values[0]

            pval_HC_list.append(pval_HC)
            pval_TRS_list.append(pval_TRS)
            pval_time_list.append(pval_time)
            pval_timexHC_list.append(pval_timexHC)
            pval_timexTRS_list.append(pval_timexTRS)
            pval_APdose_list.append(pval_APdose)
            pval_PANSSTP_list.append(pval_PANSSTP)

        # Convert p-values to a NumPy array and invert (1-pval)
        pvalues_HC = np.array(pval_HC_list)
        inv_pvals_HC = 1 - pvalues_HC

        pvalues_TRS = np.array(pval_TRS_list)
        inv_pvals_TRS = 1 - pvalues_TRS

        pvalues_time = np.array(pval_time_list)
        inv_pvals_time = 1 - pvalues_time
        
        pvalues_timexHC = np.array(pval_timexHC_list) 
        inv_pvals_timexHC = 1 - pvalues_timexHC
        
        pvalues_timexTRS = np.array(pval_timexTRS_list)
        inv_pvals_timexTRS = 1 - pvalues_timexTRS
        
        pvalues_APdose = np.array(pval_APdose_list)
        inv_pvals_APdose = 1 - pvalues_APdose
        
        pvalues_PANSSTP = np.array(pval_PANSSTP_list)
        inv_pvals_PANSSTP = 1 - pvalues_PANSSTP

        #Generate nifti files for SEED p-values using brainmask affine info
        #HC PVAL
        seedmap_vol1 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol1 = seedmap_vol1.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol1[:] = 0 # clean img
        seedmap_vol1[0,idx_GM] = inv_pvals_HC[:]
        seedmap_vol1 = seedmap_vol1.reshape(dim3d) #reshape to 3D

        # TRS PVAL
        seedmap_vol2 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol2 = seedmap_vol2.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol2[:] = 0 # clean img
        seedmap_vol2[0,idx_GM] = inv_pvals_TRS[:]
        seedmap_vol2 = seedmap_vol2.reshape(dim3d) #reshape to 3D

        # Group time
        seedmap_vol3 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol3 = seedmap_vol3.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol3[:] = 0 # clean img
        seedmap_vol3[0,idx_GM] = inv_pvals_time[:]
        seedmap_vol3 = seedmap_vol3.reshape(dim3d) #reshape to 3D

        # Group timexHC
        seedmap_vol4 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol4 = seedmap_vol4.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol4[:] = 0 # clean img
        seedmap_vol4[0,idx_GM] = inv_pvals_timexHC[:]
        seedmap_vol4 = seedmap_vol4.reshape(dim3d) #reshape to 3D

        # Group timexTRS
        seedmap_vol5 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol5 = seedmap_vol5.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol5[:] = 0 # clean img
        seedmap_vol5[0,idx_GM] = inv_pvals_timexTRS[:]
        seedmap_vol5 = seedmap_vol5.reshape(dim3d) #reshape to 3D
        
        #APdose PVAL
        seedmap_vol6 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol6 = seedmap_vol6.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol6[:] = 0 # clean img
        seedmap_vol6[0,idx_GM] = inv_pvals_APdose[:]
        seedmap_vol6 = seedmap_vol6.reshape(dim3d) #reshape to 3D
        
        #PANSSTP PVAL
        seedmap_vol7 = Vgm_vol #Vgm_vol = Vgm_nii.get_fdata()
        seedmap_vol7 = seedmap_vol7.reshape(-1, np.prod(dim3d)) #reshape to 2D
        seedmap_vol7[:] = 0 # clean img
        seedmap_vol7[0,idx_GM] = inv_pvals_PANSSTP[:]
        seedmap_vol7 = seedmap_vol7.reshape(dim3d) #reshape to 3D


        #Generate nifti objects
        seedmap_nii1 = nib.Nifti1Image(seedmap_vol1, Vgm_nii.affine)
        seedmap_nii2 = nib.Nifti1Image(seedmap_vol2, Vgm_nii.affine)
        seedmap_nii3 = nib.Nifti1Image(seedmap_vol3, Vgm_nii.affine)
        seedmap_nii4 = nib.Nifti1Image(seedmap_vol4, Vgm_nii.affine)
        seedmap_nii5 = nib.Nifti1Image(seedmap_vol5, Vgm_nii.affine)
        seedmap_nii6 = nib.Nifti1Image(seedmap_vol6, Vgm_nii.affine)
        seedmap_nii7 = nib.Nifti1Image(seedmap_vol7, Vgm_nii.affine)

        # Save as nifti files

        filename1 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_HC_1-pvals-uncorrected.nii.gz"
        
        filename2 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_TRS_1-pvals-uncorrected.nii.gz"

        filename3 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_time_1-pvals-uncorrected.nii.gz"

        filename4 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_timexHC_1-pvals-uncorrected.nii.gz"

        filename5 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_timexTRS_1-pvals-uncorrected.nii.gz"
                    
        filename6 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_APdose_1-pvals-uncorrected.nii.gz"
                    
        filename7 = "longitudinalTRT_seed-" + \
                    seednames[seed] + \
                    "_space-MNI152_dim-9110991" + \
                    "_PANSSTP_1-pvals-uncorrected.nii.gz"
                      

        if not os.path.exists(output):
            os.makedirs(output)

        seedmap_nii1.to_filename(os.path.join(output, filename1))
        seedmap_nii2.to_filename(os.path.join(output, filename2))
        seedmap_nii3.to_filename(os.path.join(output, filename3))
        seedmap_nii4.to_filename(os.path.join(output, filename4))
        seedmap_nii5.to_filename(os.path.join(output, filename5))
        seedmap_nii6.to_filename(os.path.join(output, filename6))
        seedmap_nii7.to_filename(os.path.join(output, filename7))

        print(f"\t ----- Finished seed {seed}: {seednames[seed]} ----")

        # # Create the compressed archive
        # print(f'Compressing folder {outputseed}')
        # shutil.make_archive(outputseed, 'gztar',outputseed)
        # # Remove the original folder
        # print(f'Folder compressed.\nRemoving original folder {outputseed}')
        # shutil.rmtree(outputseed)

#%% MAIN 
def main():
    os.environ["ROOTDIR"] = '/Users/brainsur/'  # seth path
    # os.environ["ROOTDIR"] = '/Users/angeles/'  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir, "Desktop/striatconnTRT")
        masks = os.path.join(workdir, "masks")
        tmp = os.path.join(workdir, "secondlevel")
        output = os.path.join(tmp, "results", "uncorrected")
    else:
        options = parse()
        tmp = options.workdir
        masks = options.mask
        output = options.outputs

    seednames = [
        'DCPutamen',
        'DorsalCaudate',
        'DRPutamen',
        'InfVentralCaudate',
        'SupVentralCaudate',
        'VRPutamen' # _space-MNI152NLin2009cAsym.nii.gz
    ]

    Vgm_nii = nib.load(os.path.join(masks, 'GrayMattermask_thalamus_space-MNI152_dim-9110991.nii.gz'))
    Vgm_vol = Vgm_nii.get_fdata()
    dim3d = Vgm_vol.shape
    Vgm_2d = Vgm_vol.reshape(-1, np.prod(dim3d)).T
    idx_GM = np.where(Vgm_2d)[0]

    variables = [
        ("HC", "pval_HC"), 
        ("TRS", "pval_TRS"), 
        ("time", "pval_time"), 
        ("timexHC", "pval_timexHC"), 
        ("timexTRS", "pval_timexTRS"), 
        ("APdose", "pval_APdose"), 
        ("PANSSTP", "pval_PANSSTP")
    ]

    for seed in range(len(seednames)):
        print(f"\t ----- Working on seed {seed}: {seednames[seed]} ----")
        pval_filepath = os.path.join(tmp, 'pvals', f'{seednames[seed]}')

        pval_lists = {name: [] for name, _ in variables}

        for voxel in range(len(idx_GM)):
            filename = f'voxel_{voxel}_{seednames[seed]}.txt'
            if not os.path.isfile(os.path.join(pval_filepath, filename)):
                for name, _ in variables:
                    pval_lists[name].append(np.nan)
            else:
                df = pd.read_csv(os.path.join(pval_filepath, filename), sep='\t')
                for name, column in variables:
                    pval_lists[name].append(df[column].values[0])

        inv_pvals = {name: 1 - np.array(pval_list) for name, pval_list in pval_lists.items()}

        for i, (name, _) in enumerate(variables):
            seedmap_vol = Vgm_vol.reshape(-1, np.prod(dim3d))
            seedmap_vol[:] = 0
            seedmap_vol[0, idx_GM] = inv_pvals[name][:]
            seedmap_vol = seedmap_vol.reshape(dim3d)

            seedmap_nii = nib.Nifti1Image(seedmap_vol, Vgm_nii.affine)
            filename = f"longitudinalTRT_seed-{seednames[seed]}_space-MNI152_dim-9110991_{name}_1-pvals-uncorrected.nii.gz"

            if not os.path.exists(output):
                os.makedirs(output)

            seedmap_nii.to_filename(os.path.join(output, filename))

        print(f"\t ----- Finished seed {seed}: {seednames[seed]} ----")

if __name__ == "__main__":
    main()


#%%

if __name__ == '__main__':
    main()








