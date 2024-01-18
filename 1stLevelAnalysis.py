
import numpy as np
import os, sys
import bids
import nibabel as nib
from scipy.stats import pearsonr
import statsmodels.api as sm
from nilearn import plotting as nplot
from nilearn import image as nimg
from nilearn.image import resample_to_img
import matplotlib.pyplot as plt
import argparse
import pandas as pd
#%%
##
def parse():

        options = argparse.ArgumentParser(description="Run 1st level analysis. Created by ...")
        options.add_argument('-p', '--participants', nargs='+',dest="participants", action='store', type=str, required=False,
                            help='id of subject or list of subjects')
        options.set_defaults(participants=None)
        options.add_argument('-w', '--workdir',dest="workdir", action='store', type=str, required=False,
                            help='the work directory for the project')
        options.set_defaults(workdir=os.environ["ROOTDIR"])
        options.add_argument('-b', '--bidsdir',dest="rawdata", action='store', type=str, required=False,
                            help='the work directory for the project')
        options.set_defaults(rawdata=os.path.join(os.environ["ROOTDIR"],"BIDS"))
        options.add_argument('-d', '--derivatives',dest="derivatives", action='store', type=str, required=True,
                            help='path to fMRIprep directory')
        options.add_argument('-o', '--outputs',dest="output", action='store', type=str, required=True,
                            help='path to fMRIprep directory')
        #print(options.parse_args())

        return options.parse_args()

#%%
def main ():
    os.environ["ROOTDIR"] = '/Users/brainsur/'  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir,"Desktop/striatconn")
        rawdata = os.path.join(workdir,"FEPtrt_prepro")
        #derivat = os.path.join(workdir,"derivatives","fmriprep-example")
        output  = os.path.join(workdir,"firstlevel")
        masks = os.path.join(workdir,"masks")
        participants = []

    else :
        options = parse()
        participants = options.participants
        workdir = options.workdir
        rawdata = options.rawdata
        derivat = options.derivatives
        output  = options.output
        script_path = os.path.dirname(__file__)
        masks = os.path.join(script_path,"masks")

    seednames = ['DCPutamen_space-MNI152_dim-9110991.nii.gz',
                 'DorsalCaudate_space-MNI152_dim-9110991.nii.gz',
                 'InfVentralCaudate_space-MNI152_dim-9110991.nii.gz',
                 'SupVentralCaudate_space-MNI152_dim-9110991.nii.gz',
                 'VRPutamen_space-MNI152_dim-9110991.nii.gz',
                 'DRPutamen_space-MNI152_dim-9110991.nii.gz']

    print('Rawdata: ', rawdata)

    #################
    Vgm_img = nib.load(os.path.join(masks,'GrayMattermask_thalamus_space-MNI152_dim-9110991.nii.gz'))
    #extract data array
    Vgm = Vgm_img.get_fdata()
    #save origianl dimensions (voxels_x, voxels_y, voxels_z)
    dim3d = Vgm.shape
    #reshape to 2D
    Vgm = Vgm.reshape(-1, np.prod(dim3d))  # -1 means auto-calculate size of dimension
    #save indexes in which Vgm == 1 (indexes for gray matter location)
    idx_GM = np.where(Vgm)[1]

    bidslayout = bids.BIDSLayout(rawdata,derivatives=False, validate = False) #With validate = True it doesn't find any subjects

    if not participants:
        participants = bidslayout.get_subjects()
    for i in  participants:
        print("Subject: ", i)
        i = i.replace("sub-", "")
        for ses in bidslayout.get_sessions(subject=i):
            print("Session: ", ses)
            #for multi-echo files, apply transform to standard space
            filesrest = bidslayout.get(subject=i,
                             session=ses,
                             extension=".nii.gz",
                             suffix="bold",
                             #regex_search=True,
                             #invalid_filters="allow"
                             )

            print("Filerest:", filesrest)
            if len(filesrest)>1:
                print(f"More than one file match bold filters. Subject {i}, ses {ses}")
                continue
            elif len(filesrest)<1:
                print(f"No file match bold filters. Subject {i}, ses {ses}")
                continue
            bold = filesrest[0]

            filesanat = bidslayout.get(subject=i,
                             session=ses,
                             extension=".nii.gz",
                             suffix="T1w",
                             #desc="preproc",
                             #regex_search=True,
                             invalid_filters="allow"
                             )

            print(f"Anat: {filesanat} ")
            if len(filesrest)>1:
                print(f"More than one file match anat filters. Subject {i}, ses {ses}")
                continue
            #anat = filesanat[0]

            #subj_i=participants.index(i)
            t = np.zeros((len(idx_GM), len(seednames)))

            #bold data
            func_img = nimg.load_img(bold)
            dim4d = func_img.shape
            N = dim4d[-1]  # number of time points

            # #reshape into 2d
            func_2d = func_img.get_fdata().reshape(-1, func_img.shape[-1]).T

            # pre-allocate space for seeds' time series
            TS_seeds = np.empty((N, len(seednames)))
            # Input new nifti for seeds
            # Compute SEEDs TS
            for k in range(len(seednames)):
                print('----- TS Seed:', seednames[k], '-----')

                seed = os.path.join(masks, seednames[k])
                #load seed image
                Vseed_img = nib.load(seed)
                Vseed = Vseed_img.get_fdata().reshape(func_2d.shape[-1])

                #for each time point, compute mean of voxels in which seed==1
                for vol in range(N):
                    Vaux = func_2d[vol, :]
                    TS_seeds[vol, k] = np.mean(Vaux[Vseed == 1])
                    del Vaux

                del Vseed

            # 1st Level Analysis: TSvox = b0 + b1*TSseed1 + b2*TSseed2 + ... + b6*TSseed6 + confounds
            print('First level analysis -----')
            print('TSvox = b0 + b1*TSseed1 + b2*TSseed2 + ... + b6*TSseed6')

            # Add a constant term to the independent variables
            design_matrix = TS_seeds
            design_matrix = sm.add_constant(design_matrix)

            for vox in range(len(idx_GM)):
                idx_fullvol = idx_GM[vox]
                # Fit the multiple linear regression model
                model = sm.OLS(func_2d[:, idx_fullvol], design_matrix)
                results = model.fit()
                # Get the t-values
                t_values = results.tvalues
                t[vox, :] = t_values[1:7]   # t[voxels,seeds]

            #SAVE tmaps as nifti files

            #Generate T-map using GMmask affine info
            tmap_img = Vgm_img
            tmap = tmap_img.get_fdata()

            for t_index in range(t.shape[1]):
                tmap = tmap.reshape(-1, np.prod(dim3d)) #reshape to 2D
                tmap[:] = 0 # clean img
                tmap[0,idx_GM] = t[:,t_index] #put t-values for each seed; seed 1 (t[:,0])
                tmap = tmap.reshape(dim3d) #reshape to 3D

                #Generate nifti object
                tmap_img = nib.Nifti1Image(tmap, tmap_img.affine)

                tmap_metadata = {
                    'task': 'rest',
                    'suffix': 'tstat',
                    'extension': '.nii.gz',
                    'session': str(ses),
                    'subject': 'sub-'+i,
                    'space': 'MNI152',
                    'seed': seednames[t_index].split("_")[0]
                }

                ##uncomment to plot tmap on top of T1 image
                #t1_img = nib.load(anat)
                #display=nplot.plot_anat(t1_img)
                #display.add_overlay(tmap_img)

                # Save each tmap as a nifti file
                filepath = os.path.join(output,tmap_metadata["subject"],"ses-"+tmap_metadata['session'],"func")
                filename = tmap_metadata["subject"] + "_" + \
                            "ses-"+tmap_metadata['session'] + "_" + \
                            "task-"+tmap_metadata['task'] + "_" + \
                            "space-"+tmap_metadata['space'] + "_" + \
                            "seed-"+tmap_metadata['seed'] + "_" + \
                            tmap_metadata["suffix"] + \
                            tmap_metadata["extension"]

                folder_subj=os.path.join(output,tmap_metadata["subject"])
                if not os.path.exists(folder_subj):
                    os.makedirs(folder_subj)

                if not os.path.exists(filepath):
                     os.makedirs(filepath)

                tmap_img.to_filename(os.path.join(filepath, filename))

            print('\t\t ----- Finished subject', i, ' Sess: ' ,ses,'----')

    #%%
if __name__ == '__main__':
    main()