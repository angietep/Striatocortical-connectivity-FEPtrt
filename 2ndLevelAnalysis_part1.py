# SECOND LEVEL PART 1: This script reads t-maps from all subjects (output of
# the first level) and saves in a tvals folder txt files with t-values of all
# subjects for each voxel of the gray matter mask.

# CONSIDERATIONS:

    # 1) The pyBIDS function to filter files based on 'invalid_filters' does not
    # work for me. Because of that, I cannot (easily) select the seed-file based
    # on file name. I'm exclusively relying on the order of the files being the
    # same for all subjects, i.e. seed=0 is DCPutamen because alphabetically it
    # is the same file, and it *should* be first for all subjects.
    # This is not the safest approach but it is fast. (I've checked it works for
    # a subset of 10 subjects)
    #
    # Example of the bids filter function which doesn't work:
    # bidslayout.get(subject=p,
                     # session=ses,
                     # run=r,
                     # extension=".nii.gz",
                     # suffix="tstat",
                     # space="MNI152NLin2009cAsym",
                     # #regex_search=True,
                     # # seed="DCPutamen", #does not work
                     # #invalid_filters="allow"
                     # )

#%%


import numpy as np
import os, sys
import bids
import nibabel as nib
from nilearn import image as nimg
import argparse
import pandas as pd
import shutil


#%%
def parse():

        options = argparse.ArgumentParser(description="Run 1st level analysis. Created by ...")
        options.add_argument('-p', '--participants', nargs='+',dest="participants", action='store', type=str, required=False,
                            help='id of subject or list of subjects')
        options.set_defaults(participants=None)
        options.add_argument('-w', '--workdir',dest="workdir", action='store', type=str, required=False,
                            help='the work directory for the project')
        options.set_defaults(workdir=os.environ["ROOTDIR"])
        #options.add_argument('-b', '--bidsdir',dest="rawdata", action='store', type=str, required=False,
        #                    help='the work directory for the project')
        #options.set_defaults(rawdata=os.path.join(os.environ["ROOTDIR"],"BIDS"))
        options.add_argument('-d', '--derivatives',dest="derivatives", action='store', type=str, required=True,
                            help='path to fMRIprep directory')
        options.add_argument('-o', '--outputs',dest="output", action='store', type=str, required=True,
                            help='path to fMRIprep directory')
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
        #workdir = os.path.join("/Volumes/TOSHIBA/striatconnTRT")
        rawdata = os.path.join(workdir,"FEPtrt_prepro") # sub-1401_ses-01_task-rest_fdpower.txt
        #rawdata = os.path.join(workdir,"../Preprocessed","FEPtrt_bids") # sub-1401_ses-01_task-rest_fdpower.txt
        masks = os.path.join(workdir,"masks")
        firstleveldir = os.path.join(workdir,"firstlevel") #  
        
        demographic = workdir
        #confounds = os.path.join(workdir,"derivatives","fmriprep")
        
        output  = os.path.join(workdir,"secondlevel")
        tmp  = os.path.join(output,"tvals")
        participants = []

    else :
        options = parse()
        participants = options.participants
        workdir = options.workdir
        #rawdata = options.rawdata
        #derivat = options.derivatives
        script_path = os.path.dirname(__file__)
        masks = os.path.join(script_path,"masks")

   
    seednames = ['DCPutamen',
                 'DorsalCaudate',
                 'DRPutamen',
                 'InfVentralCaudate',
                 'SupVentralCaudate',
                 'VRPutamen' #_space-MNI152NLin2009cAsym.nii.gz
                 ]

    print('firstlevel: ', firstleveldir)
      

    demographics=pd.read_csv(os.path.join(demographic,'subjects_finallist_DIT.csv'))
    demographics['MRI_'] = pd.to_datetime(demographics['MRI_'])
    demographics['Fecha_nacimiento'] = pd.to_datetime(demographics['Fecha_nacimiento'])
                   

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

    bidslayout = bids.BIDSLayout(firstleveldir, validate = False) #With validate = True it doesn't find any subjects
    confoundslayout = bids.BIDSLayout(rawdata, validate = False) #With validate = True it doesn't find any subjects

    if not participants:
        participants = bidslayout.get_subjects() #toma los del tsv y de algunos no tengo imágenes
        #participants_sub = ["sub-" + item for item in participants]

    for seed in range(len(seednames)):

        print(f"Seed {seednames[seed]}")
        tmap_allsubj = [] #initialize array for tmaps
        covars = []
        for p in participants:
            #print(f"Subject: {p}")
            p = p.replace("sub-", "")
            for ses in bidslayout.get_sessions(subject=p):
                
                if seed == 0: #only need to do this once
                    #Load demographics and covariates
                    dem_i = demographics[demographics.ID==float(p)] #find subject
                    dem_i = dem_i.reset_index()
                    
                    sex_i = dem_i.Sexo[float(ses)-1]
                    age_i = (dem_i.MRI_[float(ses)-1]-dem_i.Fecha_nacimiento[float(ses)-1]).days // 365 #find ses
                    group_i = dem_i.Grupo[float(ses)-1]
                    t_DIT_i = dem_i.t[float(ses)-1]
                    PANSSTP_i = dem_i.PANSS_TP_[float(ses)-1]
                    APdose_i = dem_i.APdose_[float(ses)-1]
                    
                    confounds_ents = {}
                    confounds_ents["subject"] = p
                    confounds_ents["session"] = ses
                    confounds_ents["task"] = 'rest'
                    #confounds_ents['desc'] = "confounds"
                    confounds_ents['suffix'] = "fdpower"
                    confounds_ents['extension'] = ".txt"
                    confounds = confoundslayout.get(return_type='file', **confounds_ents, invalid_filters="allow")[0]
                   
                    # sub-1401_ses-01_task-rest_fdpower.txt
                    noiseconfounds_df = pd.read_csv(confounds, sep='\t', header=None, names=["framewise_displacement"])
                    fdmean_i = noiseconfounds_df.framewise_displacement.mean()

                    covars.append((p, group_i, ses, t_DIT_i, PANSSTP_i, APdose_i, age_i, sex_i, fdmean_i))
                
                
                #print(f"Session: {ses}")
                tmap_sixseeds = bidslayout.get(subject=p,
                                 session=ses,
                                 extension=".nii.gz",
                                 suffix="tstat",
                                 space="MNI152",
                                 #regex_search=True,
                                 #seed="DCPutamen", #does not work
                                 #invalid_filters="allow"
                                 )
                if len(tmap_sixseeds) < 1:
                    continue
                #print(f"Subject: {p}, Session: {ses}, Seed: {tmap_sixseeds[seed].filename.split('seed-')[1].split('_')[0]}")

                #Load volume
                tmap_bids =tmap_sixseeds[seed]
                tmap_nii = nimg.load_img(tmap_bids)

                # reshape into 2d
                tmap_2d = tmap_nii.get_fdata().reshape(-1, np.prod(dim3d)).T
                tmap_gm = tmap_2d[idx_GM]

                tmap_allsubj.append(tmap_gm)
                #aux_list.append((p,ses,r))
             

        # Stack the arrays in the list horizontally
        tmap_allsubj = np.hstack(tmap_allsubj).T #dimensions: subj x voxels_GM

        outputseed = os.path.join(tmp,f'{seednames[seed]}')

        if not os.path.exists(output):
            os.mkdir(output)
        if not os.path.exists(tmp):
            os.mkdir(tmp)
        if not os.path.exists(outputseed):
            os.mkdir(outputseed)

        if seed == 0:
            df = pd.DataFrame(covars, columns=['ID', 'group','ses','t_DIT', 'PANSS_TP', 'APdose', 'age', 'sex', 'fdmean' ])
            df = df.dropna()
            df.to_csv(os.path.join(output,'df_covars.csv'), index=False)

        for voxel in range(len(idx_GM)):
            voxel_data = tmap_allsubj[:, voxel]  # Extract voxel ts (t-vals for each subject)
            file_name = f'voxel_{voxel}_{seednames[seed]}.txt'
            #You need to put these in different folders.
            if not os.path.exists(os.path.join(tmp,seednames[seed])):
                os.mkdir(os.path.join(tmp,seednames[seed]))
            file_path = os.path.join(tmp,seednames[seed], file_name)

           # Open the file in write mode ('w') to overwrite the existing content or create a new file
            with open(file_path, 'w') as file:
               for value in voxel_data:
                   file.write(f'{value:.6f}\n')    # Save voxel_data as a floating-point number with 6 decimal places

            #np.savetxt(os.path.join(tmp,file_name), voxel_data, fmt='%f')
            #print(f'Saved {file_name}')

        #print(f'Compressing folder {outputseed}')
        # Create the compressed archive
        #shutil.make_archive(outputseed, 'gztar',outputseed)
        # Remove the original folder
        #print(f'Folder compressed.\nRemoving original folder {outputseed}')
        #shutil.rmtree(outputseed)

#%%

if __name__ == '__main__':
    main()
