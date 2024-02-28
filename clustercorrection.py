#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:56:57 2024

@author: angeles
"""

# CLUSTER CORRECTION WITH AFNI
# 1) Read residuals from firstlevel analysis for each subject
# 2) Compute 3dFWHMx to get acf parameters
# 3) save output to csv
# 4) take mean of acf parameters 
# 5) call 3dclustsim to obtain size for clusters correction


#%%

import os, sys
import bids
import argparse
import pandas as pd

import subprocess



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
    #os.environ["ROOTDIR"] = '/Users/brainsur/Desktop'  # seth path
    os.environ["ROOTDIR"] = '/Volumes/TOSHIBA'  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir, "striatconnTRT")
        masks = os.path.join(workdir, "masks")
        firstleveldir = os.path.join(workdir, "firstlevel")

        output = workdir
        participants = []

    else:
        options = parse()
        participants = options.participants
        workdir = options.workdir
        # rawdata = options.rawdata
        # derivat = options.derivatives
        script_path = os.path.dirname(__file__)

    print('firstlevel: ', firstleveldir)

    #################
    GMmask = os.path.join(masks,'GrayMattermask_thalamus_space-MNI152_dim-9110991.nii.gz')
    # read vol
    # Vgm_vol = Vgm_nii.get_fdata()
    # save origianl dimensions (voxels_x, voxels_y, voxels_z)
    # dim3d = Vgm_vol.shape
    # reshape to 2D
    # Vgm_2d = Vgm_vol.reshape(-1, np.prod(dim3d)).T  # -1 means auto-calculate size of dimension
    # save indexes in which Vgm == 1 (indexes for gray matter location)
    # idx_GM = np.where(Vgm_2d)[0]

    # With validate = True it doesn't find any subjects
    bidslayout = bids.BIDSLayout(firstleveldir, validate=False)
    # confoundslayout = bids.BIDSLayout(rawdata, validate = False) #With validate = True it doesn't find any subjects

    if not participants:
        # toma los del tsv y de algunos no tengo imágenes
        participants = bidslayout.get_subjects()
        # participants_sub = ["sub-" + item for item in participants]

    # Define the output DataFrame
    output_df = pd.DataFrame(
        columns=['subject', 'ses', 'acf_x', 'acf_y', 'acf_z'])

    for p in participants:
        # print(f"Subject: {p}")
        p = p.replace("sub-", "")
        for ses in bidslayout.get_sessions(subject=p):
            # print(f"Session: {ses}")
            residuals = bidslayout.get(subject=p,
                                       session=ses,
                                       extension=".nii.gz",
                                       suffix="residuals",
                                       space="MNI152",
                                       # regex_search=True,
                                       # seed="DCPutamen", #does not work
                                       # invalid_filters="allow"
                                       )
            if len(residuals) < 1:
                continue
            print(f"Subject: {p}, Session: {ses}, Residuals: {residuals}")

            # Replace 'your_command_here' with the actual 3dclustsim command and its arguments
            #afni_command = ['/Users/brainsur/abin/3dFWHMx',
            afni_command = ['3dFWHMx',
                            '-mask', GMmask,
                            '-input', residuals[0].path]

            try:
                result = subprocess.run(
                    afni_command, check=True, capture_output=True, text=True)
                # Print the output if needed
                # Extract the three numbers from the result.stdout
                values = [float(val) for val in result.stdout.split()[4:7]]

                # Append the results to the DataFrame
                # output_df = output_df.append({'subject': p,'ses': ses, 'acf_x': values[0], 'acf_y': values[1], 'acf_z': values[2]}, ignore_index=True)
                output_df = pd.concat([output_df, pd.DataFrame({'subject': [p], 'ses': [ses], 'acf_x': [
                                      values[0]], 'acf_y': [values[1]], 'acf_z': [values[2]]})], ignore_index=True)

            except subprocess.CalledProcessError as e:
                # Print the error if the command fails
                print(f"Error: {e}")
                print(f"Command output: {e.output}")

    # Save the DataFrame to a CSV file
    output_csv_path = os.path.join(output,'clustcorrection_acfparameters.csv')
   
    output_df.to_csv(output_csv_path, index=False)
     
    acf_x = output_df['acf_x'].mean()
    acf_y = output_df['acf_y'].mean()
    acf_z = output_df['acf_z'].mean()
    
    # Replace 'your_command_here' with the actual 3dclustsim command and its arguments
    #afni_command = ['/Users/brainsur/abin/3dClustSim',
    afni_command = ['3dClustSim',
                    '-acf', str(acf_x), str(acf_y), str(acf_z), 
                    '-nxyz', '91','109','91', 
                    '-dxyz', '2','2','2',
                    '-athr', '0.05',
                    '-pthr', '0.001']

    try:
        result = subprocess.run(
            afni_command, check=True, capture_output=True, text=True)
        
        print(result.stdout)
        # Extract the three numbers from the result.stdout
       
    except subprocess.CalledProcessError as e:
        # Print the error if the command fails
        print(f"Error: {e}")
        print(f"Command output: {e.output}")


            
#            3dClustSim -acf 0.38942 4.91292 11.9489 -nxyz 91 109 91 -dxyz 2 2 2 -athr 0.05 -pthr 0.001
            

#%%

if __name__ == '__main__':
    main()
