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
    os.environ["ROOTDIR"] = '/Users/brainsur/Desktop'  # seth path
    #os.environ["ROOTDIR"] = '/Volumes/TOSHIBA'  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir, "striatconnTRT")
        masks = os.path.join(workdir, "masks")
        firstleveldir = os.path.join(workdir, "firstlevel")
        secondlevel = os.path.join(workdir, "secondlevel")
        participants = []

    else:
        options = parse()
        participants = options.participants
        workdir = options.workdir
        secondlevel = options.secondlevel

        script_path = os.path.dirname(__file__)

    print('firstlevel: ', firstleveldir)
    sample = pd.read_csv(os.path.join(workdir,'cleansample_covars.csv'))
    intermediate_output = os.path.join(secondlevel,'clustcorrection_acfparameters_p01.csv')
    
    
    if not(os.path.exists(intermediate_output)):
        print("Working on first level data to generate dataframe with 3dFWHMx results- RUN AGAIN WHEN FINISHED")
        #################
        GMmask = os.path.join(masks,'GrayMattermask_thalamus_space-MNI152_dim-9110991.nii.gz')
 
        bidslayout = bids.BIDSLayout(firstleveldir, validate=False)
    
        if not participants:
            participants = bidslayout.get_subjects()
    
        # Define the output DataFrame
        output_df = pd.DataFrame(
            columns=['subject', 'ses', 'acf_x', 'acf_y', 'acf_z'])
    
        for p in participants:
            if p in sample['ID'].values:
                #print(f"Subject: {p}")
                p = p.replace("sub-", "")
                for ses in bidslayout.get_sessions(subject=p):
                    if int(ses) in sample.ses[sample['ID']==p].values:
                        #print(f"Session: {ses}")
                    
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
            
                        #afni_command = ['3dFWHMx',
                        afni_command = ['/Users/brainsur/abin/3dFWHMx',
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
       
        output_df.to_csv(intermediate_output, index=False)
        
        
    
    else:
        print("CSV with 3dFWHMx data found - not reading first level data")
        output_df = pd.read_csv(intermediate_output)
        
        acf_x = output_df['acf_x'].mean()
        acf_y = output_df['acf_y'].mean()
        acf_z = output_df['acf_z'].mean()
    
        #afni_command = ['3dClustSim',
        afni_command = ['/Users/brainsur/abin/3dClustSim',
                        '-acf', str(acf_x), str(acf_y), str(acf_z), 
                        '-nxyz', '91','109','91', 
                        '-dxyz', '2','2','2',
                        '-athr', '0.05',#'0.05/6', #Bonferroni corrected for 6 seeds
                        '-pthr', '0.01'] #0.001
    
        try:
            result = subprocess.run(
                afni_command, check=True, capture_output=True, text=True)
            
            print(result.stdout)
            with open(os.path.join(secondlevel,'3dClustSim_output_p01.txt'), 'w') as file:
                file.write(result.stdout)
    
           
        except subprocess.CalledProcessError as e:
            # Print the error if the command fails
            print(f"Error: {e}")
            print(f"Command output: {e.output}")
    
    
                

#%%

if __name__ == '__main__':
    main()
