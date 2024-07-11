#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 14:55:48 2024

@author: brainsur
"""

# This script reads txt-files with p-values of each voxel,
# and generates nifti file with 1-pval on each voxel of the gm mask.

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

def main():
    os.environ["ROOTDIR"] = '/Users/brainsur/'  # seth path
    # os.environ["ROOTDIR"] = '/Users/angeles/'  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir, "Desktop/striatconnTRT")
        masks = os.path.join(workdir, "masks")
        tmp = os.path.join(workdir, "secondlevel")
        output = os.path.join(tmp, "results", "supp_uncorrected")
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
        ("APdose", "pval_APdose", "beta_APdose"), 
        ("PANSSTP", "pval_PANSSTP", "beta_PANSSTP")
    ]

    for seed in range(len(seednames)):
        print(f"\t ----- Working on seed {seed}: {seednames[seed]} ----")
        pval_filepath = os.path.join(tmp, 'pvals', f'{seednames[seed]}')

        pval_lists = {name: [] for name, _ , _ in variables}
        beta_lists = {name: [] for name, _ , _ in variables}

        for voxel in range(len(idx_GM)):
            filename = f'voxel_{voxel}_{seednames[seed]}.txt'
            if not os.path.isfile(os.path.join(pval_filepath, filename)):
                for name, _ , _ in variables:
                    pval_lists[name].append(np.nan)
                    beta_lists[name].append(0)
            else:
                df = pd.read_csv(os.path.join(pval_filepath, filename), sep='\t')
                for name, columnp , columnb in variables:
                    pval_lists[name].append(df[columnp].values[0])
                    if df[columnp].values[0] < 0.05:
                        beta_lists[name].append(df[columnb].values[0])
                    else: beta_lists[name].append(0)

        #inv_pvals = {name: 1 - np.array(pval_list) for name, pval_list in pval_lists.items()}

        for i, (name, _ , _) in enumerate(variables):
            seedmap_vol = Vgm_vol.reshape(-1, np.prod(dim3d))
            seedmap_vol[:] = 0
            seedmap_vol[0, idx_GM] = np.array(beta_lists[name][:])
            seedmap_vol = seedmap_vol.reshape(dim3d)

            seedmap_nii = nib.Nifti1Image(seedmap_vol, Vgm_nii.affine)
            filename = f"longitudinalTRT_seed-{seednames[seed]}_space-MNI152_dim-9110991_{name}_betas-for-uncorrected-pvals-under05.nii.gz"

            if not os.path.exists(output):
                os.makedirs(output)

            seedmap_nii.to_filename(os.path.join(output, filename))

        print(f"\t ----- Finished seed {seed}: {seednames[seed]} ----")


#%%

if __name__ == '__main__':
    main()








