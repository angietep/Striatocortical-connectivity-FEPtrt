 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:24:09 2024

@author: brainsur
"""

# This script reads txt files for each voxel and compute
# mixed model

  # TO RUN IN CLI (CALL one seed per terminal tab)
  # conda activate my-env 
  # cd ~/Desktop/striatconnTRT/secondlevel/tvals/InfVentralCaudate 
  # find . -maxdepth 1 -type f -name "*.txt" | xargs -n 1 -P 10 python ~/Desktop/GitHub_repos/Striatocortical-connectivity-FEPtrt/suppFig1.py -c ~/Desktop/striatconnTRT/cleansample_covars.csv -o ~/Desktop/striatconnTRT/secondlevel/pvals_suppFig1 -t ~/Desktop/striatconnTRT/secondlevel/tvals -v

    
  
#%%


import os, sys
import statsmodels.formula.api as smf
import argparse
import pandas as pd

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning

# Suppress convergence warnings
warnings.filterwarnings('ignore', category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="statsmodels")


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
        #print(options.parse_args())

        return options.parse_args()

#%%
def main ():
    os.environ["ROOTDIR"] = '/Users/brainsur/Desktop/'  # seth path
    # os.environ["ROOTDIR"] = '/Users/angeles/'  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir, "striatconnTRT")
        tvalpath = os.path.join(workdir, "secondlevel", "tvals")
        covariates = os.path.join(workdir, 'cleansample_covars.csv')
        output = os.path.join(workdir, "secondlevel","pvals_suppFig1")

        voxel = "voxel_33874_DorsalCaudate.txt"
        #seedname = 'DorsalCaudate'
       

    else:
        options = parse()
        voxel = options.voxel
        covariates = options.covariates
        tvalpath = options.tvalpath
        output = options.output

      
    df = pd.read_csv(covariates, header=0)

   
    seedname = voxel.split("_")[-1]
    seedname = seedname.split(".")[0] 
   
    if not os.path.exists(os.path.join(output,f"{seedname}")):
        os.makedirs(os.path.join(output, f"{seedname}"))       
    
    
    file_name = voxel
    file_path = os.path.join(
        tvalpath, f"{seedname}", file_name)

    voxel_t = pd.read_csv(
        file_path, header=None, names=['voxel_t'])
                
    if voxel_t['voxel_t'].isnull().any() | voxel_t['voxel_t'].isnull().all():
        print(f"Skipping voxel {voxel} due to NaN values in 'voxel_t'")
        sys.exit()
        
    df['voxel_t'] = voxel_t

    model_formula = 'voxel_t ~ age + sex ' 
    mixed_model = smf.mixedlm(model_formula, df, groups=df['ID'])
    result = mixed_model.fit()
    
    pval = result.pvalues['Intercept']
    beta = result.params['Intercept']
   
    df = df.drop(columns=['voxel_t'])

    filevals = [beta, pval]

    result_path = os.path.join(output, f"{seedname}", file_name)

    with open(os.path.join(result_path), 'w') as file:
        # Write column headers
        file.write(
            "beta_intercept\tpval_intercept\n")
        # Write p-values in different columns
        for val in filevals:
            if val == filevals[-1]:
                file.write(f"{val}")
            else:
                file.write(f"{val}\t")  # Write p-values
    print(f'Saved {result_path}')



#%%

if __name__ == '__main__':
    main()







