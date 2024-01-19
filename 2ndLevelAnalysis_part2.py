# SECOND LEVEL PART 2: This script reads txt files for each voxel, compute
# mixed models, compare them and save a new txt-file with the p-value of such
# comparison (chi-squared of the likelihood ratio test)

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

    # 2) I'm ignoring convergence warnings from MixedLM models to avoid
    # having as many warning messages as voxels for which we are fitting models.
    # I think it's okay to do it in this context, since we don't care about the
    # actual output of the models, but only about which model fits best.


    # 3) For voxels where not all subjects have values from the firstlevel (voxel_t has Nan)
    # I'm skipping the model and not saving a voxel.txt file - in 3rd part assing p = nan
    
    # TO RUN IN CLI
    #  conda activate my-env 
    # cd ~/Desktop/striatconn/secondlevel/tmp/InfVentralCaudate 
    # find . -maxdepth 1 -type f -name "*.txt" | xargs -n 1 -P 10 python ~/Desktop/code/Striatocortical-connectivity-FEPtrt/2ndLevelAnalysis_part2.py -c ~/Desktop/striatconn/secondlevel/df_covars.csv -o ~/Desktop/striatconn/secondlevel/ -t ~/Desktop/striatconn/secondlevel/tmp -v

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
    os.environ["ROOTDIR"] = '/Users/brainsur/'  # seth path
    # os.environ["ROOTDIR"] = '/Users/angeles/'  # seth path
    rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        workdir = os.path.join(rootdir, "Desktop/striatconn")
        tvalpath = os.path.join(workdir, "secondlevel", "tmp")
        covariates = os.path.join(workdir, "secondlevel", "df_covars.csv")
        output = os.path.join(workdir, "secondlevel")

        voxel = "voxel_33874_DorsalCaudate.txt"
        #seedname = 'DorsalCaudate'
        #'DCPutamen'          
         #
          #    'DRPutamen',
           #   'InfVentralCaudate',
            #  'SupVentralCaudate',
             # 'VRPutamen' 

    else:
        options = parse()
        voxel = options.voxel
        covariates = options.covariates
        tvalpath = options.tvalpath
        output = options.output

      
    df = pd.read_csv(covariates, header=0)
   
    seedname = voxel.split("_")[-1]
    seedname = seedname.split(".")[0] 
   
    if not os.path.exists(os.path.join(output,f"{seedname}" ,'pvals')):
        os.makedirs(os.path.join(output, f"{seedname}",'pvals'))       
    
            
    print(f"{voxel} out of 137035")
    
    file_name = voxel
    file_path = os.path.join(
        tvalpath, f"{seedname}", file_name)

    voxel_t = pd.read_csv(
        file_path, header=None, names=['voxel_t'])
                
    if voxel_t['voxel_t'].isnull().any() | voxel_t['voxel_t'].isnull().all():
        print(f"Skipping voxel {voxel} due to NaN values in 'voxel_t'")
        sys.exit()
        
    df['voxel_t'] = voxel_t

    model_formula = 'voxel_t ~ age + sex + fdmean + t_DIT+ group + t_DIT:group + (1|ID)'
    mixed_model = smf.mixedlm(model_formula, df, groups=df['ID'])
    result = mixed_model.fit()

    pval_interaction = result.pvalues['t_DIT:group[T.nonTRT]']
    beta_interaction = result.params['t_DIT:group[T.nonTRT]']
    pval_DIT = result.pvalues['t_DIT']
    beta_DIT = result.params['t_DIT']
    pval_Group = result.pvalues['group[T.nonTRT]']
    beta_Group = result.params['group[T.nonTRT]']

    df = df.drop(columns=['voxel_t'])

    filevals = [beta_interaction, pval_interaction, beta_DIT, pval_DIT, beta_Group, pval_Group]

    result_path = os.path.join(output, f"{seedname}",'pvals', file_name)

    with open(os.path.join(result_path), 'w') as file:
        # Write column headers
        file.write(
            "beta_interaction\tpval_interaction\tbeta_DIT\tpval_DIT\tbeta_Group\tpval_Group\n")
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





