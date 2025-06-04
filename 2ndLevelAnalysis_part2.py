# SECOND LEVEL PART 2: This script reads txt files for each voxel, compute
# mixed model

  # TO RUN IN CLI (CALL one seed per terminal tab)
  # conda activate my-env 
  # cd ~/Desktop/striatconnTRT/secondlevel/tvals/InfVentralCaudate 
  # find . -maxdepth 1 -type f -name "*.txt" | xargs -n 1 -P 10 python ~/Desktop/GitHub_repos/Striatocortical-connectivity-FEPtrt/2ndLevelAnalysis_part2.py -c ~/Desktop/striatconnTRT/cleansample_covars.csv -o ~/Desktop/striatconnTRT/secondlevel/pvals -t ~/Desktop/striatconnTRT/secondlevel/tvals -v


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

    # 3) For voxels where not all subjects have values from the firstlevel (voxel_t has Nan)
    # I'm skipping the model and not saving a voxel.txt file - in 3rd part assing p = nan
    
  
#%%


import os, sys
import statsmodels.formula.api as smf
import argparse
import pandas as pd
import numpy as np

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
        output = os.path.join(workdir, "secondlevel","pvals")

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

    model_formula = 'voxel_t ~ age + sex + APdose + PANSS_TP + fdmean + HC + TRS + t_DIT + t_DIT:HC + t_DIT:TRS' 
    mixed_model = smf.mixedlm(model_formula, df, groups=df['ID'])
    result = mixed_model.fit()
    
    pval_HC = result.pvalues['HC']
    beta_HC = result.params['HC']
    
    pval_TRS = result.pvalues["TRS"]
    beta_TRS = result.params["TRS"]
    
    pval_time = result.pvalues['t_DIT']
    beta_time = result.params['t_DIT']

    pval_timexHC = result.pvalues['t_DIT:HC']
    beta_timexHC = result.params['t_DIT:HC']
   
    pval_timexTRS = result.pvalues['t_DIT:TRS']
    beta_timexTRS = result.params['t_DIT:TRS'] 
    
    pval_APdose = result.pvalues['APdose']
    beta_APdose = result.params['APdose']
    
    pval_PANSSTP = result.pvalues['PANSS_TP']
    beta_PANSSTP = result.params['PANSS_TP']

    # CONTRAST 1: Baseline TRS - HC
    param_names = result.params.index.tolist()
    contrast_baseline = np.zeros(len(param_names)-1)
    contrast_baseline[param_names.index('TRS')] = 1
    contrast_baseline[param_names.index('HC')] = -1
    contrast_result_baseline = result.t_test(np.atleast_2d(contrast_baseline))
    beta_contrast_baseline = contrast_result_baseline.effect.item()
    pval_contrast_baseline = contrast_result_baseline.pvalue.item()
    
    # CONTRAST 2: Slope (TRS*time - HC*time)
    contrast_slope = np.zeros(len(param_names)-1)
    contrast_slope[param_names.index('t_DIT:TRS')] = 1
    contrast_slope[param_names.index('t_DIT:HC')] = -1
    contrast_result_slope = result.t_test(np.atleast_2d(contrast_slope))
    beta_contrast_slope = contrast_result_slope.effect.item()
    pval_contrast_slope = contrast_result_slope.pvalue.item()
    
    # Clean up voxel column
    df = df.drop(columns=['voxel_t'])
    
    # Combine all outputs
    filevals = [
        beta_HC, pval_HC,
        beta_TRS, pval_TRS,
        beta_time, pval_time,
        beta_timexHC, pval_timexHC,
        beta_timexTRS, pval_timexTRS,
        beta_APdose, pval_APdose,
        beta_PANSSTP, pval_PANSSTP,
        beta_contrast_baseline, pval_contrast_baseline,
        beta_contrast_slope, pval_contrast_slope
    ]
    
    # Output path
    result_path = os.path.join(output, f"{seedname}", file_name)
    
    # Write results
    with open(os.path.join(result_path), 'w') as file:
        file.write(
            "beta_HC\tpval_HC\tbeta_TRS\tpval_TRS\tbeta_time\tpval_time\t"
            "beta_timexHC\tpval_timexHC\tbeta_timexTRS\tpval_timexTRS\t"
            "beta_APdose\tpval_APdose\tbeta_PANSSTP\tpval_PANSSTP\t"
            "beta_TRSvsHC\tpval_TRSvsHC\tbeta_timexTRSvsHC\tpval_timexTRSvsHC\n"
        )
        file.write("\t".join(f"{val}" for val in filevals) + "\n")
    
    print(f'Saved {result_path}')
    
    




#%%

if __name__ == '__main__':
    main()





