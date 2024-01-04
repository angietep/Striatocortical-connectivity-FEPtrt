# SECOND LEVEL PART 2: This script reads txt files for each voxel, compute
# mixed models, compare them and save a new txt-file with the p-value of such
# comparison (chi-squared of the likelihood ratio test)

# CONSIDERATIONS:

    # TSV files with preprocessing confounds need to be in a different folder
    # than firstlevel outputs. If not, pyBIDS will read subjects list from
    # TSV folder (and there are some subjects for who we have tsv folder but
    # no firstlevel output)

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

    # 3) We're comparing a quadratic model (age + age2) to an intercept-only model.
    # I think it's okay but the most typical thing to do is to go step-by-step, i.e.,
    # comparing quadratic to linear and linear to intercept. Also, if the
    # 'grountruth' were to be a linear association, we might be missing some
    # power by adding the quadratic term, aren't we? Anyhow, we should
    # explain why we choose a quadratic approach...

    # 4) There are some subjects with an apparent mismatch of information. For example,
    # subject ID 30 has fmri data for sessions 1 and 3, but demographic data for
    # sessions 1 and 2. As the code is today, we are losing information from fmri
    # session 3 due to this.

    # 5) For some voxels I got an error LinAlgError from the models (singular matrix).
    # This are the voxels for which not all subjects had values from the firstlevel,
    # due either to FOV issues or misalignment. I'm assigning p-val = nan in those voxels.

    # 6) I don't know why there are some subjects who have more than one run for the same
    # session. We are using all runs available.
#%%

import numpy as np
import os, sys
#import bids
#import nibabel as nib
#from scipy.stats import pearsonr
from scipy.stats import chi2
import statsmodels.api as sm
import statsmodels.formula.api as smf
#from nilearn import plotting as nplot
#from nilearn import image as nimg
#from nilearn.image import resample_to_img
#import matplotlib.pyplot as plt
import subprocess
import argparse
import pandas as pd

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from scipy.linalg import LinAlgError


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

        options.add_argument('-t', '--tarpath',dest="tarpath", action='store', type=str, required=False,
                             help='the work directory for the tarfile')
        # options.set_defaults(workdir=os.environ["ROOTDIR"])
        # options.add_argument('-b', '--bidsdir',dest="rawdata", action='store', type=str, required=False,
        #                     help='the work directory for the project')
        # options.set_defaults(rawdata=os.path.join(os.environ["ROOTDIR"],"BIDS"))
        # options.add_argument('-d', '--derivatives',dest="derivatives", action='store', type=str, required=True,
        #                     help='path to fMRIprep directory')
        options.add_argument('-o', '--output',dest="output", action='store', type=str, required=True,
                            help='path to output directory')
        #print(options.parse_args())

        return options.parse_args()

#%%
def main ():
    #os.environ["ROOTDIR"] = '/Users/angeles/'  # seth path
    #rootdir = os.environ["ROOTDIR"]
    if hasattr(sys, "ps1"):
        options = {}
        #workdir = os.environ["ROOTDIR"]
        voxel = "./voxel_128025_DCPutamen.txt"
        tarpath = "/data/NIMH_scratch/zugmana2/angeles/INPD_secondlevel-2/tmp/DCPutamen.tar.gz"
        covariates = "/data/zugmana2/INPD/derivatives/angeles/INPD_secondlevel-2/tmp/df_covars.csv"
        output = "/data/zugmana2/INPD/derivatives/angeles/INPD_secondlevel-2/test1"
        #rawdata = os.path.join(workdir,"fmriprep_examples")
        #masks = os.path.join(rawdata,"masks")
        #firstleveldir = os.path.join(workdir,"INPD_striatal-connectivity") #  "INPD_subsample") #
#        confounds = os.path.join(workdir,"INPD_tsvs")
#        demographic = firstleveldir
        #output  = os.path.join(firstleveldir,"../INPD_secondlevel")
        #tmp  = os.path.join(output,"tmp")
#        participants = []

    else :
        options = parse()
        voxel = options.voxel
        covariates = options.covariates
#        participants = options.participants
        # workdir = options.workdir
        # rawdata = options.rawdata
 #       derivat = options.derivatives
        tarpath = options.tarpath
        output  = options.output


    # seednames = ['DCPutamen',
    #              'DorsalCaudate',
    #              'DRPutamen',
    #              'InfVentralCaudate',
    #              'SupVentralCaudate',
    #              'VRPutamen' #_space-MNI152NLin2009cAsym.nii.gz
    #              ]

    #print('firstlevel: ', firstleveldir)

    #################
    #Vgm_nii = nib.load(os.path.join(masks,'GrayMattermask_thalamus_space-MNI152NLin2009cAsym.nii.gz'))
    #read vol
    #Vgm_vol = Vgm_nii.get_fdata()
    #save origianl dimensions (voxels_x, voxels_y, voxels_z)
    #dim3d = Vgm_vol.shape
    #reshape to 2D
    #Vgm_2d = Vgm_vol.reshape(-1, np.prod(dim3d)).T  # -1 means auto-calculate size of dimension
    #save indexes in which Vgm == 1 (indexes for gray matter location)
    #idx_GM = np.where(Vgm_2d)[0]

    #df = pd.read_csv(os.path.join(tmp,'df_covars.csv'), header=0)
    df = pd.read_csv(covariates,header=0)
    #p_values_list=[]
    if not os.path.exists(os.path.join(output,'pvals')):
        os.makedirs(os.path.join(output,'pvals'))
    os.chdir(os.path.dirname(tarpath))
    subprocess.run(["tar","-zxvf",tarpath,voxel])
    #this should be a function... Each voxel should run in this function independently. I'll edit a bit.
    #for seed in range(len(seednames)):
    #    for voxel in range(len(idx_GM)):
    #print(f"Voxel {voxel} out of {len(idx_GM)} for seed {seed}: {seednames[seed]}")
    try:
        file_name = f'voxel : {voxel}'
        voxel_t = pd.read_csv(voxel, header=None, names=['voxel_t'])

        df['voxel_t'] = voxel_t

        model_1 = sm.MixedLM.from_formula(
            #"voxel_t ~ 1 + age + np.square(age) + sex + fdmean",
            "voxel_t ~ 1 + age + sex + fdmean",
            data=df,
            groups=df["site"],
            re_formula="1",
            vc_formula={"ID": "0 + C(ID)"},
            missing='drop'  # Omit rows with missing values
        )

        model_0 = sm.MixedLM.from_formula(
            "voxel_t ~ 1 + sex + fdmean ",
            data=df,
            groups=df["site"],
            re_formula="1",
            vc_formula={"ID": "0 + C(ID)"},
            missing='drop'  # Omit rows with missing values
        )

        # Fit models
        result_1 = model_1.fit()
        result_0 = model_0.fit()

        # View the model summary
        #print(result.summary())
        pval_age = result_1.pvalues['age']
        #pval_age2 = result_1.pvalues['np.square(age)']
        pval_age2 = np.nan
        beta_age = result_1.params['age']
        #beta_age2 = result_1.params['np.square(age)']
        beta_age2 = np.nan
        #COMPARE THE TWO MODELS
        # Calculate the likelihood ratio test (LRT) statistic
        lr_statistic = 2 * (result_1.llf - result_0.llf)
        # Calculate the degrees of freedom
        dof = result_1.df_modelwc - result_0.df_modelwc
        # Calculate the p-value using the chi-squared distribution
        pval_modelvsnull = 1 - chi2.cdf(lr_statistic, dof)

        df = df.drop(columns=['voxel_t'])

        filevals = [pval_modelvsnull, beta_age, pval_age, beta_age2, pval_age2]


        file_name = os.path.basename(voxel)
        file_path = os.path.join(output,'pvals',file_name)


        with open(os.path.join(file_path), 'w') as file:
            # Write column headers
            file.write("pval_modelvsnull\tbeta_age\tpval_age\tbeta_age2\tpval_age2\n")
            # Write p-values in different columns
            for val in filevals:
                if val == filevals[-1]:
                    file.write(f"{val}")
                else:
                    file.write(f"{val}\t")  # Write p-values
        print(f'Saved {file_path}')

    except LinAlgError as e:  # Handle the LinAlgError (singular matrix)
        # Handle the LinAlgError exception
        print(f"LinAlgError occurred: {e} at voxel {voxel}. P-val = NaN")
        p_value = np.nan
        #p_values_list.append(p_value)

        df = df.drop(columns=['voxel_t'])

        file_name = os.path.basename(voxel)
        file_path = os.path.join(output,'pvals',file_name)

        with open(file_path, 'w') as file:
           file.write(str(p_value))
        print(f'Saved {file_name}')
        subprocess.run(["rm","-v",voxel])

#%%

if __name__ == '__main__':
    main()





