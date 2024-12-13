# -*- coding: utf-8 -*-
"""
Author: Francesco Dell'Aira
Created on: June 10, 2024
Last Edit: June 24, 2024
Python version: 3.9.13

This script represents an asset associated with the research work by 
F. Dell'Aira and C. I. Meier (2024), titled: "Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urban Watersheds".

Tutorial 3: here we show how to test the proposed urbanization metric HCIU (Dell'Aira 
& Meier, 2024) against the traditional percentage of total impervious area (TIA). 
This script produces figures in support of the results presented in Dell'Aira and 
Meier (2024).


REFERENCES
F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urban Watersheds". Hydrology and Earth System Sciences.

DEPENDENCIES: 
1) Numpy Python library (https://numpy.org/)
2) Pandas Python library (https://pandas.pydata.org/)
3) HCIUfit_test_tools module ( HCIUfit_test_tools.py file provided along with this script)
4) HCIUplot_tools module ( HCIUplot_tools.py file provided along with this script)

    
HCIUfit_test_tools and  HCIUplot_tools modules have their own dependencies. See details in their corresponding .py files. 

                        
CREDITS:
This scrip use functions from the open-source Numpy and Pandas Python libraries. 

https://numpy.org/
https://pandas.pydata.org/
"""

#%%==============IMPORTS==================
import pandas as pd
import numpy as np
from HCIUfit_test_tools import fit_and_test_fullDataset, fit_and_test_Kfold, pred_errors
from HCIUplot_tools import plot_FullDataset_errors, kfold_error_boxplots


#%%=========== TEST HCIU(n) ===============

# For considerinf HCIU(n), we set weight_str equal to 'n'
weight_str = 'n'

# Here we choose what case study to plot, by setting the corresponding label.
# Choose one among 'EPAE', 'VA', and 'MO'
caseStudy = 'EPAE'   


# organize all basin information following the format of Table A1 in Dell'Aira & Meier (2024)
# file location:
input_dir = 'C:\\Users\\User1\\Data_folder\\'

# name of the file with the table (in this example, it is an excel file, but you can change the data format, as long as
# you load it into Python as a Pandas dataframe and follow the format as in Table A1)
input_tab = 'Data.xlsx'

# names of key dataframe columns (basin id, area, TIA, HCIU, case study labels, and flood quantiles, respectively)
gauge_id_str = 'gauge_id'
A_str = 'A'
TIA_str = 'TIA'
HCIU_str = 'HCIU(' + weight_str +')'
caseStudy_str = 'Case Study'
Qstr = [ 'Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500']


# load dataframe with basin data (with same format as Table A1 in Dell'Aira & Meier, 2024). 
# Specifically, the table must inlcude, for each basin: basin id, area, TIA, HCIU, and flood quantiles.
df = pd.read_excel(input_dir+input_tab, 
                   converters = {gauge_id_str : str},
                   index_col = 0)


# location where you want to save fitted model parameters and model predictions
output_dir = 'C:\\Users\\user1\\ModelErrors_folder\\'

# fit regional peak-flow models Q~A+TIA and Q~A+HCIU on the full dataset, then 
# calculate and save model predictions
fit_and_test_fullDataset (df, caseStudy,
                          weight_str, # 'n' or 'CN'
                          caseStudy_str,
                          A_str, TIA_str, HCIU_str, gauge_id_str,
                          output_dir,
                          Qstr = Qstr)

# the search directory for next function must be the same directory where we have
# just saved model predictions 
search_dir = output_dir

# we set this parameter to 'FullDataset' in order to calculate errors of the models fitted on the full dataset. 
dataset_type = 'FullDataset'

# calculate R2 and adjusted R2 of model predictions obtained when fitting the regional
# equation on the full dataset
df_error = pred_errors (caseStudy, # MO, VA, or EPAE
                        dataset_type, # 'FullDataset' or 'Kfold'
                        weight_str, 
                        search_dir, 
                        Qstr = Qstr)

# plot and compare adjusted R2 for the Q~A+TIA and Q~A+HCIU models
plot_FullDataset_errors (df_error, 
                         caseStudy,
                         weight_str, 
                         error_measure = 'R2_adj',
                         Qstr = Qstr,
                         fig_hsize = 7, fig_vsize=10)


# below we set the number of folds (in a k-fold validation framework) that we 
# considered in our work, to blind test Q~A+TIA and Q~A+HCIU models using different 
# random subsamples into a training and a test set. 

v_kf = [3, 4, 5]

# below we set the number of random seeds that we considered in our work, to 
# repeat the k-fold validation 10 times, each time with a different random smpling 
# of the basins into a range of training and test sets. This is to avoid any 
# potential bias in the comparisons of the two models associated with a sinlge sampling. 

v_rs = np.linspace(1, 10, 10)


# repeat the call to the fit_and_test_Kfold function for each combination of 
# the number of folds and the random seed. 
for n_KF in v_kf:
    for random_seed in v_rs :
        
        # fit regional peak-flow models Q~A+TIA and Q~A+HCIU on each training and 
        # test sets obtained following a the k-fold validation approach, for given 
        # number of folds and random seed to generate the samplings; then 
        # calculate and save model predictions (here we are saving in the same folder 
        # as when we fitted and tested the models on the full dataset.)
        fit_and_test_Kfold (df, caseStudy, 
                            weight_str, # 'n' or 'CN'
                            caseStudy_str, 
                            A_str, TIA_str, HCIU_str, gauge_id_str, 
                            n_KF, int(random_seed), 
                            output_dir,
                            Qstr = Qstr)

# calculate R2 and adjusted R2 of model predictions obtained when fitting and testing 
# regional models Q~A+TIA and Q~A+HCIU on the different training and test sets obtained 
# in the k-fold validation framework. We use the same function pred_errors previosuly 
# called, but this time we specify that we want to work with the k-fold predictions
# by setting dataset_type keyword to 'Kfold'. The function returns two dataframes 
# with error measures calculated on the training and test set, respectively.
dataset_type = 'Kfold'
_, df_error_test = pred_errors (caseStudy, # MO, VA, or EPAE
                                dataset_type, # 'FullDaset' or 'Kfold'
                                weight_str, 
                                search_dir, 
                                Qstr = Qstr)

# plot the boxplots of model test errors obtained within the k-fold validation 
kfold_error_boxplots (df_error_test, 
                      caseStudy,
                      weight_str, 
                      error_measure = 'R2_adj',
                      Qstr = Qstr,
                      fig_hsize=10, fig_vsize=20)
   


#%%=========== TEST HCIU(CN) ===============

# In this section we repeat exactly the same calculations and produce the same plots as 
# above, but considering the HCIU(CN) formulation instead (Dell'Aira & Meier, 2024). 

# Please refer to comments in the section above for details about the following lines . 

weight_str = 'CN'
HCIU_str = 'HCIU(' + weight_str +')'

fit_and_test_fullDataset (df, caseStudy,
                          weight_str, # 'n' or 'CN'
                          caseStudy_str,
                          A_str, TIA_str, HCIU_str, gauge_id_str,
                          output_dir,
                          Qstr = Qstr)

search_dir = output_dir
dataset_type = 'FullDataset'
df_error = pred_errors (caseStudy, # MO, VA, or EPAE
                        dataset_type, # 'FullDataset' or 'Kfold'
                        weight_str, 
                        search_dir, 
                        Qstr = Qstr)

plot_FullDataset_errors (df_error, 
                         caseStudy,
                         weight_str, 
                         error_measure = 'R2_adj',
                         Qstr = Qstr,
                         fig_hsize = 7, fig_vsize=10)


v_kf = [3, 4, 5]
v_rs = np.linspace(1, 10, 10)

for n_KF in v_kf:
    for random_seed in v_rs :
        
        fit_and_test_Kfold (df, caseStudy, 
                            weight_str, # 'n' or 'CN'
                            caseStudy_str, 
                            A_str, TIA_str, HCIU_str, gauge_id_str, 
                            n_KF, int(random_seed), 
                            output_dir,
                            Qstr = Qstr)
        
dataset_type = 'Kfold'
_, df_error_test = pred_errors (caseStudy, # MO, VA, or EPAE
                                dataset_type, # 'FullDaset' or 'Kfold'
                                weight_str, 
                                search_dir, 
                                Qstr = Qstr)
        
kfold_error_boxplots (df_error_test, 
                      caseStudy,
                      weight_str, 
                      error_measure = 'R2_adj',
                      Qstr = Qstr,
                      fig_hsize=10, fig_vsize=20)


