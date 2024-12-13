# -*- coding: utf-8 -*-
'''
Author: Francesco Dell'Aira
Created on: Mar 7, 2022
Last Edit: June 23, 2024
Python version: 3.9.13

Suite of function for testing the new urbanization metric HCIU proposed by 
Dell'Aira & Meier (2024) against the traditional total impervious area (TIA). 

These functions represent an asset associated with the research work by 
F. Dell'Aira and C. I. Meier (2024), titled: "Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urban Watersheds".

REFERENCES
F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urban Watersheds". Hydrology and Earth System Sciences.

DEPENDENCIES: 
1) Numpy Python library (https://numpy.org/)
2) Pandas Python library (https://pandas.pydata.org/)
3) Scikit-learn Python library (https://scikit-learn.org/stable/index.html)
                           
CREDITS:
These scripts use functions from the open-source Numpy, Pandas, and Scikit-learn Python libraries. 

https://numpy.org/
https://pandas.pydata.org/
https://scikit-learn.org/stable/index.html


'''


import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
import os

#%%=====================INPUT subset definition===============================

def fit_and_test_fullDataset (df, caseStudy,
                              weight_str, # 'n' or 'CN'
                              caseStudy_str,
                              A_str, TIA_str, HCIU_str, gauge_id_str,
                              output_dir,
                              Qstr = [ 'Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500']
                              ):
    '''
   This function performs two tasks:
   1)  fit the regional linear equation below on (Q_T, A, TIA) and (Q_T, A, HCIU) data:
              Q_T ~ log(A) + U
        where Q_T : flood quantile
              A : basin drainage area
              U : generic urbanization metric (TIA and HCIU are alternatively considered here)
              
        on the basin data fed to the function through the argument df
        
   2) saves flood predictions obtained from the two (Q_T, A, TIA) and (Q_T, A, HCIU) models,
       as well as model parameters, into excel files:
           1) Q_full_regress_*weight_str*.xlsx  (flood predictions)
           2) params_A_TIA_models.xlsx  (parameters of A-TIA models)
           3) params_A_HCIU_*weight_str*_models.xlsx' (parameters of A-HCIU models)
    
    

    Parameters
    ----------
    df : pandas dataframe 
        dataframe with area, TIA, and HCIU of each basins, as well as their
        predicted flood quantiles, and their case study labels. Each row corresponds to a 
        basin, each column to a parameter. See Table A1 in Dell'Aira and Meier (2024)
        for an example on how df is organized. 
    caseStudy : str
        label (either 'EPAE', 'MO', or 'VA') of the case study you want to consider.
    weight_str : str
        eiher 'n' or 'CN', depending on whether you test HCIU(n) or HCIU(CN).
    caseStudy_str : str
        name of the column with case study lables in the df pandas dataframe.
    A_str : str
        name of the column with basin drainage areas in the df pandas dataframe.
    TIA_str : str
        name of the column with basin TIA values in the df pandas dataframe.
    HCIU_str : str
        name of the column with basin HCIU values in the df pandas dataframe.
    gauge_id_str : str
        name of the column with basin IDs in the df pandas dataframe.
    output_dir : str
        directory where you want to save the ouutput excel files with peak-flow
        equation predictions, considering the regional equation given in 
        Dell'Aira & Meier (2024): Q_T ~ A + U, where Q_T is a flood quantile, 
        A is basin area, and U is a generic urbanization metric (TIA or HCIU here)
    Qstr : list of strings, optional
        names of the columns with flood quantiles associated with different return periods
        in the df pandas dataframe. The default is [ 'Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500'].

    Returns
    -------
    excel files : files saved in the hard-disk of the computer
        excel files with model parameters and predictions of the regional equation
        for the range of flood quantiles.

    '''
    
    df = df[df[caseStudy_str]==caseStudy]
    log_feats = [A_str] + Qstr
    Qlogstr=[ el+'log10' for el in Qstr]
    for var in log_feats:
        df[var+'log10'] = np.log10(df[var])
    x1_feats = [A_str+'log10' , TIA_str ] 
    x2_feats = [A_str+'log10' , HCIU_str ] 
    
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        
    x1=df.get(x1_feats)
    x2=df.get(x2_feats)
  
            
    LR = LinearRegression()
  
    v_params1 = []
    v_params2 = []
       
        
    df_Qpreds_fullData = pd.DataFrame()
        
    for ix in range(len(Qstr)):
        y = df.get([Qlogstr[ix]])
        Q = df.get([Qstr[ix]])

        model1 = LR.fit(x1,y)
        params1 = np.append(model1.intercept_, model1.coef_)
        v_params1.append(params1)
        y_pred1 = model1.predict(x1)
        Q_pred1 = 10**y_pred1
        
 
        model2 = LR.fit(x2,y)
        params2 = np.append(model2.intercept_, model2.coef_)
        v_params2.append(params2)
        y_pred2 = model2.predict(x2)
    
        Q_pred2 = 10**y_pred2
  
        if ix == 0:
            df_Qpreds_fullData[gauge_id_str] = df[gauge_id_str]
    
        df_Qpreds_fullData[Qstr[ix]] = Q
        df_Qpreds_fullData[Qstr[ix]+'_A_TIA'] = Q_pred1
        df_Qpreds_fullData[Qstr[ix]+'_A_HCIU'] = Q_pred2
     
    fold_name = 'FullDataset_'+caseStudy+'\\'
    if not os.path.isdir(output_dir+fold_name):
        os.mkdir(output_dir+fold_name)
        
    dfparams1 = pd.DataFrame(v_params1, index=Qstr)
    dfparams1.to_excel(output_dir+fold_name+'params_A_TIA_models'+'.xlsx')
    
    dfparams2 = pd.DataFrame(v_params2, index=Qstr)
    dfparams2.to_excel(output_dir+fold_name+'params_A_HCIU_'+weight_str+'_models.xlsx')
    
    df_Qpreds_fullData.to_excel(output_dir+fold_name+'Q_full_regress_'+weight_str+'.xlsx')
    return 0
    
    
def count_nvars (model_str):
    '''
    Count the number of explanatory variables based on the id string of the regional 
    model (either '_A_TIA' or '_A_HCIU'). Used within the steps for calculating
    the adjusted R square. 

    Parameters
    ----------
    model_str : str
        either '_A_TIA' or '_A_HCIU'.

    Returns
    -------
    kk : int
        number of explanatory vatiables.

    '''
      
    kk = 0
    
    for ch in model_str:
        if ch == '_':
            kk+=1
    
    return kk
            

def fit_and_test_Kfold (df, caseStudy, 
                        weight_str, # 'n' or 'CN'
                        caseStudy_str, 
                        A_str, TIA_str, HCIU_str, gauge_id_str, 
                        n_KF, random_seed, 
                        output_dir,
                        Qstr = [ 'Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500']
                        ):
    '''
    Like the fit_and_test_fullDataset function, but repeated following a
    k-fold validation framework. You can set the number of folds (through the 
    parameter n_KF), as well as the random seed to generate the k-fold resampling. 
    Different random seeds correspond to different resamplings. 
    Resampling are shuffled, and made in a way that the union of the different test 
    sets produces the full dataset without duplicates. 

    Parameters
    ----------
    df : pandas dataframe 
        dataframe with area, TIA, and HCIU of each basins, as well as their
        predicted flood quantiles, and their case study labels. Each row corresponds to a 
        basin, each column to a parameter. See Table A1 in Dell'Aira and Meier (2024)
        for an example on how df is organized. 
    caseStudy : str
        label (either 'EPAE', 'MO', or 'VA') of the case study you want to consider.
    weight_str : str
        eiher 'n' or 'CN', depending on whether you test HCIU(n) or HCIU(CN).
    caseStudy_str : str
        name of the column with case study lables in the df pandas dataframe.
    A_str : str
        name of the column with basin drainage areas in the df pandas dataframe.
    TIA_str : str
        name of the column with basin TIA values in the df pandas dataframe.
    HCIU_str : str
        name of the column with basin HCIU values in the df pandas dataframe.
    gauge_id_str : str
        name of the column with basin IDs in the df pandas dataframe.
    n_KF : int
        number of k-folds.
    random_seed : int
        random seed.
    output_dir : str
        directory where you want to save the ouutput excel files with peak-flow
        equation predictions, considering the regional equation given in 
        Dell'Aira & Meier (2024): Q_T ~ A + U, where Q_T is a flood quantile, 
        A is basin area, and U is a generic urbanization metric (TIA or HCIU here)
    Qstr : list of strings, optional
        names of the columns with flood quantiles associated with different return periods
        in the df pandas dataframe. The default is [ 'Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500'].

    Returns
    -------
    excel files : files saved in the hard-disk of the computer
        excel files with predictions of the regional equation
        for the range of flood quantiles.


    '''
    
    df = df[df[caseStudy_str]==caseStudy]
    log_feats = [A_str] + Qstr
    Qlogstr=[ el+'log10' for el in Qstr]
    for var in log_feats:
        df[var+'log10'] = np.log10(df[var])
    x1_feats = [A_str+'log10' , TIA_str ] 
    x2_feats = [A_str+'log10' , HCIU_str ] 
    
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        
    x1=df.get(x1_feats)
    x2=df.get(x2_feats)
  
    n = len(x1)
        
    LR = LinearRegression()
    kf = KFold(n_splits=n_KF, shuffle=True, random_state=random_seed)   
           

    df_Qpreds_test = pd.DataFrame()
    df_Qpreds_train = pd.DataFrame()

    for ix in range(len(Qstr)):
        y = df.get([Qlogstr[ix]])
        count=0
        for train_index, test_index in kf.split(range(n)):
            xtrain1, ytrain = x1.iloc[train_index], y.iloc[train_index]
            xtrain2, ytrain = x2.iloc[train_index], y.iloc[train_index]

            
            xtest1, ytest = x1.iloc[test_index], y.iloc[test_index]
            xtest2, ytest = x2.iloc[test_index], y.iloc[test_index]


            Q_train = 10**ytrain
            Q_test = 10**ytest
            
            #========================== A + TIA ==============================
            
            cmodel1 = LR.fit(xtrain1,ytrain)
            train_pred1 = cmodel1.predict(xtrain1)
            test_pred1 = cmodel1.predict(xtest1)
            
            Q_train_pred1 = 10**train_pred1
            Q_test_pred1 = 10**test_pred1
    
            
            
            #========================== A + HCIU =============================
            cmodel2 = LR.fit(xtrain2,ytrain)
            train_pred2 = cmodel2.predict(xtrain2)
            test_pred2 = cmodel2.predict(xtest2)
            
            Q_train_pred2 = 10**train_pred2
            Q_test_pred2 = 10**test_pred2
            
            
            basin_ids = list(df.iloc[test_index][gauge_id_str].values)
            Q_test_vals = [val[0] for val in Q_test.values]
            Q_test_pred1_vals = [val[0] for val in Q_test_pred1]
            Q_test_pred2_vals = [val[0] for val in Q_test_pred2]
           
            
            
            if count>0: # some folds may have different numbers of test basins as compared to previous folds, 
                        # if the size of the full dataset is not an exact multiple of the number of k-folds n_KF.
                        # Hence, add nans to columns to make the length of all columns the same 
                while  len(Q_test_vals) < len(df_Qpreds_test):
                    Q_test_vals = Q_test_vals + [np.nan]
                    Q_test_pred1_vals = Q_test_pred1_vals + [np.nan]
                    Q_test_pred2_vals = Q_test_pred2_vals + [np.nan]
                    basin_ids = basin_ids + [np.nan]
                
            df_Qpreds_test[gauge_id_str+'_'+Qstr[ix]+'_KF'+str(count+1)] = basin_ids
            df_Qpreds_test[Qstr[ix]+'_KF'+str(count+1)] = Q_test_vals
            df_Qpreds_test[Qstr[ix]+'_A_TIA_KF'+str(count+1)] = Q_test_pred1_vals
            df_Qpreds_test[Qstr[ix]+'_A_HCIU_KF'+str(count+1)] = Q_test_pred2_vals
            
            
            
            basin_ids = list(df.iloc[train_index][gauge_id_str].values)
            Q_train_vals = [val[0] for val in Q_train.values] 
            Q_train_pred1_vals = [val[0] for val in Q_train_pred1]
            Q_train_pred2_vals = [val[0] for val in Q_train_pred2]
          
            
            
            if count==0: # some folds may have different numbers of train basins as compared to previous folds, 
                         # if the size of the full dataset is not an exact multiple of the number of k-folds n_KF.
                         # Hence, add nans to columns to make the length of all columns the same 
                while  len(Q_train_vals) < len(train_index)+n_KF:
                    Q_train_vals = Q_train_vals + [np.nan]
                    Q_train_pred1_vals = Q_train_pred1_vals + [np.nan]
                    Q_train_pred2_vals = Q_train_pred2_vals + [np.nan]
                    basin_ids = basin_ids + [np.nan]
                
            if count>0: # some folds may have different numbers of train basins as compared to previous folds, 
                        # if the size of the full dataset is not an exact multiple of the number of k-folds n_KF.
                        # Hence, add nans to columns to make the length of all columns the same 
                while  len(Q_train_vals) < len(df_Qpreds_train):
                    Q_train_vals = Q_train_vals + [np.nan]
                    Q_train_pred1_vals = Q_train_pred1_vals + [np.nan]
                    Q_train_pred2_vals = Q_train_pred2_vals + [np.nan]
                    basin_ids = basin_ids + [np.nan]
                
            df_Qpreds_train[gauge_id_str+'_'+Qstr[ix]+'_KF'+str(count+1)]= basin_ids
            df_Qpreds_train[Qstr[ix]+'_KF'+str(count+1)] = Q_train_vals
            df_Qpreds_train[Qstr[ix]+'_A_TIA_KF'+str(count+1)] = Q_train_pred1_vals
            df_Qpreds_train[Qstr[ix]+'_A_HCIU_KF'+str(count+1)] = Q_train_pred2_vals
                      
    
            count=count+1
        
        
        
    fold_name = 'Fold_'+caseStudy+'_KF'+str(n_KF)+'_rs'+str(random_seed)+'\\'
    
    if not os.path.isdir(output_dir+fold_name):
        os.mkdir(output_dir+fold_name)
    
    df_Qpreds_test.to_excel(output_dir+fold_name+'Q_Kfold'+'_test_'+weight_str+'.xlsx')
    df_Qpreds_train.to_excel(output_dir+fold_name+'Q_Kfold'+'_train_'+weight_str+'.xlsx')
    return 0 


def get_Kfold_params(root):
    '''
    this function retrieves the number of k-folds and the randon seed from 
    the string corresponding to the folder name where model predictions are saved. 
    Folder name must match the format followed by the fit_and_test_Kfold function
    when saving excel files with model predictions. 

    Parameters
    ----------
    root : str
        directory.

    Returns
    -------
    KF, rs : int, int
        number of k-folds and random seed that were used to obtain model results
        saved in the root folder.

    '''
    
    
    for ix in range(len(root)-2):
        
        if root[ix:ix+2]=='KF':
            
            ixx = ix+2
            KF = ''
            while root[ixx].isnumeric():
                
                KF = KF + root[ixx]
                ixx+=1
                
        
        if root[ix:ix+2]=='rs':
            
            ixx = ix+2
            rs = ''
            while root[ixx].isnumeric():
                
                rs = rs + root[ixx]
                ixx+=1
                if not ixx<len(root):
                    break
        
        
    
    return int(KF), int(rs)

def pred_errors (caseStudy, # MO, VA, or EPAE
                 dataset_type, # 'FullDataset' or 'Kfold'
                 weight_str,  
                 search_dir, 
                 Qstr = ['Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500']):
    '''
    read the excel files generated either by fit_and_test_fullDataset or fit_and_test_Kfold, 
    function, and calculate an error metric. 

    Parameters
    ----------
    caseStudy : str
        label (either 'EPAE', 'MO', or 'VA') of the case study you want to consider.
   dataset_type : str
        either 'FullDataset' or 'Kfold', depending on whether you want to calculate
        errors associated with the models fitted on the full dataset (by using
        the fit_and_test_fullDataset function), or else the range of models fitted 
        within the k-fold validation framework (by using the fit_and_test_Kfold function).
    weight_str : str
        eiher 'n' or 'CN', depending on whether you test HCIU(n) or HCIU(CN).
    search_dir : str
        location where you told the fit_and_test_fullDataset or fit_and_test_Kfold function 
        to save the excel files with model predictions. 
    Qstr : list of strings, optional
        names of the columns of the flood quantiles. The default is ['Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500'].

    Returns
    -------
    df_comp : pandas dataframe (returned if dataset_type=='FullDataset')
        dataframe with error metrics associated with the predictions of the 
        Q~A+TIA and Q~A+HCIU models, for the range of distinct flood quantiles.
    
    df_comp_train, df_comp_test : two pandas dataframes (returned if dataset_type=='Kfold')
        dataframes with error metrics associated with the predictions of the 
        Q~A+TIA and Q~A+HCIU models fitted several times on different training and test 
        sets, following a k-fold validation procedure. 

    '''
    
  
    ix = 0
    
    
    if dataset_type == 'FullDataset':
        
        columns = ['model'] + ['R2_'+str(cc) for cc in Qstr] + ['R2_adj_'+str(cc) for cc in Qstr] 
        df_comp = pd.DataFrame(columns=columns)
        
        
        
        input_dir = search_dir + 'FullDataset_'+ caseStudy + '\\'
        if os.path.isdir(input_dir):
            converters = {'gauge_id':str }
            cfile = 'Q_full_regress_' + weight_str + '.xlsx'
            tab = pd.read_excel(input_dir+cfile, index_col=(0), converters=converters)
            

            for model in ['_A_TIA', '_A_HCIU']:
                
                vr2 = []
                vr2_adj = []
                
                for QQ in Qstr:
                    
                    actual_label = QQ
                    model_label = QQ + model
                    
                    act, pred = tab[actual_label].copy(), tab[model_label].copy()
                    act =act.dropna()
                    pred = pred.dropna()
                    
                    r2 = r2_score(act, pred)
                    nn = len(act)
                    kk = count_nvars (model)
                    r2_adj = 1 - (1 - r2) * (nn - 1) / (nn - kk - 1)
                   
                    vr2.append(r2)
                    vr2_adj.append(r2_adj)
                  
            
                        
                
                row = [model] + vr2 + vr2_adj 
                row = pd.DataFrame([row], columns=columns, index=[ix])
                df_comp = pd.concat([df_comp, row], axis=0)
                ix +=1 
                
        return df_comp
                
                                  
                                    
    elif dataset_type == 'Kfold':
              
        columns = ['KF', '#KF', 'rs', 'model'] + ['R2_'+str(cc) for cc in Qstr] + ['R2_adj_'+str(cc) for cc in Qstr]
        
        df_comp_train = pd.DataFrame(columns=columns)
        df_comp_test = pd.DataFrame(columns=columns)
            
        
            
        count = 0
        
        for root, cdir, files in os.walk(search_dir):
            
            count+=1
           
            if count>1: # skip root directory == search_dir
            
                if caseStudy in root and not 'FullDataset' in root:
                    
                    KF, random_seed = get_Kfold_params(root)
                    os.chdir(root)
                    cfile1 = 'Q_Kfold_test_' + weight_str + '.xlsx'
                    cfile2 = 'Q_Kfold_train_' + weight_str + '.xlsx'
                    
                    converters = {'gauge_id'+'_'+str(QQ)+'_KF'+str(cKF):str for QQ in Qstr for cKF in range(KF) }
                    tab_test = pd.read_excel(cfile1, index_col=(0), converters=converters)
                    tab_train = pd.read_excel(cfile2, index_col=(0), converters=converters)
                    
                    
                    for cKF in range(1,KF+1):
                        for model in ['_A_TIA', '_A_HCIU']:
                            
                            vr2 = []
                            vr2_adj = []
                                                         
                            for QQ in Qstr:
                                
                                actual_label = str(QQ) + '_KF' + str(cKF) 
                                model_label = str(QQ) + model + '_KF' + str(cKF)   
                                
                                act, pred = tab_test[actual_label].copy(), tab_test[model_label].copy()
                                
                                act = act.dropna()
                                pred = pred.dropna()
                                
                                                                    
                                r2 = r2_score(act, pred)
                                nn = len(act)
                                kk = count_nvars (model)
                                r2_adj = 1 - (1 - r2) * (nn - 1) / (nn - kk - 1)
                                
                                vr2.append(r2)
                                vr2_adj.append(r2_adj)
                               
                                
                    
                            
                                                        
                        
                            row_test = [KF] + [cKF] + [int(random_seed)] + [model] + vr2 + vr2_adj 
                            row_test = pd.DataFrame([row_test], columns=columns, index=[ix])  
                            df_comp_test = pd.concat([df_comp_test, row_test], axis=0)
            
                    
                    
                            vr2 = []
                            vr2_adj = []
                                                           
                            for QQ in Qstr:
                                
                                actual_label = str(QQ) + '_KF' + str(cKF)  
                                model_label = str(QQ) + model + '_KF' + str(cKF)  
                                
                                act, pred = tab_train[actual_label].copy(), tab_train[model_label].copy()
                                
                                act =act.dropna()
                                pred = pred.dropna()
                                
                                                                      
                                r2 = r2_score(act, pred)
                                
                                nn = len(act)
                                kk = count_nvars (model)
                                r2_adj = 1 - (1 - r2) * (nn - 1) / (nn - kk - 1)
                                
                               
                                
                               
                                vr2.append(r2)
                                vr2_adj.append(r2_adj)
                               
                                                        
                            row_train = [KF] + [cKF] + [int(random_seed)] + [model] + vr2 + vr2_adj 
                            row_train = pd.DataFrame([row_train], columns=columns, index=[ix])
                            df_comp_train = pd.concat([df_comp_train, row_train], axis=0)
                            ix +=1 
        
        else: 
            print('dataset_type keyword not recognized. Use either "FullDataset" or "Kfold".')
                    
        return     df_comp_train, df_comp_test         

