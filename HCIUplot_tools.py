# -*- coding: utf-8 -*-
'''
Author: Francesco Dell'Aira
Created on: Nov 9, 2023
Last Edit: June 23, 2024
Python version: 3.9.13

Suite of function for plotting the results presented in Dell'Aira & Meier (2024) 

These functions represent an asset associated with the research work by 
F. Dell'Aira and C. I. Meier (2024), titled: "Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urbanizing Watersheds".

REFERENCES
F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urbanizing Watersheds".

DEPENDENCIES: 
1) Numpy Python library (https://numpy.org/)
2) Matplotlib Python library (https://matplotlib.org/)
                           
CREDITS:
These scripts use functions from the open-source Numpy and Matplotlib Python libraries. 

https://numpy.org/
https://matplotlib.org/


'''


#%%===============IMPORTS==========================
import numpy as np
import matplotlib.pyplot as plt

#%%===============FUNCTIONS===============

def get_model_vars (model_str):
    '''
    get the names of the model variables from the model identification strings  
    '_A_TIA' and '_A_HCIU'
    
    Parameters
    ----------
    model_str : str
        either '_A_TIA' or '_A_HCIU'
    
    Returns
    -------
    v_vars : list of str
        list of model variable names
    
    '''
    v_vars = []
    var_name = ''
    for el in model_str:
        if el == '_':
            if var_name != '':
                v_vars.append(var_name)
                var_name = ''
        else:
            var_name = var_name+el
    if not var_name in v_vars:
        v_vars.append(var_name)
    return v_vars


         
                                        
 
def plot_FullDataset_errors (df_error, 
                            caseStudy,
                            weight_str, 
                            error_measure = 'R2_adj',
                            Qstr = ['Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500'],
                            fig_hsize = 7, fig_vsize=10):
    '''
    plot and compare error metrics for the Q~A+TIA and Q~A+HCIU models, calculated by
    the pred_errors function (from the HCIUfit_test_tools module), when its keyword 
    argument dataset_type is set to 'FullDataset'. 

    Parameters
    ----------
    df_error : pandas dataframe
        output of the pred_errors function (from the HCIUfit_test_tools module), 
        when its keyword argument dataset_type is set to 'FullDataset'. 
    caseStudy : str
        label (either 'EPAE', 'MO', or 'VA') of the case study you want to consider.
    weight_str : str
        eiher 'n' or 'CN', depending on whether you test HCIU(n) or HCIU(CN).
    error_measure : str, optional
        Error metric that you want to plot (choose among "R2" and "R2_adj"). The default is 'R2_adj'.
    Qstr : list of str, optional
        names of the flood quantile columns. The default is ['Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500'].
    fig_hsize : int, optional
        horizontal dimension of the output plot. The default is 7.
    fig_vsize : int, optional
        horizontal dimension of the output plot. The default is 10.

    Returns
    -------
    Plot : Plot in the console

    '''
    plt.figure(dpi=300, figsize=(fig_hsize, fig_vsize))
    if weight_str == 'n':
        v_marker = ['x', '^']
        color = 'blue'
    else:
        v_marker = ['x', 'o']
        color = 'orange'
    v_model = ['_A_TIA', '_A_HCIU']
    x_ticks = np.array([i+1 for i in range(len(Qstr))])
    x_labels = Qstr
    
    ix = 0
    for model in v_model:
        row = df_error [df_error['model']==model].copy()
        
        columns = [ (error_measure + '_' + QQ) for QQ in Qstr]
        model_vars = get_model_vars(model)
        if model_vars[1] == 'HCIU':
            model_vars[1] = model_vars[1]+'('+weight_str+')'
        label = '$Q_T$~$' + model_vars[0] + '$+$' + model_vars[1] + '$'
        if model == '_A_TIA':
            plt.plot( x_ticks, row[columns].values[0],  label=label,  marker=v_marker[ix], markersize=15, c='k')
        else:
            plt.scatter( x_ticks, row[columns].values[0],  label=label,  marker=v_marker[ix], s=150, c=color)
        ix +=1 
        
        
    plt.xticks(ticks=x_ticks, labels=x_labels, fontsize=18)
    plt.yticks(fontsize=18)
    plt.title(caseStudy+' '+'Full Dataset, TIA vs. HCIU('+weight_str+')')
    
    if caseStudy == 'MO':
        plt.legend(fontsize=20, loc='lower right')
    else:
        plt.legend(fontsize=20, loc='lower left')
        
    if error_measure=='R2':
        plt.ylabel('$R^2$', fontsize=24)
    
    elif error_measure=='R2_adj':
        plt.ylabel('$R^2_{adj}$', fontsize=24)
        
    plt.grid()
    
    
def kfold_error_boxplots (df_test_error, 
                          caseStudy,
                          weight_str, 
                          error_measure = 'R2_adj',
                          Qstr = ['Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500'],
                          fig_hsize=10, fig_vsize=20):
    
    '''
    plot and compare boxplots of error metrics for the Q~A+TIA and Q~A+HCIU models, 
    calculated by the pred_errors function (from the HCIUfit_test_tools module), 
    when its keyword argument dataset_type is set to 'FullDataset'. 

    Parameters
    ----------
    df_test_error : pandas dataframe
        2nd output of the pred_errors function (from the HCIUfit_test_tools module), 
        when its keyword argument dataset_type is set to 'Kfold'. 
    caseStudy : str
        label (either 'EPAE', 'MO', or 'VA') of the case study you want to consider.
    weight_str : str
        eiher 'n' or 'CN', depending on whether you test HCIU(n) or HCIU(CN).
    error_measure : str, optional
        Error metric that you want to plot (choose among "R2" and "R2_adj"). The default is 'R2_adj'.
    Qstr : list of str, optional
        names of the flood quantile columns. The default is ['Q2', 'Q5', 'Q10', 'Q25', 'Q50', 'Q100', 'Q500'].
    fig_hsize : int, optional
        horizontal dimension of the output plot. The default is 10.
    fig_vsize : int, optional
        horizontal dimension of the output plot. The default is 20.

    Returns
    -------
    Plot : Plot in the console

    '''
    
    
    plot_order = ['_A_TIA', '_A_HCIU']
    v_model = ['_A_TIA', '_A_HCIU']
    boxplot_width = 0.8  # 0.8
    subgroup_dist = (len(plot_order)+2)*boxplot_width
    tick_register = []
    label_register = []
    
    
    colors = ['lightblue', 'lightgreen', 'lightyellow']
    
    
    fig, ax = plt.subplots(dpi=300, figsize=(fig_hsize, fig_vsize))  
    
    label_register = []
    
    QQix = -1
    for QQ in Qstr:
        QQix +=1
        pick = [error_measure+'_'+str(QQ)]
        col_count = 0
        for model in v_model:
            model_vars = get_model_vars(model)
            if model_vars[1] == 'HCIU':
                model_vars[1] = model_vars[1]+'('+weight_str+')' 
            label = '$'+QQ+'$~$' + model_vars[0] + '$+$' + model_vars[1] + '$'
            df_sel = df_test_error[df_test_error['model'] == model]
            position = [i for i in range(len(plot_order)) if plot_order[i]==model]
            position1 = [position[0] + QQix*subgroup_dist]
            tick_register.append(position1[0])
            bplot = ax.boxplot(df_sel[pick], vert=True, positions=position1 , widths = boxplot_width, patch_artist=True) 
            col_count+subgroup_dist
            label_register.append(label)
            bplot['boxes'][0].set_facecolor(colors[col_count])
            col_count +=1
            
    # plt.legend(fontsize=10)
    if error_measure=='R2':
        plt.ylabel('$R^2$', fontsize=24)
    
    elif error_measure=='R2_adj':
        plt.ylabel('$R^2_{adj}$', fontsize=24)
    
   
    ax.set_title(caseStudy+' Test, TIA vs. HCIU(' + weight_str + ')', fontsize=30)
    plt.grid()    
    
    
    plt.yticks(fontsize=18)
    
    
    ax.set_xticks(tick_register)
    ax.set_xticklabels(label_register, rotation=90, fontsize=15)
    if caseStudy == 'EPAE':
        ax.set_ylim(-0.5, 1)   # -1.5, 1
    elif caseStudy == 'MO':
        ax.set_ylim(-2.5, 1)
    elif caseStudy == 'VA':
        ax.set_ylim(0, 1)
        
    return 0 
    
        
        
            
            
            
