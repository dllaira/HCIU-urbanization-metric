# -*- coding: utf-8 -*-
"""          
Author: Francesco Dell'Aira
Created on: May 2, 2024
Last Edit: June 13, 2024
Python version: 3.9.13

Suite of functions for calculating the normalized connectivity index (normalized HCI), and 
the hydrologic-connectivity-based index of urbanization HCIU, both proposed by 
Dell'Aira and Meier (2024)

These functions represent an asset associated with the research work by 
F. Dell'Aira and C. I. Meier (2024), titled: "Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urban Watersheds". Hydrology and Earth System Sciences.


REFERENCES
F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urban Watersheds".

DEPENDENCIES: 
1) Numpy Python library (https://numpy.org/)
2) Rasterio Python library (https://rasterio.readthedocs.io/en/stable/intro.html)
3) StreamNetAnalysis_tools module ( StreamNetAnalysis_tools.py file provided along with this script)
4) HydroConnectAnalysis_tools module ( HydroConnectAnalysis_tools.py file provided along with this script)
5) HydroRasterAnalysis_tools module ( HydroRasterAnalysis_tools.py file provided along with this script)
6) os Default Python library (https://docs.python.org/3/library/os.html)
                           
CREDITS:
These scripts use functions from the open-source Numpy and Rasterio Python libraries. 

https://numpy.org/
https://rasterio.readthedocs.io/en/stable/intro.html

"""

#%%==================IMPORTS===================
from StreamNetAnalysis_tools import pourPoint_mapping
from HydroConnectAnalysis_tools import downslope_comp_f
from HydroRasterAnalysis_tools import D8_WFlowAcc_f
import os
import numpy as np
import rasterio as rio


#%%================FUNCTIONS====================
    
    
def normHCI_f(flowDir_dir, flowDir_file,
              flowAcc_dir, flowAcc_file, 
              W_dir, W_file,
              S_dir, S_file, 
              StrmNet_dir, StrmNet_file, 
              out_dir, out_name, 
              cell_size, W_imp):
    '''
    Calculate the raster of Normalized Hydrologic Connectity Indices proposed by  
    Dell'Aira and Meier (2024)
    
    REFERENCES
    F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
    A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
    Characterizing Urbanizing Watersheds".

    Parameters
    ----------
    flowDir_dir : str
        location of D8 Flow Direction raster.
    flowDir_file : str
        file name of D8 Flow Direction raster. Direction labels must follow the 
        ESRI convention: 
            
                    :-----:-----:-----:
                    : 32  : 64  : 128 :
                    :-----:-----:-----:
                    : 16  :cell :  1  :
                    :-----:-----:-----:
                    :  8  :  4  :  2  :
                    :-----:-----:-----:
        
    More info about the ESRI convention available at:      
    https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-direction.htm
    
               
    flowAcc_dir : str
        location of Flow Accumulation raster.
    flowAcc_file : str
        filename of Flow Accumulation raster.
    W_dir : str
        location of the raster of Weights for calculating HCI (Dell'Aira and Meier, 2024).
    W_file : str
        filename of the raster of Weights for calculating HCI (Dell'Aira and Meier, 2024).
    S_dir : str
        location of the raster of Slopes.
    S_file : str
        filename of the raster of Slopes.
    StrmNet_dir : str
        location of the raster of the stream network (with stream net cells labeled as 1, 
        0 otherwsie).
    StrmNet_file : str
        filename of the raster of the stream network (with stream net cells labeled as 1, 
        0 otherwsie).
    out_dir : str
        location where you want to save the output Normalized Connectivity raster.
    out_name : str
        file name of the output Normalized Connectivity raster.
    cell_size : int or float 
        Size (in meters) of a DEM (or Flow Accumulation) cell.
    W_imp : float
        Weight associated with fully developed cells in the totally impervious benchmark 
        (see Dell'Aira and Meier, 2024, for more details).

    Returns
    -------
    Normalized Connectivity raster : file saved in the hard-disk
        See Dell'Aira and Meier (2024) for the definition of Normalized Connectivity.

    '''

    

    aDdn = downslope_comp_f (flowAcc_dir, flowAcc_file, 
                             flowDir_dir, flowDir_file,
                             S_dir, S_file, 
                             W_dir, W_file,  
                             StrmNet_dir, StrmNet_file, 
                             cell_size, 
                             all_imp = 0)
    
    aDdn_imp = downslope_comp_f (flowAcc_dir, flowAcc_file, 
                                 flowDir_dir, flowDir_file,
                                 S_dir, S_file, 
                                 W_dir, W_file,  
                                 StrmNet_dir, StrmNet_file, 
                                 cell_size, 
                                 W_imp=W_imp, all_imp = 1)
    
    D8_WFlowAcc_f(flowDir_dir, flowDir_file, W_dir, W_file, out_dir, 'WFlowAcc')
    rWFlowAcc = rio.open(out_dir + 'WFlowAcc.tif')
    aWFlowAcc = rWFlowAcc.read(1)
    rFlowAcc = rio.open(flowAcc_dir + flowAcc_file)
    aFlowAcc = rFlowAcc.read(1) 
    aWbar = aFlowAcc.copy()
    aWbar[np.where(aFlowAcc>0)] = aWFlowAcc[np.where(aFlowAcc>0)] / aFlowAcc[np.where(aFlowAcc>0)]
    rWFlowAcc.close()
    os.remove( out_dir + 'WFlowAcc.tif')
    aRASTER = rFlowAcc.read()
    aNormHCI = aFlowAcc.copy()
    aNormHCI[np.where(aFlowAcc>0)] = aWbar[np.where(aFlowAcc>0)]*aDdn_imp[np.where(aFlowAcc>0)]/aDdn[np.where(aFlowAcc>0)]/W_imp
    # correct any local numerical errors 
    aNormHCI[np.where(aFlowAcc>0)] = np.where(aNormHCI[np.where(aFlowAcc>0)]>1, 1, aNormHCI[np.where(aFlowAcc>0)])
    aNormHCI[np.where(aFlowAcc>0)] = np.where(aNormHCI[np.where(aFlowAcc>0)]<0, 0, aNormHCI[np.where(aFlowAcc>0)])
    aRASTER[0,:,:] = aNormHCI
    with rio.open(out_dir+out_name+'.tif', 'w', **rFlowAcc.profile) as dst:
        dst.write(aRASTER)
        
    return 0
    





def HCIU_f(flowDir_dir, flowDir_file, 
           flowAcc_dir, flowAcc_file,
           normHCI_dir, normHCI_file, 
           StrmNet_dir, StrmNet_file,
           StrmTrvlT_dir, StrmTrvlT_file, 
           cell_size):
    '''
    Calculate the Hydrologic-Connectivity-based Index of Urbanization (HCIU) proposed by  
    Dell'Aira and Meier (2024), from the raster of Normalized Hydrologic Connectivity 
    Indices (Dell'Aira and Meier, 2024)
    
    REFERENCES
    F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
    A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
    Characterizing Urbanizing Watersheds".

    Parameters
    ----------
    flowDir_dir : str
        location of D8 Flow Direction raster.
    flowDir_file : str
        file name of D8 Flow Direction raster. Direction labels must follow the 
        ESRI convention: 
            
                    :-----:-----:-----:
                    : 32  : 64  : 128 :
                    :-----:-----:-----:
                    : 16  :cell :  1  :
                    :-----:-----:-----:
                    :  8  :  4  :  2  :
                    :-----:-----:-----:
        
    More info about the ESRI convention available at:      
    https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-direction.htm
    flowAcc_dir : str
        location of Flow Accumulation raster.
    flowAcc_file : str
        file name of Flow Accumulation raster.
    normHCI_dir : str
        location of the raster of Normalized Hydrologic Connectivity Indices.
    normHCI_file : str
        file name of the raster of Normalized Hydrologic Connectivity Indices.
    StrmNet_dir : str
        location of the raster of the Stream Network (stream cells labeled as 1, 0 otherwise).
    StrmNet_file : str
        file name of the raster of the Stream Network (stream cells labeled as 1, 0 otherwise).
    StrmTrvlT_dir : str
        location of the raster of the Stream Network Travel Distances.
    StrmTrvlT_file : str
        file name of the raster of the Stream Network Travel Distances.
    cell_size : int or float
        size (in meters) of a DEM (or Flow Accumulation) cell.

    Returns
    -------
    HCIU : float
        HCIU value calculated for the basin.

    '''
    rFlowDir = rio.open(flowDir_dir + flowDir_file)
    rFlowAcc = rio.open(flowAcc_dir + flowAcc_file)
    rNormHCI = rio.open(normHCI_dir + normHCI_file)
    rStrmNet = rio.open(StrmNet_dir + StrmNet_file)
    rStrmTrvlT = rio.open(StrmTrvlT_dir + StrmTrvlT_file)
    aStrmNet = rStrmNet.read(1)
    aStrmTrvlT = rStrmTrvlT.read(1)
    aFlowAcc = rFlowAcc.read(1)
    aNormHCI = rNormHCI.read(1)
    aPourPointRow, aPourPointCol = pourPoint_mapping (rFlowDir, rFlowAcc, rStrmNet)
    maxrow, maxcol =  aFlowAcc.shape 
    aPourPointRow[np.where(  (aStrmNet>0)   | (aPourPointRow>maxrow-1) | (aPourPointCol>maxcol-1)  )]=-1
    aPourPointCol[np.where(  (aStrmNet>0)   | (aPourPointRow>maxrow-1) | (aPourPointCol>maxcol-1)  )]=-1
    vrow, vcol = np.where(  (aNormHCI>0) & (aPourPointRow>=0) & (aPourPointCol>=0) )
    rows, cols = aPourPointRow[vrow, vcol] , aPourPointCol[vrow, vcol]
    NormHCIvals = aNormHCI[vrow, vcol]
    introws = [int(val) for val in rows]
    intcols = [int(val) for val in cols]
    dists = aStrmTrvlT[introws, intcols] + cell_size
    minCELLw = 0.5
    weights = 1 - minCELLw * (dists-cell_size)/(max(dists)-cell_size) 
    HCIU = np.sum(  np.array(NormHCIvals)[np.where( (np.array(NormHCIvals)>=0) & (np.array(NormHCIvals)<=1 ))]  *  np.array(weights)[np.where( (np.array(NormHCIvals)>=0) & (np.array(NormHCIvals)<=1 )) ]) / np.sum( np.array(weights)[np.where( (np.array(NormHCIvals)>=0) & (np.array(NormHCIvals)<=1 )) ])
    return HCIU


    
    


