# -*- coding: utf-8 -*-
"""          
Author: Francesco Dell'Aira
Created on: May 2, 2024
Last Edit: June 13, 2024
Python version: 3.9.13

Suite of functions for calculating the Hydrologic Connectivity Index HCI 

These functions represent an asset associated with the research work by 
F. Dell'Aira and C. I. Meier (2024), titled: "Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urban Watersheds".

These functions implement a modified formulation of the Connectivity Index,
as compared to the traditional one (see, e.g., Hooke et al., 2021, for a 
review of uses of the traditional formulation of the Connectivity Index). 
Find information about HCI in the work: 
    F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
    A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
    Characterizing Urban Watersheds". Hydrology and Earth System Sciences.

DEPENDENCIES: 
1) Numpy Python library (https://numpy.org/)
2) Rasterio Python library (https://rasterio.readthedocs.io/en/stable/intro.html)
3) HCIUutils module ( HCIUutils.py file provided along with this script)
4) HydroRasterAnalysis_tools module ( HydroRasterAnalysis_tools.py 
   file provided along with this script)
5) os default Python library (https://docs.python.org/3/library/os.html)
                           
CREDITS:
These scripts use functions from the open-source Numpy and Rasterio Python libraries. 

https://numpy.org/
https://rasterio.readthedocs.io/en/stable/intro.html

"""

#%%==================IMPORTS===================
from HCIUutils import get_direction, find_mode
from HydroRasterAnalysis_tools import D8_WFlowAcc_f
import os
import numpy as np
import rasterio as rio
from rasterio.enums import Resampling


#%%================FUNCTIONS====================


def HCI_f(S_dir, S_file, W_dir, W_file, flowDir_dir, flowDir_file, flowAcc_dir, flowAcc_file, StrmNet_dir, StrmNet_file, out_dir, out_name, cell_size): 
    '''
    This function calculates the raster of Hydrologic Connectivity Indices,
    using the formulation proposed by F. Dell'Aira and C. I. Meier (2024).
    
    F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
    A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
    Characterizing Urbanizing Watersheds".
    
    
    NOTES:
        All input raster maps must have the same size. 
        
    Parameters
    ----------
    S_dir : str
        Location of the raster of Slopes.
    S_file : str
        File name of the raster of Slopes.
    W_dir : str
        Location of the raster of Weights.
    W_file : str
        File name of the raster of Weights.
    FlowDir_dir : str
        Location of the Flow Direction raster.
    FlowDir_file : str
        file name of the Flow Direction raster.
    FlowAcc_dir : str
        Location of the Flow Accumulation raster.
    FlowAcc_file : str
        file name of the Flow Accumulation raster.
    StrmNet_dir : str
        Location of the Stream Network raster (marking stream network cells 
        as 1; 0 otherwise).
    StrmNet_file : str
        File name of the Stream Network raster (marking stream network cells 
        as 1; 0 otherwise).
    out_dir : str
        location where you want to save the HCI raster.
    out_name : str
        name that you want to assign to the HCI raster file.
    cell_size : float or int
        size (in meters) of one cell in the Flow Accumulation raster.

    Returns
    -------
    raster of hydrologic connectivity indices : file saved in the hard-disk
       

    '''
    
    aDup = upslope_comp_f(W_dir, W_file, 
                          S_dir, S_file, 
                          flowDir_dir,flowDir_file, 
                          flowAcc_dir, flowAcc_file, 
                          cell_size)
       
    aDdn = downslope_comp_f (flowAcc_dir, flowAcc_file, 
                             flowDir_dir, flowDir_file,
                             S_dir, S_file, 
                             W_dir, W_file,  
                             StrmNet_dir, StrmNet_file, 
                             cell_size)
    
    rFlowAcc = rio.open(flowAcc_dir+flowAcc_file)
    aFlowAcc = rFlowAcc.read(1)
    aHCI = aFlowAcc.copy()
    aHCI[np.where(aFlowAcc>0)] = aDup[np.where(aFlowAcc>0)] / aDdn[np.where(aFlowAcc>0)]

    aRaster = rFlowAcc.read()
    with rio.open(out_dir+out_name+'.tif', 'w', **rFlowAcc.profile) as dst:
        aRaster[0,:,:] = aHCI
        dst.write(aRaster)
    return 0

def Wraster_f (input_dir, input_file, input_dir2, input_file2, out_dir, out_name, Wflag):
    
    '''
    This function creates a raster of weights W as a function of a custom input 
    parameter dependent on LULC types, such as n or CN (the latter also changes with
    hydrologic soilg group).  
    Depending on what parameter you provide as an input, set the Wflag variable 
    accordingly (Wflag=1 for n, Wflag=2 for CN). You can also define your custom function 
    in the space marked as *** define here your custom weighting function using 
    the anonymous function sintax (see examples above) ***, and then set
    Wflag=3 to use your custom function. 
    
    Irrespective of what kind of W you use, the output W raster will be 
    automatically resampled to match the resolution of the Flow Accumulation
    raster, and its nodata value will be replaced with the nodata value of the 
    Flow Accumulation raster. 
    
    Parameters
    ----------
    input_dir : str
        location of raster map of a parameter function of LULC, such as n or CN
    input_file : str
        filename of raster map of a parameter function of LULC, such as n or CN
    input_dir2 : str
        location of Flow Accumulation raster
    input_file2 : str
        file name of Flow Accumulation raster 
    out_dir : str
        output location for raster of weights W
    out_name : str
        filename for output raster of weights W  
    Wflag : int, flag variable
        Wflag=1 -> W=1-n ; Wflag=2 -> W=CN/100, Wflag=3 -> W=W (or custom func)
    
    Returns
    ----------
    Raster of Weights : file in the hard-disk
        Raster of Weights for calculating HCI. The raster is resampled
        to match the resolution of the DEM and Flow Accumulation raster, 
        for streamlining computations related with HCI. 
    '''
    Wflag = int(Wflag)
    if Wflag == 1:
        Wf = lambda n : 1-n
    elif Wflag == 2:
        Wf = lambda CN : CN/100
    elif Wflag == 3:
        # *** define here your custom weighting function using the anonymous 
        #     function sintax (see examples above) ***
        Wf = lambda W : W  # Either keep this line as it is, if your input   
                           # raster already contains the weights and you only 
                           # want to resample the raster, OR  
                           # suitably change the function if your weights are 
                           # a function of the cell values in the input raster. 
    else:
        Wflag = input('Flag value is not valid. Set either 1, 2, or 3 if you want to compute W=1-n, W=CN/100, or W=W (or use your own custom function), respectively.')
    output_str = out_dir + out_name + 'notres.tif'  
    rLULCparam = rio.open(input_dir+input_file)
    no_data = rLULCparam.nodata
    rFlowAcc = rio.open(input_dir2 + input_file2)
    aFlowAcc = rFlowAcc.read(1)
    FlowAcc_nodata = rFlowAcc.nodata
    aLULCparam = rLULCparam.read(1)
    aW = np.zeros(aLULCparam.shape)# aLULCparam.copy()  # initialization
    # if aW.dtype != aFlowAcc.dtype:
    #     aW.dtype = aFlowAcc.dtype
    aW[np.where(aLULCparam != no_data)] = Wf(aLULCparam[np.where(aLULCparam != no_data)])
    aW[np.where(aLULCparam == no_data)] = FlowAcc_nodata  # change no data value for consistency
    out_meta = rLULCparam.meta.copy()
    out_meta.update({"driver": "GTiff",
                     "dtype": aFlowAcc.dtype,
                     "nodata": FlowAcc_nodata,
                    })
    with rio.open(output_str, "w", **out_meta) as dest:
        dest.write(aW, 1)             
    src = rio.open(output_str)
    data = src.read(1, 
            out_shape=(
                aFlowAcc.shape[0],
                aFlowAcc.shape[1]
                ),
            resampling=Resampling.nearest
            )
    src.close()
    data[np.where(aFlowAcc<0)] = rFlowAcc.nodata # clean any error at the boundaries due to resampling
    data[np.where( (data<0) & (aFlowAcc>0) )] = find_mode(data[np.where(data>=0)]) # clean any error at the boundaries due to resampling
    out2_meta = rFlowAcc.meta.copy()
    output2_str = out_dir + out_name + '.tif' 
    with rio.open(output2_str, "w", **out2_meta) as dest:
        dest.write(data,1)   
    os.remove(output_str)
    return 0 



def upslope_comp_f(W_dir, W_file, S_dir, S_file, 
                   flowDir_dir, flowDir_file, 
                   flowAcc_dir, flowAcc_file,  
                   cell_size):
    '''
    Calculate the upslope component of HCI (Dell'Aira & Meier, 2024)
    
    
    F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
    A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
    Characterizing Urbanizing Watersheds".

    Parameters
    ----------
    W_dir : str
        location of the raster of weights.
    W_file : str
        file name of the raster of weights.
    S_dir : str
        location of the raster of slopes.
    S_file : str
        file name of the raster of slopes.
    flowDir_dir : str
        location of the flow direction raster.
    flowDir_file : str
        filename of the flow direction raster.
    flowAcc_dir : str
        location of the flow accumulation raster.
    flowAcc_file : str
        file name of the flow accumulation raster.
    cell_size : int or float
        size (in meters) of a DEM cell.

    Returns
    -------
    aDup : numpy array, float
        raster of the upslope components of HCI.

    '''

    D8_WFlowAcc_f(flowDir_dir, flowDir_file, W_dir, W_file, W_dir, 'WFlowAcc')
    rWFlowAcc = rio.open(W_dir + 'WFlowAcc.tif')
    aWFlowAcc = rWFlowAcc.read(1)
    rFlowAcc = rio.open(flowAcc_dir + flowAcc_file)
    aFlowAcc = rFlowAcc.read(1) 
    aWbar = aFlowAcc.copy()
    aWbar[np.where(aFlowAcc>0)] = aWFlowAcc[np.where(aFlowAcc>0)] / aFlowAcc[np.where(aFlowAcc>0)]
    D8_WFlowAcc_f(flowDir_dir, flowDir_file, S_dir, S_file, S_dir, 'SFlowAcc')
    rSFlowAcc = rio.open(S_dir + 'SFlowAcc.tif')
    aSFlowAcc = rSFlowAcc.read(1)
    aSbar = aFlowAcc.copy()
    aSbar[np.where(aFlowAcc>0)] = aSFlowAcc[np.where(aFlowAcc>0)] / aFlowAcc[np.where(aFlowAcc>0)]
    aAreaAcc = np.where(aFlowAcc>0, aFlowAcc*cell_size**2, aFlowAcc)
    aDup = aFlowAcc.copy()
    aDup[np.where(aFlowAcc>0)] = aSbar[np.where(aFlowAcc>0)]/100*aWbar[np.where(aFlowAcc>0)]*np.sqrt(aAreaAcc[np.where(aFlowAcc>0)])
    rWFlowAcc.close()
    rSFlowAcc.close()
    os.remove(S_dir + 'SFlowAcc.tif')
    os.remove(W_dir + 'WFlowAcc.tif')
    return aDup
    
    
    


def downslope_comp_f (FlowAcc_dir, FlowAcc_file, 
                      FlowDir_dir, FlowDir_file,
                      S_dir, S_file, W_dir, W_file,  
                      StrmNet_dir, StrmNet_file, cell_size, 
                      W_imp=0.99,  W_water=1, 
                      divergent_flows_dist=10**7, 
                      hillSlope_minimumSlope = 0.005,
                      all_imp = 0):
    '''
    Calculate the downslope component of HCI (Dell'Aira & Meier, 2024)
    
    
    F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
    A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
    Characterizing Urbanizing Watersheds".

    Parameters
    ----------
    FlowAcc_dir : str
        location of flow accumulation raster.
    FlowAcc_file : str
        file name of flow accumulation raster.
    FlowDir_dir : str
        location of flow direction raster.
    FlowDir_file : str
        file name of flow direction raster.
    S_dir : str
        location of raster of slopes.
    S_file : str
        file name of raster of slopes.
    W_dir : str
        location of raster of HCI weights.
    W_file : str
        filename of raster of HCI weights.
    StrmNet_dir : str
        location of raster of Stream Network (marking river cells with 1; 0 otherwise).
    StrmNet_file : str
        file name of raster of Stream Network (marking river cells with 1; 0 otherwise).
    cell_size : int or float
        size (in meters) of a DEM cell.
    W_imp : float, optional
        Weight associated with highly developed cells of the fully
        impervious benchmark (see Dell'Aira & Meier 2024). 
        The default is 0.99.
    W_water : float, optional
        Weight W assigned to stream (water) cells, used when flowpaths approach the
        stream network. The default is 1.
    divergent_flows_dist : float, optional
        It should be a large value (the default is 10**7). 
        A properly defined Flow Direction raster typically should not 
        have diverging flows. However, in the event that a flow path leaves the basin  
        by crossing the water divide somewhere that is not the outlet, hence never
        reaching the stream network and then the outlet, those cells belonging to that
        diverging flow path will get such large value, to indicate that they are 
        far away from the stream (virtually with an infinite travel distance).
    hillSlope_minimumSlope : float, optional
        Impose a minimum point slope at each basin cell, to avoid 
        numerical errors when calculating the downslope component
        (where slope is at the denominator). The default value is 
        0.005, suggested by Borselli et al. (2008). 
    all_imp : int, optional
        Flag variable. Set to 1 if you want to calculate the downslope
        component of the fully developed benchmark (Dell'Aira & Meier 2024).
        The default is 0.

    Returns
    -------
    aWTDst : numpy array, float
        Downslope component of HCI

    '''
    rFlowAcc = rio.open(FlowAcc_dir+FlowAcc_file)
    rS = rio.open(S_dir + S_file)
    rW = rio.open(W_dir + W_file)
    rFlowDir = rio.open(FlowDir_dir + FlowDir_file)
    rStrmNet = rio.open(StrmNet_dir + StrmNet_file)
    aFlowAcc_nodata = rFlowAcc.nodata
    n_row, n_col = rFlowAcc.shape
    aStreamNet = rStrmNet.read(1)
    aFlowDir = rFlowDir.read(1)
    aFlowAcc = rFlowAcc.read(1)
    aS = rS.read(1)
    aW = rW.read(1)
    if all_imp:
        aW[np.where(aFlowAcc>0)] = W_imp
    aWTDst = aFlowAcc.copy() 
    aAlreadyMapped = aFlowAcc.copy() 
    aWTDst[np.where(aWTDst>0)] = divergent_flows_dist  # initialization
    aWTDst[np.where(aStreamNet==1)] = 0  
    aW[np.where(aStreamNet==1)] = W_water  
    aAlreadyMapped[np.where(aFlowAcc>0)] = 0  
    Null_UAA_cells = np.where(aFlowAcc==1)
    for ix in range(len(Null_UAA_cells[0])):
        c_row, c_col = Null_UAA_cells[0][ix], Null_UAA_cells[1][ix]
        riverNotReached = 1 
        simultaneous_time_update = []
        while riverNotReached:
            if aAlreadyMapped[c_row, c_col] == 1:
                c_WDst = aWTDst[c_row, c_col]
                n_listRows = len(simultaneous_time_update)
                for i in range(n_listRows):
                    simultaneous_time_update[i][1] = simultaneous_time_update[i][1] + c_WDst
                riverNotReached = 0
                for c_listRow in simultaneous_time_update:
                    crow, ccol = c_listRow[0][0], c_listRow[0][1]
                    aAlreadyMapped[crow, ccol] = 1
                    aWTDst[crow, ccol] = c_listRow[1]  
            elif aAlreadyMapped[c_row, c_col] == 0:
                orientation, next_row, next_col = get_direction (aFlowDir[c_row, c_col], c_row, c_col)
                adj_next_row, adj_next_col = max(0, min(next_row, (n_row-1) ) ), max(0, min(next_col, (n_col-1) ) )
                if adj_next_row != next_row or adj_next_col != next_col:  #we are going out the basin
                    riverNotReached = 0 
                    for c_listRow in simultaneous_time_update:
                        crow, ccol = c_listRow[0][0], c_listRow[0][1]
                        aAlreadyMapped[crow, ccol] = -1
                        aWTDst[crow, ccol] = divergent_flows_dist   
                elif orientation > 0:
                    next_UAA = aFlowAcc[next_row, next_col] 
                    if next_UAA == aFlowAcc_nodata:
                        riverNotReached = 0
                        for c_listRow in simultaneous_time_update:
                            crow, ccol = c_listRow[0][0], c_listRow[0][1]
                            aAlreadyMapped[crow, ccol] = -1
                            aWTDst[crow, ccol] = divergent_flows_dist 
                    else:
                        c_dist = cell_size*orientation
                        cW = min ( max (aW[c_row, c_col], 0.05) , W_imp) # lower and upper bounds to avoid numerical errors 
                        nextW = min ( max (aW[next_row, next_col], 0.05) , W_imp)  
                        c_Slope = min ( max( aS[c_row, c_col] /100  , hillSlope_minimumSlope) , 10 ) # lower and upper bounds to avoid numerical errors 
                        next_Slope = min ( max ( aS[next_row, next_col] /100 ,  hillSlope_minimumSlope), 10 )     
                        c_WDst = (c_dist/2/c_Slope/cW+ c_dist/2/next_Slope/nextW)/3600        
                        n_listRows = len(simultaneous_time_update)
                        if n_listRows>0:
                            for i in range(n_listRows):
                                simultaneous_time_update[i][1] = simultaneous_time_update[i][1] + c_WDst
                        c_listRow = [[c_row, c_col], c_WDst]
                        simultaneous_time_update.append(c_listRow)    
                        c_row, c_col = next_row, next_col # move on
                else: #orientation == -1 #sink or divergence 
                    riverNotReached = 0 
                    for c_listRow in simultaneous_time_update:
                        crow, ccol = c_listRow[0][0], c_listRow[0][1]
                        aAlreadyMapped[crow, ccol] = -1
                        aWTDst[crow, ccol] = divergent_flows_dist      
            else: #aAlreadyMapped[c_cell_indices[0], c_cell_indices[1]] == -1
                riverNotReached = 0 
                for c_listRow in simultaneous_time_update:
                    crow, ccol = c_listRow[0][0], c_listRow[0][1]
                    aAlreadyMapped[crow, ccol] = -1
                    aWTDst[crow, ccol] = divergent_flows_dist
            if aStreamNet[c_row, c_col] == 1:  # river reached
                riverNotReached = 0
                for c_listRow in simultaneous_time_update:
                    crow, ccol = c_listRow[0][0], c_listRow[0][1]
                    aAlreadyMapped[crow, ccol] = 1
                    aWTDst[crow, ccol] = c_listRow[1]  
    return aWTDst
