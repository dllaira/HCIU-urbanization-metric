# -*- coding: utf-8 -*-
'''
Author: Francesco Dell'Aira
Created on: May 2, 2024
Last Edit: June 13, 2024
Python version: 3.9.13

Suite of functions for performing traditional DEM hydrologic analysis operations 
using existing open-source Python software, or implementing literature algorithms. 

These functions represent an asset associated with the research work by 
F. Dell'Aira and C. I. Meier (2024), titled: "Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urbanizing Watersheds".

DEPENDENCIES: 
1) Numpy Python library (https://numpy.org/)
2) WhiteBox Python library (https://whitebox.readthedocs.io/en/latest/readme.html#)
3) Rasterio Python library (https://rasterio.readthedocs.io/en/stable/intro.html)
4) HCIUutils module ( HCIUutils.py file provided along with this script)
                           
CREDITS:
These scripts use functions from the open-source Numpy,  
WhiteBox, and Rasterio Python libraries. 

https://numpy.org/
https://whitebox.readthedocs.io/en/latest/readme.html#
https://rasterio.readthedocs.io/en/stable/intro.html

The D8_WFlowAcc_f function in this module implements the flow accumulation 
algorithm proposed by Zhou et al. (2019)

REFERENCES:
    
    Zhou, G., Dong, W., and Wei, H.,  2019. A fast and simple algorithm for calculating flow accumulation 
    matrices from raster digital elevation models. Abstract of the International Cartographic Association 
    (https://doi.org/10.5194/ica-abs-1-434-2019)
    

'''

#%%==================IMPORTS===================
from HCIUutils import convert_to_tif, get_file_ext, get_direction
from whitebox_workflows import WbEnvironment
import rasterio as rio
import numpy as np
#%%===============GLOBAL VARS==================

global wbe
wbe = WbEnvironment()


#%%================FUNCTIONS====================

def fill_dem_f (input_dir, input_file, out_dir, out_name):
    '''
    This function takes as input a raw (i.e., with possible presence of sinks)  
    Digital Elevation Model (DEM), and returns a  depression- (or sink-)free 
    DEM. The depression-free DEM is saved in the hard-disk of the computer, in 
    the directory specified by the user.
    

    Parameters
    ----------
    input_dir : str
        Location of the input raw DEM.
    input_file : str
        Name of the input raw DEM file.
    out_dir : str
        Location of the output depression-free DEM.
    out_name : str
        Name of the output depression-free DEM file. 

    Returns
    -------
    Sink-free DEM : file in the hard-disk
    
    NOTES:
    Because this script uses functions from the WhiteBox Python module
    (https://whitebox.readthedocs.io/en/latest/readme.html#), 
    supported raster data formats depend on compatibility with that module.  
    If this function is unsuccessful with reading the input DEM, by default it will 
    automatically try to convert the original file into a compatible, GeoTIFF file (i.e., 
    with ".tif" exension), and will save it in the same folder as the original 
    input file (without overwriting). 
    

    
    '''
    
    try: 
        wbDEM = wbe.read_raster(input_dir+input_file)
    except: # if format of the input file is not compatible, convert to GeoTIFF format to ensure compatibility
        ext, len_ext = get_file_ext(input_file)
        convert_to_tif(input_dir, input_file)
        wbDEM = wbe.read_raster(input_dir+input_file[0:-len_ext]+'.tif')
        
    wbfDEM = wbe.fill_depressions(wbDEM)
    wbe.write_raster(wbfDEM, out_dir+out_name+'.tif', compress=False) 
    
    
    return 0

def Slope_f (input_dir, input_file, out_dir, out_name):
    '''
    This function takes as input a Digital Elevation Model (DEM), and returns
    a raster of slopes, expressed as percentages. The latter is saved in 
    the hard-drive of the computer, in the directory specified by the user. 
    
    
    NOTES:
    Point slopes expressed in percent
    
    In the context of HCIU computations, it is recommended feeding a 
    sink-free DEM, to be consistent with computations related with 
    flow direction and flow accumulation raster maps (which strictly 
    require using a sink-free DEM)
    
    Because this script uses functions from the WhiteBox Python module
    (https://whitebox.readthedocs.io/en/latest/readme.html#), 
    supported raster data formats depend on compatibility with that module.  
    If this function is unsuccessful with reading the input DEM, by default it will 
    automatically try to convert the original file into a compatible, GeoTIFF file (i.e., 
    with ".tif" exension), and will save it in the same folder as the original 
    input file (without overwriting). 


    Parameters
    ----------
    input_dir : str
        Location of the input DEM.
    input_file : str
        Name of the input DEM file.
    out_dir : str
        Location of the output raster of Slopes.
    out_name : str
        Name of the output raster of Slopes.

    Returns
    -------
    Raster of Slopes : file in the hard-disk
    
    '''
    try: 
        wbDEM = wbe.read_raster(input_dir+input_file)
    except: # if format of the input file is not compatible, convert to GeoTIFF format to ensure compatibility
        ext, len_ext = get_file_ext(input_file)
        convert_to_tif(input_dir, input_file)
        wbDEM = wbe.read_raster(input_dir+input_file[0:-len_ext]+'.tif')
    wbSlope = wbe.slope(wbDEM, units='percent')
    wbe.write_raster(wbSlope, out_dir+out_name+'.tif', compress=False) 
    return 0
        
def D8FlowDir_f (input_dir, input_file, out_dir, out_name):
    '''
    This function takes as input a sink-free DEM and returns a raster of 
    flow directions (using D8 algorithm). The output raster is saved in
    the hard-drive of the computer, in the directory specified by the user. 
    
    
    NOTES:
    It is strongly recommeneded only feeding a preprocessed, sink-free DEM, 
    not a raw DEM
        
    The output flow direction follows the ESRI convention for direction labels:
        
        :-----:-----:-----:
        : 32  : 64  : 128 :
        :-----:-----:-----:
        : 16  :cell :  1  :
        :-----:-----:-----:
        :  8  :  4  :  2  :
        :-----:-----:-----:
        
    More info about the ESRI convention available at:      
    https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-direction.htm
    
    Because this script uses functions from the WhiteBox Python module
    (https://whitebox.readthedocs.io/en/latest/readme.html#), 
    supported raster data formats depend on compatibility with that module.  
    If this function is unsuccessful with reading the input DEM, by default it will 
    automatically try to convert the original file into a compatible, GeoTIFF file (i.e., 
    with ".tif" exension), and will save it in the same folder as the original 
    input file (without overwriting). 
    

    Parameters
    ----------
    input_dir : str
        Location of the input depression-free DEM.
    input_file : str
        Name of the depression-free DEM input file.
    out_dir : str
        Location of the output D8-Flow Direction raster.
    out_name : str
        Name of the output D8-Flow Direction raster file. 

    Returns
    -------
    Raster of D8 Flow Directions : file in the hard-disk
    
    '''
    try: 
        wbDEM = wbe.read_raster(input_dir+input_file)
    except: # if format of the input file is not compatible, convert to GeoTIFF format to ensure compatibility
        ext, len_ext = get_file_ext(input_file)
        convert_to_tif(input_dir, input_file)
        wbDEM = wbe.read_raster(input_dir+input_file[0:-len_ext]+'.tif')
    wbFlowDir = wbe.d8_pointer(wbDEM, esri_pointer=True )
    wbe.write_raster(wbFlowDir, out_dir+out_name+'.tif', compress=False) 
    
def D8FlowAcc_f(input_dir, input_file, out_dir, out_name, include_log_transform=0):
    '''
    This function takes as input a D8-Flow Direction raster, and gives as 
    output a raster of Flow Accumulation, which is saved in the hard-drive 
    of the computer, in the directory specified by the user. 
    If "include_log_transform" is set to 1, this function also returns a 
    log-transformed Flow Accumulation raster, which may be optionally used 
    to better visualize cells with high numbers of draining cells upstream. 

    NOTES:
    Use input Flow Direction Raster that uses the ESRI convention for directions
    (see below)
        
        :-----:-----:-----:
        : 32  : 64  : 128 :
        :-----:-----:-----:
        : 16  :cell :  1  :
        :-----:-----:-----:
        :  8  :  4  :  2  :
        :-----:-----:-----:
        
    More info about the ESRI convention available at:      
    https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-direction.htm
    

    Parameters
    ----------
    input_dir : str
        Location of the input D8-Flow Direction raster.
    input_file : str
        Name of the input D8-Flow Direction raster.
    out_dir : str
        Location of the output Flow Accumulation raster.
    out_name : str
        Name of the output Flow Accumulation raster file. 
    include_log_transform : int 
        Flag variable: if 1, also generate log-transformed Flow Accumulation raster

    Returns
    -------
    Raster of Flow Accumulation : file in the hard-disk
    (Optional) Raster of log-transformed Flow Accumulation : file in the hard-disk
    
    '''
    wbFlowDir = wbe.read_raster(input_dir+input_file)
    wbFlowAcc = wbe.d8_flow_accum(wbFlowDir, input_is_pointer=True, esri_pntr=True, log_transform=False, out_type='cells')
    wbe.write_raster(wbFlowAcc, out_dir+out_name+'.tif', compress=False) 
    if include_log_transform==1:
        wbFlowAccLog = wbe.d8_flow_accum(wbFlowDir, input_is_pointer=True, esri_pntr=True, log_transform=True, out_type='cells')
        wbe.write_raster(wbFlowAccLog, out_dir+out_name+'_log.tif', compress=False) 
    
    return 0



def D8_WFlowAcc_f(flowDir_dir, flowDir_file, W_dir, W_file, out_dir, out_name):
    '''
    This function computes a weighted flow accumulation raster map, starting 
    from a D8-FLow Direction raster (using the ESRI convention for direction labels) 
    and a raster of weights with same size as the Flow Direction raster.
    This function implements the flow accumulation algorithm proposed by Zhou et al. (2019),
    slightly modified to allow for the computation of a *weighted* flow accumulation raster. 
    
    Zhou, G., Dong, W., and Wei, H.,  2019. A fast and simple algorithm for calculating flow accumulation 
    matrices from raster digital elevation models. Abstract of the International Cartographic Association 
    (https://doi.org/10.5194/ica-abs-1-434-2019)
    
    Parameters
    ----------
    flowDir_dir : str
        location of the input D8 flow direction raster.
    flowDir_file : str
        file name of the input D8 flow direction raster.
    W_dir : str
        location of the raster of weights.
    W_file : str
        file name of the raster of weights.

    Returns
    -------
    Weighted Flow Accumulation raster map : file in the hard-disk
    
    NOTES:
    The flow direction raster and the raster of weights must have the same dimensions. 
    
    ESRI convention for flow directions labels: 
        
                
        :-----:-----:-----:
        : 32  : 64  : 128 :
        :-----:-----:-----:
        : 16  :cell :  1  :
        :-----:-----:-----:
        :  8  :  4  :  2  :
        :-----:-----:-----:
        
    More info about the ESRI convention available at:      
    https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-direction.htm

    '''
    
    rFlowDir = rio.open(flowDir_dir + flowDir_file)
    rW = rio.open(W_dir + W_file)
    
    aFlowDir = rFlowDir.read(1)
     
    aW = rW.read(1)
    aWFlowAcc = aW.copy()
    
    nrow, ncol = aFlowDir.shape
    
    NIDP = aFlowDir.copy()
    NIDP[np.where(aFlowDir>-10**3)] = 0 # initialization
    vrows, vcols = np.where(aFlowDir>-10**3)
    for ix in range(len(vrows)):
        i, j = vrows[ix], vcols[ix]
        _, ii, jj = get_direction(aFlowDir[i,j], i, j)
        if NIDP[ii,jj] >=0:
            NIDP[ii,jj] =  NIDP[ii,jj]+1
    vrows, vcols = np.where(NIDP>=0)        
    for ix in range(len(vrows)):
        i, j = vrows[ix], vcols[ix]
        if NIDP[i,j] == 0 :
            path_completed = 0
            mi, mj = i, j
            while not path_completed:
                _, ii, jj = get_direction(aFlowDir[mi,mj], mi, mj)
                if aWFlowAcc[ii,jj] > -10**3:
                    aWFlowAcc[ii,jj] = aWFlowAcc[ii,jj] + aWFlowAcc[mi,mj]
                    mi, mj = ii, jj
                    if NIDP[mi, mj] > 1:
                        NIDP[mi, mj] = NIDP[mi, mj]-1
                        path_completed = 1
                else:
                    path_completed = 1
                    
    with rio.open(out_dir+out_name+'.tif', 'w', **rW.meta) as dst:
        dst.write(aWFlowAcc,1)   
    return 0
    

    
    
    