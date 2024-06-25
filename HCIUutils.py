# -*- coding: utf-8 -*-
'''
Author: Francesco Dell'Aira
Created on: May 2, 2024
Last Edit: June 13, 2024
Python version: 3.9.13

Suite of auxiliary functions for the computations of the connectivity index,
normalized connectivity index, and HCIU (Dell'Aira & Meier 2024) 

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
2) Rasterio Python library (https://rasterio.readthedocs.io/en/stable/intro.html)
                           
CREDITS:
These scripts use functions from the open-source Numpy and Rasterio Python libraries. 

https://numpy.org/
https://rasterio.readthedocs.io/en/stable/intro.html


'''
#%%==================IMPORTS===================
import numpy as np
import rasterio as rio
#%%===============FUNCTIONS=================
def get_direction (dirLabel, c_row, c_col):
    '''
    Parameters
    ----------
    dirLabel : int
        FlowDir cell value, using D8 algorithm, and ESRI convention for 
        direction labels.
    c_row : int
        row index of current FlowDir raster cell.
    c_col : int
        column index of current FlowDir raster cell.

    Returns
    -------
    orientation : float
        1 or sqrt(2), depending on cell travel orientaion (straight vs. diagonal).
    next_row : int
        row index of adjacent (destination) FlowDir raster cell
    next_col : int
        column index of adjacent (destination) FlowDir raster cell

    NOTES: ESRI convention consists of the following direction labels:
        
        :-----:-----:-----:
        : 32  : 64  : 128 :
        :-----:-----:-----:
        : 16  :cell :  1  :
        :-----:-----:-----:
        :  8  :  4  :  2  :
        :-----:-----:-----:
            
    More info about ESRI convention at:      
    https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-direction.htm
    '''
    
    if dirLabel == 128 or dirLabel == 2 or dirLabel == 8 or dirLabel == 32:
        orientation = np.sqrt(2)  #'diagonal'
    elif dirLabel == 1 or dirLabel == 4 or dirLabel == 16 or dirLabel == 64:
        orientation = 1 # 'straight'
    else:
        orientation = -1 # 'unrecognized'
        
    if dirLabel == 128:
        next_row = c_row-1
        next_col = c_col+1
    elif dirLabel == 1:
        next_row = c_row
        next_col = c_col+1
    elif dirLabel == 2: 
        next_row = c_row+1
        next_col = c_col+1
    elif dirLabel == 4: 
        next_row = c_row+1
        next_col = c_col
    elif dirLabel == 8: 
        next_row = c_row+1
        next_col = c_col-1 
    elif dirLabel == 16: 
        next_row = c_row
        next_col = c_col-1
    elif dirLabel == 32: 
        next_row = c_row-1
        next_col = c_col-1
    elif dirLabel == 64:
        next_row = c_row-1
        next_col = c_col
    else:
        next_row, next_col = -1, -1
        
    return orientation, next_row, next_col


def convert_to_tif(data_dir, file_name):
    '''
    convert raster data with given format to a GeoTIFF raster
    

    Parameters
    ----------
    data_dir : str
        location of the input raster file.
    file_name : str
        file name of the input raster.

    Returns
    -------
    the same input raster, but in GeoTIFF format : file saved in the hard-disk

    '''
    _, ext_len = get_file_ext(file_name)
    with rio.open(data_dir+file_name) as src:
        profile = src.profile
        profile.update(driver='GTiff')  # Set desired options
        with rio.open(data_dir+file_name[0:-ext_len]+'.tif', 'w', **profile) as dst:
            dst.write(src.read())
            
def get_file_ext(input_file):
    '''
    get file extension and the length of the extension string. 
    E.g., for the '.tif' extension, the length of the string is 4.

    Parameters
    ----------
    input_file : str
        file name.

    Returns
    -------
    str
        extension str.
    int
        lenght of extension str.

    '''
    end = len(input_file)-1
    scroll = 0
    ext = ''
    while input_file[end-scroll] !='.':
        ext = ext+input_file[end-scroll] 
        scroll +=1
    return '.'+ext[::-1], scroll+1

def look_neighborhood (rRASTER, c_x, c_y, neighborhood_size=1):
    
    '''
    Given a raster (rRASTER), and a pair of coorinates (latitude c_x and longitude c_y)
    of a location within the raster, retrieve the cell values and the 
    corresponding coordinates and indices within the (2z+1)x(2z+1) square 
    (neighborhood) centered at (c_x, c_y), where z is the "neighborhood_size" parameter
    
    Parameters
    ----------
    rRASTER : raster file opened with the rasterio library.
    c_x : (float) longitude.
    c_y : (float) latitude.
    neighborood_size: (integer), optional (default is 1)
        how many cells to take, along the horizontal and vertical 
        direction, from the center cell identified by [c_x, c_y]
        
    Returns
    -------
    neigh_cells : list of floats
        list of cell values in the (2z+1)x(2z+1) neighborhood 
    neigh_locs : list of floats
        list of coordinates of the cells in the (2z+1)x(2z+1) neighborhood 
    neigh_indices : list of ints
        list of rows and columns of the cells in the (2z+1)x(2z+1) neighborhood 
    
    
    '''
    c_row, c_col = rRASTER.index(c_x, c_y)
    n_row, n_col = rRASTER.shape
    #(for brevity, we will refer to 'neighborood_size' as 'z')    
    z = neighborhood_size
    neigh_cells = []
    neigh_locs = []
    neigh_indices = []
    aRASTER = rRASTER.read()
    nodata = rRASTER.nodata
    #this if statement checks that at least part of the (2z+1)x(2z+1) square  
    #falls within the 2D array (i.e., that there is a partial overlap at least)
    if (((c_row-z)<=n_row) and ((c_row+z) >= 0) and ((c_col-z)<=n_col) and ((c_col+z) >= 0)): 
        #nested for loops across every cell in the neighborhood overlapping the 2D array
        for i in range( max((c_row-z), 0) , min((c_row+z)+1, n_row) ):
            for j in range( max((c_col-z), 0) , min ((c_col+z)+1, n_col) ):
                if aRASTER[0, i, j] != nodata:
                    neigh_cells.append(aRASTER[0, i, j])
                    c_x, c_y = rRASTER.xy(i, j)
                    neigh_indices.append([i,j])
                    neigh_locs.append([c_x, c_y])
    
    return neigh_cells, neigh_locs, neigh_indices



def LULC_to_Manning (input_dir, input_file, out_dir, out_name, convtab, LULC_col, Manning_col):
    '''
    Generate a raster map of Manning coefficients from a raster map of  
    land-use/land-cover (LULC) categories, using the conversion table provided 
    with the "convtab" argument. Specify "LULC_col" and "Manning_col" parameters as 
    the names of the convtab columns associated with LULC categories and Manning   
    values, respectively.

    Parameters
    ----------
    input_dir : str
        Location  of LULC raster map.
    input_file : str
        File name of LULC raster map
    out_dir : str
         Location  of raster map of Manning's coefficients .
    out_name : str
        Name of raster map of Manning's coefficients.
    convtab : Pandas DataFrame
        Conversion table, with each row associated with a different LULC 
        category and its corresponding Manning's coefficient (n).
    LULC_col : str
        convtab column name of LULC categories.
    Manning_col : TYPE
        convtab column name of n values.

    Returns
    -------
    Raster map of Manning's roughness coefficients' : file in the hard disk

    '''     
        
    rLULC = rio.open(input_dir+input_file)
    aLULC = rLULC.read(1)
    aManning = np.zeros(aLULC.shape)
        
    for ix in convtab.index:
        LULClabel = convtab[LULC_col][ix]
        n =  convtab[Manning_col][ix]
        aManning[np.where(aLULC==LULClabel)] = n
    
    profile = rLULC.meta
    profile.update({'dtype': 'float32'})
    with rio.open(out_dir+out_name+'.tif', 'w', **profile) as dst:
        dst.write(aManning,1) 
  
    return 0


def find_mode(sample):
    '''
    find the mode in the list of numbers in the input "sample"
    variable
    
    Parameters
    ----------
    sample : list of numbers (int or float)
    
    Return
    ------
    mode of the sample : number (int or float)
    '''
    
    counts = {}
    for value in sample:
        counts[value] = counts.get(value, 0) + 1
    return max(counts, key=counts.get)



