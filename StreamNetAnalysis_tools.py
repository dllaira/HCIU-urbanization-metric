# -*- coding: utf-8 -*-
"""          
Author: Francesco Dell'Aira
Created on: May 2, 2024
Last Edit: June 13, 2024
Python version: 3.9.13

Suite of functions for outlining the stream network of a watershed
from the flow accumulation raster, either using a fixed threshold 
for the number of upstream draining cells, or else considering 
actual headwater locations. 

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
2) Geopandas Python library (https://geopandas.org/en/stable/getting_started/introduction.html)
3) Rasterio Python library (https://rasterio.readthedocs.io/en/stable/intro.html)
4) Shapely Python library (https://pypi.org/project/shapely/)
5) HCIUutils module ( HCIUutils.py file provided along with this script)

                           
CREDITS:
These scripts use functions from the open-source Numpy,  
Geopandas, Shapely, and Rasterio Python libraries. 

https://numpy.org/
https://geopandas.org/en/stable/getting_started/introduction.html
https://pypi.org/project/shapely/
https://rasterio.readthedocs.io/en/stable/intro.html

"""

#%%==================IMPORTS===================
import HCIUutils
import numpy as np
import rasterio as rio
from geopandas import GeoDataFrame
import geopandas as gpd
from shapely.geometry import Point


#%%================FUNCTIONS====================

def StreamNet_f (input_dir1, input_file1, input_dir2, input_file2, out_dir, out_name, out2_dir, out2_name, cell_size, headwaters_dir = '', headwaters_file = '', headwaters=[], step = 10000, logstep=0):
    '''
    This function returns 1) a raster with 0, 1, and nodata cells, 
    where 1 indicates stream network cells, 0, indicates other basin cells, 
    and the remaining cells are nodata (taken from the Flow Accumulation 
    raster); and 2) a raster where each stream network cell displays its 
    travel distance to the outlet, measured along the stream network. 
    
    
    You can use this function in three alternative ways: 
        1) provide the list of coordinates of the headwater locations, based on, 
           e.g., the available stream-network flow lines available for your basin 
           (if applicable), or other data sources. List of locations follows the 
           format:
               headwaters = [ (lon1, lat1), (lon2, lat2), ..., (lonN, latN)] 
        
        2) same as above, but instead of a list, provide input shapefile with 
           points that you intend to use as headwaters for your basin (use 
           headwaters_dir & headwaters_file variables)
               
        3) traditional way: set a minimum threshold (step variable) for the 
           number of upstream drainage cells. This is recommended if you do not 
           have information about the stream network of your basin
    
    NOTES:
    Use input Flow Direction Raster that uses ESRI convention for directions
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

            
    If you use the traditional way:
    Sometimes it is easier to identify vallyes, tributaries, and main channel 
    from the log-transformed FLow Accumulation raster. For convenience, you can 
    optionally set the step variable based on values observed from the 
    log-transformed Flow Accumulation raster, instead of considering values 
    from the regular Flow Accumulation raster. If you wish to do so, set the 
    flag-variable logstep to 1. Note that the input FLow Accumulation raster 
    must be the non-transformed either way. 
    
    If you want to provide the headwater locations (either by a list of coordinates, or by feeding 
    the shapefile with headwater points):
    If headwater coordinates do not exactly correspond to Flow Accumulation pixels with large numbers
    of upstream drainage cells, for instance because flow lines and DEM data come from different sources, 
    you can use the "snap_headwaters_to_highFlowAcc_cells" function (provided with the 
    StreamNetAnalysis_tools.py module) to snap those headwater locations to neighboring Flow Accumulation 
    cells with large numbers of upstream draining pixels.
         

    Parameters
    ----------
    input_dir1 : str
        location of D8 flow direction raster (using ESRI convention).
    input_file1 : str
        name of D8 flow direction raster file. 
    input_dir2 : str
        location of Flow Accumulation raster.
    input_file2 : str
        name of FLow Accumulation raster file.
    out_dir : str
        location for the output raster identifying the stream network.
    out_name : str
        name of output raster identifying the stream network.
    headwaters_dir : str, optional
        location of the shapefile with headwater locations
    headwaters_file : str, optional
        name of the shapefile with headwater locations (given as points)
    headwaters : list, optional
        list of headwater locations that you want to enforce. The stream net
        will be drawn starting from those points. Follow the format:
            headwaters = [ (lon1, lat1), (lon2, lat2), ..., (lonN, latN)] 
    step : int, optional
        Traditional approach: set a unique minimum threshold 
        value for the number of upstream cells, to draw the stream network. 
        The default is 10000 cells.
    logstep : int, flag variable, optional
        (to be used in combination with step variable) set this flag to 1 if 
        you provide the step variable based on readings from the log-tranformed 
        FLow Accumulation raster. By default, this functionality is not used. 
    Returns
    -------
    Raster map of the basin with stream-setwork and non-stream-network cells
    (1 and 0 cell values, respectively): file in the hard-disk.
    
    Raster map of the stream cell distances to the outlet, as measured along 
    the stream network : file in the hard-disk.

    '''
    rFlowAcc = rio.open(input_dir2+input_file2)
    rFlowDir = rio.open(input_dir1+input_file1)
    if len(headwaters)==0 and headwaters_file=='':
        if logstep:
            step = int(np.floor(np.e**step))
        headwaters = get_headwaters(rFlowAcc, rFlowDir, step=int(step))
    elif len(headwaters)==0 and headwaters_file != '':
        headwaters = get_headwaters(rFlowAcc, rFlowDir, headwaters_dir=headwaters_dir, headwaters_file=headwaters_file)
    DEMstreamNet(rFlowDir, rFlowAcc, headwaters, out_dir, out_name)
    AlongStrmNet_dist(rFlowDir, rFlowAcc, headwaters, cell_size, out2_dir, out2_name)
    return 0



def pourPoint_mapping (rFlowDir, rFlowAcc, rStreamNet):
    
    '''
    Parameters
    ----------
    rFlowDir : D8 flow direction raster opened with the rasterio module.
        Flow Direction labels must follow the ESRI convention:
            
                :-----:-----:-----:
                : 32  : 64  : 128 :
                :-----:-----:-----:
                : 16  :cell :  1  :
                :-----:-----:-----:
                :  8  :  4  :  2  :
                :-----:-----:-----:
                    
         More info about the ESRI convention available at:      
         https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-direction.htm

        
    rFlowAcc : flow accumulation raster opened with the rasterio module.
               
    rStreamNet: raster of stream network opened with the rasterio module 
        this raster indicates stream network cells with 1's, and with 0's all other basin cells 
    
    Returns
    -------
    aDestStrCellrow : numpy array 
        each raster cell has the row number of the stream network cell (pour point) where that cell 
        will eventually drain 
    
    aDestStrCellcol : numpy array 
        each raster cell has the column number of the stream network cell (pour point) where that cell 
        will eventually drain 
        
    '''
    
    
    aFlowAcc_nodata = rFlowAcc.nodata
    n_row, n_col = rFlowAcc.shape
    
    aStreamNet = rStreamNet.read(1)
    aFlowDir = rFlowDir.read(1)
    aFlowAcc = rFlowAcc.read(1)
    
    aDestStrCellrow = aFlowAcc.copy()
    aDestStrCellrow[np.where(aFlowAcc>0)] = -1
    aDestStrCellcol = aDestStrCellrow.copy()
    
    aAlreadyMapped = aFlowAcc.copy() #to keep track of cells already  mapped 
    aAlreadyMapped[np.where(aFlowAcc>0)] = 0  # initialization: no cells are mapped
    
    Null_UAA_cells = np.where(aFlowAcc==1)
    
    for ix in range(len(Null_UAA_cells[0])):
        
        c_row, c_col = Null_UAA_cells[0][ix], Null_UAA_cells[1][ix]
        
        riverNotReached = 1 
        simultaneous_time_update = []
        
        while riverNotReached:
            
            if aAlreadyMapped[c_row, c_col] == 1:
            
      
                riverNotReached = 0
                for c_listRow in simultaneous_time_update:
                    crow, ccol = c_listRow[0], c_listRow[1]
                    aAlreadyMapped[crow, ccol] = 1
                    aDestStrCellrow[crow, ccol] = aDestStrCellrow[c_row, c_col]
                    aDestStrCellcol[crow, ccol] = aDestStrCellcol[c_row, c_col]
                
            elif aAlreadyMapped[c_row, c_col] == 0:
                
                orientation, next_row, next_col = HCIUutils.get_direction (aFlowDir[c_row, c_col], c_row, c_col)
                adj_next_row, adj_next_col = max(0, min(next_row, (n_row-1) ) ), max(0, min(next_col, (n_col-1) ) )

                if adj_next_row != next_row or adj_next_col != next_col:  #we are going out the basin
                    riverNotReached = 0 
                    for c_listRow in simultaneous_time_update:
                        crow, ccol = c_listRow[0], c_listRow[1]
                        aAlreadyMapped[crow, ccol] = -1
                        aDestStrCellrow[crow, ccol] = -1
                        aDestStrCellcol[crow, ccol] = -1
                        
                                   
                
                elif orientation > 0:
         
                    next_UAA = aFlowAcc[next_row, next_col]

                    if next_UAA == aFlowAcc_nodata:
                        riverNotReached = 0 
                        for c_listRow in simultaneous_time_update:
                            crow, ccol = c_listRow[0], c_listRow[1]
                            aAlreadyMapped[crow, ccol] = -1
                            aDestStrCellrow[crow, ccol] = -1
                            aDestStrCellcol[crow, ccol] = -1

                    else:
                        
                        

                        c_listRow = [c_row, c_col]
                        
                        simultaneous_time_update.append(c_listRow)
                        
                        c_row, c_col = next_row, next_col # move on
                    
     
                
                else: #orientation == -1 #sink in the DEM 
                    riverNotReached = 0 
                    for c_listRow in simultaneous_time_update:
                        crow, ccol = c_listRow[0], c_listRow[1]
                        aAlreadyMapped[crow, ccol] = -1
                        aDestStrCellrow[crow, ccol] = -1
                        aDestStrCellcol[crow, ccol] = -1
                        
                        
                    
            else: #aAlreadyMapped[c_cell_indices[0], c_cell_indices[1]] == -1
                riverNotReached = 0 
                for c_listRow in simultaneous_time_update:
                    crow, ccol = c_listRow[0], c_listRow[1]
                    aAlreadyMapped[crow, ccol] = -1
                    aDestStrCellrow[crow, ccol] = -1
                    aDestStrCellcol[crow, ccol] = -1
                      
                    
                
            if aStreamNet[c_row, c_col] == 1:  # river reached
                riverNotReached = 0
                for c_listRow in simultaneous_time_update:
                    crow, ccol = c_listRow[0], c_listRow[1]
                    aAlreadyMapped[crow, ccol] = 1
                    aDestStrCellrow[crow, ccol] = c_row
                    aDestStrCellcol[crow, ccol] = c_col
               
    return aDestStrCellrow, aDestStrCellcol

     
        
def DEMstreamNet(rFlowDir, rFlowAcc, headwaters, out_dir, out_name):
    '''
    Outline the stream network in the Flow Accumulation raster, considering cells with high numbers 
    of upstream draining cells, by starting from the headwater locations provided with the "headwaters"
    argument

    Parameters
    ----------
    rFlowDir : D8 Flow Direction raster opened with the rasterio module
        The D8 Flow direction raster must follow the ESRI convention for labeling directions:
            
                :-----:-----:-----:
                : 32  : 64  : 128 :
                :-----:-----:-----:
                : 16  :cell :  1  :
                :-----:-----:-----:
                :  8  :  4  :  2  :
                :-----:-----:-----:
                    
         More info about the ESRI convention available at:      
         https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-direction.htm
            
    rFlowAcc : Flow Accumulation raster opened with the rasterio module
        
    headwaters : list
        must use the following format:
        headwaters = [(lon1, lat1), (lon2, lat2), ..., (lonN, latN)]
    out_dir : str
        location where you want to save the output stream network raster.
    out_name : str
        name of the output stream network raster. 

    Returns
    -------
    stream network raster : .tif file saved in the hard-disk
        Raster where stream network cells are labeled as 1, 0 otherwise.

    '''
    aFlowDir = rFlowDir.read(1)
    aFlowAcc = rFlowAcc.read(1)
    maxFlowAcc = np.max(aFlowAcc)
    n_row, n_col = rFlowAcc.shape
    aStreamNet = aFlowAcc.copy()
    aStreamNet[np.where(aFlowAcc>0)] = 0 # initialization 
    for headwater in headwaters:
        c_row, c_col = rFlowAcc.index(headwater[0], headwater[1])
        path_not_completed = 1
        path_indices = []
        while path_not_completed:   
            orientation, next_row, next_col = HCIUutils.get_direction (aFlowDir[c_row, c_col], c_row, c_col)
            if orientation > 0:
                aStreamNet[c_row, c_col] = 1
                if aFlowAcc[next_row, next_col] == maxFlowAcc: # just hit the outlet
                    path_not_completed = 0
                    aStreamNet[next_row, next_col] = 1
                elif aStreamNet[next_row, next_col]==1:  # just hit part of already-mapped stream net
                    path_not_completed = 0  
                else:
                    c_row, c_col = next_row, next_col
                    path_indices.append([next_row, next_col])
                    aStreamNet[c_row, c_col] = 1     
            else:   
                path_not_completed = 0 
                for pair in path_indices:
                    trow, tcol = pair[0], pair[1]
                    aStreamNet[trow, tcol] = 0
    out_meta = rFlowAcc.meta
    save_str = out_dir+out_name+'.tif'
    aRaster = rFlowAcc.read()
    aRaster[0,:,:] = aStreamNet
    with rio.open(save_str, 'w', **out_meta) as dest:
        dest.write(aRaster)  
    return 0



def AlongStrmNet_dist(rFlowDir, rFlowAcc, headwaters, cell_size, out_dir, out_name, flow_speed=1):
    '''
    Determine the distance (in meters) or the travel time to the outlet from each stream network cell

    Parameters
    ----------
    rFlowDir : D8 Flow Direction raster opened with the rasterio module
        Flow Direction labels must follow the ESRI convention:
            
                :-----:-----:-----:
                : 32  : 64  : 128 :
                :-----:-----:-----:
                : 16  :cell :  1  :
                :-----:-----:-----:
                :  8  :  4  :  2  :
                :-----:-----:-----:
                    
         More info about the ESRI convention available at:      
         https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-direction.htm
         
    rFlowAcc : Flow Accumulation raster opened with the rasterio module
        
    headwaters : list
        must use the following format:
        headwaters = [(lon1, lat1), (lon2, lat2), ..., (lonN, latN)]
    cell_size : int or float
        size (in meters) of a DEM (or Flow Direction, or Accumulation) cell.
    out_dir : str
        location where you want to save the output raster of stream network travel distances.
    out_name : str
        name of the output raster of stream network travel distances.
    flow_speed : float, optional
        set a constant travel speed if you want to get the travel time (associated with that speed) 
        instead of the travel distance to the outlet. The default is 1.

    Returns
    -------
    raster of stream network travel distances (or times) : file saved in the hard-disk
        The raster will have the travel distances/times to the outlet at each stream network cell, and 
        0s at all other basin cells.

    '''
    aFlowDir = rFlowDir.read(1)
    aFlowAcc = rFlowAcc.read(1)
    aRaster = rFlowAcc.read()
    maxFlowAcc = np.max(aFlowAcc)
    n_row, n_col = rFlowAcc.shape
    aStrTrvlT = aFlowAcc.copy()
    aStrTrvlT[np.where(aStrTrvlT>0)] = 0 # initialization
    for headwater in headwaters:
        c_row, c_col = rFlowAcc.index(headwater[0], headwater[1])
        path_not_completed = 1
        path_indices = []
        while path_not_completed:   
            orientation, next_row, next_col = HCIUutils.get_direction (aFlowDir[c_row, c_col], c_row, c_col)
            if orientation > 0:
                time_dist = cell_size*orientation / flow_speed
                for pair in path_indices:
                    trow, tcol = pair[0], pair[1]
                    aStrTrvlT[trow, tcol] = aStrTrvlT[trow, tcol] + time_dist
                if aFlowAcc[next_row, next_col] == maxFlowAcc: # just hit the outlet
                    path_not_completed = 0
                    for pair in path_indices:
                        trow, tcol = pair[0], pair[1]
                        aStrTrvlT[trow, tcol] = aStrTrvlT[trow, tcol] + cell_size/flow_speed 
                elif aStrTrvlT[next_row, next_col]>0:  # just hit part of already-mapped stream net
                    path_not_completed = 0
                    time_dist = aStrTrvlT[next_row, next_col]
                    for pair in path_indices:
                        trow, tcol = pair[0], pair[1]
                        aStrTrvlT[trow, tcol] = aStrTrvlT[trow, tcol] + time_dist # add the extra travel time to the outlet
                else:
                    c_row, c_col = next_row, next_col
                    path_indices.append([next_row, next_col])
            else:   # orientation < 0:
                path_not_completed = 0 
                for pair in path_indices:
                    trow, tcol = pair[0], pair[1]
                    aStrTrvlT[trow, tcol] = 0
                
    out_meta = rFlowAcc.meta
    save_str = out_dir+out_name+'.tif'
    aRaster[0,:,:] = aStrTrvlT
    with rio.open(save_str, 'w', **out_meta) as dest:
        dest.write(aRaster)
    return 0
        

def get_headwaters(rFlowAcc, rFlowDir, step=10000, headwaters_dir='', headwaters_file=''):
    '''
    Get headwater locations from a point shapefile, or else from a Flow Accumulation raster (dependig
    on the keyword argument provided), in the latter case by considering a fixed minimum threshold 
    (step) for the number of upstream draining cells that are necessary to identify stream network pixels. 
 
    Parameters
    ----------
    rFlowAcc : Flow Accumulation raster opened with the rasterio module
        
    rFlowDir : D8 Flow Direction raster opened with the rasterio module
        Flow Direction labels must follow the ESRI convention:
            
                :-----:-----:-----:
                : 32  : 64  : 128 :
                :-----:-----:-----:
                : 16  :cell :  1  :
                :-----:-----:-----:
                :  8  :  4  :  2  :
                :-----:-----:-----:
                    
         More info about the ESRI convention available at:      
         https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/flow-direction.htm
         
    step : int, optional
        Minimum threshold on the number of upstream draining cells to identify stream network pixels
        from the Flow Accumulation raster. The default is 10000.
    headwaters_dir : str, optional
        Location of the shapefile with headwater locations. The default is ''.
    headwaters_file : TYPE, optional
        Filename of the shapefile with headwater locations. The default is ''.

    Returns
    -------
    headwaters : list
        list of headwater locations, using the following format:
            headwaters = [(lon1, lat1), (lon2, lat2), ..., (lonN, latN)]

    '''
    
    if headwaters_file == '':
        aFlowAcc = rFlowAcc.read(1)
        aFlowDir = rFlowDir.read(1)
        aStreamNet0 = aFlowAcc.copy()
        aStreamNet0[np.where(aFlowAcc>0)] = 0  # initialization
        aStreamNet0 = np.where(aFlowAcc>=step, 1, aStreamNet0)
        aMapped = aStreamNet0.copy()
        vrow, vcol = np.where(aStreamNet0==1)
        for i in range(len(vrow)):
            row, col = vrow[i], vcol[i] 
            _, next_row, next_col = HCIUutils.get_direction (aFlowDir[row, col], row, col)
            aMapped[next_row, next_col] = np.min( (aMapped[next_row, next_col], -1) ) 
        vrow, vcol = np.where(aMapped==1)
        pheadwaters = rFlowAcc.xy( vrow, vcol )
        headwaters = [ (pheadwaters[0][i] , pheadwaters[1][i]) for i in range(len(pheadwaters[0])) ] 
    else:
        pheadwaters = gpd.read_file(headwaters_dir+headwaters_file)
        if pheadwaters.crs != rFlowAcc.crs:
            pheadwaters.to_crs(rFlowAcc.crs)
        headwaters = []
        for point in pheadwaters['geometry']:
            headwaters.append( (point.xy[0][0], point.xy[1][0]) )
        
    return headwaters





def snap_headwaters_to_highFlowAcc_cells (FlowAcc_dir, FlowAcc_file, out_dir, out_name, headwaters_dir = '', headwaters_file = '', headwaters=[], neighborhood_size=3):
    '''
    If available headwater coordinates do not exactly correspond to Flow Accumulation pixels with large numbers
    of upstream drainage cells, for instance because flow lines and DEM data come from different sources, 
    you can use this function to snap those headwater locations to neighboring Flow Accumulation 
    cells with large numbers of upstream draining pixels.

    Parameters
    ----------
    FlowAcc_dir : str
        location of Flow Accumulation raster.
    FlowAcc_file : str
        file name of Flow Accumulation raster.
    out_dir : str
        location of output shapefile with snapped headwater locations.
    out_name : str
        name of output shapefile with snapped headwater locations.
    headwaters_dir : str, optional
        location of shapefile with raw headwater coordinates. The default is ''.
    headwaters_file : str, optional
        filename of shapefile with raw headwater coordinates. The default is ''.
    headwaters : list, optional
        list of raw headwater locations, given using the following format:
            headwaters = [(lon1, lat1), (lon2, lat2), ..., (lonN, latN)]
            The default is [].
    neighborhood_size : int, optional
        for a neighborhood_size of z, the neighborhood will be a (2z+1)x(2z+1)-cell grid centered on 
        the raw headwater location. The default is 3.

    Returns
    -------
    snapped_headwaters : list 
        list of snapped coordinates, in the following format:
            snapped_headwaters = [(lon1', lat1'), (lon2', lat2'), ..., (lonN', latN')]
    shapefile of snapped headwaters : file saved in the hard-disk
        shapefile with points associated to each snapped headwater location. 

    '''
    rFlowAcc = rio.open(FlowAcc_dir+FlowAcc_file)
    snapped_headwaters = []
    IDs = []
    count=1
    
    if len(headwaters)==0 and headwaters_file != '':
        headwaters = get_headwaters(rFlowAcc, 'none', headwaters_dir=headwaters_dir, headwaters_file=headwaters_file)
          
    for head in headwaters:
        lon, lat = head[0], head[1]
        cell_vals, cell_locs, cell_indices = HCIUutils.look_neighborhood (rFlowAcc, lon, lat, neighborhood_size=neighborhood_size)
        ix = np.where(np.array(cell_vals)==np.max(np.array(cell_vals)))[0][0]
        snapped_headwaters.append(cell_locs[ix])
        IDs.append(count)
        count+=1
        
    gdf = GeoDataFrame(list(IDs), geometry=[Point(shead[0], shead[1]) for shead in snapped_headwaters], crs=rFlowAcc.crs, columns=['ID'])
    gdf.to_file(out_dir+out_name, driver='ESRI Shapefile')
    return snapped_headwaters
    
