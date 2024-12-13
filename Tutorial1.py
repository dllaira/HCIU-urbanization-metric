# -*- coding: utf-8 -*-
"""
Author: Francesco Dell'Aira
Created on: June 14 2024
Last Edit: June 22, 2024
Python version: 3.9.13

This script represents an asset associated with the research work by 
F. Dell'Aira and C. I. Meier (2024), titled: "Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urban Watersheds".

Tutorial 1: here we show how to calculate the Hydrologic-Connectivity-based Index of Urbanization 
(HCIU; Dell'Aira and Meier, 2024), starting from the DEM and LULC map of your watershed, 
using the suite of functions provided along with the article.
In this Tutorial 1 we will consider the HCIU(n) formulation, and we will draw the stream network 
from custom headwater locations. See Tutorial 2 (Tutorial2.py) for an example on how to consider 
the HCIU(CN) formulation instead, and draw the stream network using the traditional approach, i.e., 
adopting a fixed minimum threshold for the number of upstream draining cells in the Flow 
Accumulation Raster. The latter approach is recommended if you do not have information about 
the headwater locations for your basin. 


REFERENCES
F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urban Watersheds". Hydrology and Earth System Sciences.

DEPENDENCIES: 
1) Numpy Python library (https://numpy.org/)
2) Pandas Python library (https://pandas.pydata.org/)
3) StreamNetAnalysis_tools module ( StreamNetAnalysis_tools.py file provided along with this script)
4) HydroConnectAnalysis_tools module ( HydroConnectAnalysis_tools.py file provided along with this script)
5) HydroRasterAnalysis_tools module (HydroRasterAnalysis_tools.py file provided along with this script)
6) HCIUutils module (HCIUutils.py file provided along with this script)
7) PyHCIU module (PyHCIU.py file provided along with this script)
    
StreamNetAnalysis_tools, HydroConnectAnalysis_tools, HydroRasterAnalysis_tools, HCIUutils, and 
PyHCIU modules have their own dependencies. See details in their corresponding .py files. 

                        
CREDITS:
This scrip use functions from the open-source Numpy and Pandas Python libraries. 

https://numpy.org/
https://pandas.pydata.org/
"""
#           IMPORT SECTION
import HCIUutils as UTLS
import HydroRasterAnalysis_tools as HRA
import StreamNetAnalysis_tools as SNA
import HydroConnectAnalysis_tools as HCA
import PyHCIU 
import pandas as pd
import numpy as np


#%%===========PRELIMINARY DEM PROCESSING==============

# 1) Generate sink-less DEM from raw basin DEM

# directory of the raw basin DEM (make sure you "close" the directory by adding "\\" at the end)
dem_dir = 'C:\\Users\\User1\\DEM_folder\\' #example
# directory where you want to save the output, sink-free DEM 
fill_dem_dir = 'C:\\Users\\User1\\SinkFreeDEM_folder\\'
# name of the raw DEM file 
dem_file = 'DEM.tif'
# name that you want to give to the output, sink-free basin DEM 
# (no need to include the extension, which will automatically be ".tif")
fill_dem_name = 'fDEM'

# command to generate sink-less DEM
HRA.fill_dem_f (dem_dir, dem_file, 
                fill_dem_dir, fill_dem_name)

# 2) Generate Flow Direction and Flow Accumulation raster maps 

# name of the sink-free DEM file previously generated 
fill_dem_file = fill_dem_name + '.tif'

# directory where you want to save the output Flow Direction raster
flowDir_dir = 'C:\\Users\\User1\\FlowDir_folder\\'

# name that you want to assign to the output Flow Direction raster
# (no need to include the extension, which will automatically be ".tif")
flowDir_name = 'FlowDir'

# command to generate the Flow Direction raster map using the D8 algorithm
# the output Flow Direction raster follows the ESRI convention for labeling flow directions
HRA.D8FlowDir_f (fill_dem_dir, fill_dem_file, 
                 flowDir_dir, flowDir_name)

# name of the Flow Direction raster file previously generated 
flowDir_file = flowDir_name + '.tif'

# directory where you want to save the output Flow Accumulation raster
flowAcc_dir = 'C:\\Users\\User1\\FlowAcc_folder\\'

# name that you want to assign to the output Flow Accumulation raster
# (no need to include the extension, which will automatically be ".tif")
flowAcc_name = 'FlowAcc' 

# command to generate the Flow Accumulation raster map from the D8 Flow Direction raster map
HRA.D8FlowAcc_f(flowDir_dir, flowDir_file, 
                flowAcc_dir, flowAcc_name)

# name of the Flow Accumulation raster file previously generated 
flowAcc_file = flowAcc_name + '.tif'


# 3) Generate raster of Slopes

# directory where you want to save the output raster of Slopes 
S_dir = 'C:\\Users\\User1\\Slopes_folder\\'

# name that you want to assign to the output raster of Slopes
S_name = 'Slopes' 

HRA.Slope_f (fill_dem_dir, fill_dem_file, 
             S_dir, S_name)

# file name of the output raster of Slopes
S_file = S_name + '.tif'


#%%===========OUTLINE THE STREAM NETWORK==============
# In this Tutorial 1 we want to use the headwaters from official blue 
# lines to outline the stream network on the Flow Accumulation raster. 
# See Tutorial 2 for an example of how to draw flow lines based on a fixed 
# minimum threshold for the number of upstream draining cells instead (traditional 
# approach).

# Below we show how to feed the headwater locations from an input shapefile; however, 
# if you wish, you can alternatively pass headwater locations manually, simply providing  the
# "headwaters" parameter (instead of "headwaters_dir" and "headwaters_file") to the
# "snap_headwaters_to_highFlowAcc_cells" function, using the following format for defining
# that parameter and for the call to the function:
#
#   headwaters = [(lon1, lat1), (lon2, lat2), ..., (lonN, latN)]    
#  
#   SNA.snap_headwaters_to_highFlowAcc_cells(flowAcc_dir, flowAcc_file,
#                                            snapped_headwaters_dir, snapped_headwaters_file, 
#                                            headwaters = headwaters,
#                                            neighborhood_size=3)


# location of the shapefile of the points with basin headwater locations 
headwaters_dir = 'C:\\Users\\User1\\headwaters_folder\\'
# name of the shapefile of the points with basin headwater locations
headwaters_file = 'headwaters.shp'

# It may be the case that the flow lines (and associated headwaters) do not perfectly 
# match the emerging network of Flow Accumulation cells with high numbers of draining cells, 
# for example if the DEM and the stream network data come from different sources. 
# To account for this, it may be preferable to first "snap" the available headwater locations 
# to nearby high-value Flow Accumulation cells. To do so, we preliminarly use the following  
# command, by setting the size of the area around the raw headwater (i.e., its "neighborhood") 
# within which you want to search for the Flow Accumulation cell with large value (i.e., the 
# corresponding headwater in the Flow Accumulation raster). If we use a value of the
# "neighborhood_size" parameter of, say, 3, the function below will search within a (3+1)-by-(3+1)-cell
# area centered on the raw headwater location. It is suggested choosing the value of "neighborhood_size" 
# parameter based on the resolution of your DEM. 

# location where you want to save the snapped headwaters
snapped_headwaters_dir = 'C:\\Users\\User1\\snapped_headwaters_folder\\'
snapped_headwaters_file = 'snapped_headwaters.shp'
SNA.snap_headwaters_to_highFlowAcc_cells(flowAcc_dir, flowAcc_file,
                                         snapped_headwaters_dir, snapped_headwaters_file, 
                                         headwaters_dir = headwaters_dir, headwaters_file = headwaters_file,
                                         neighborhood_size=3)

# Next, we run a command that produces two rasters maps: 1) a raster with 0 and 1 cells,   
# where 1 indicates stream network cells, while 0, indicates other basin cells; and  
# 2) a raster where each stream network cell contains its travel distance to the outlet, 
# measured along the stream network (non-stream cells are labeled as 0s)


# location where to save the raster of stream network cells
StrmNet_dir = 'C:\\Users\\User1\\StreamNetwork_folder\\'
# name that you want to assign to the raster file of stream network cells for the basin
StrmNet_name = 'StreamNetwork_raster'
# location where to save the raster of stream network travel distances
StrmNet_TravelDist_dir = 'C:\\Users\\User1\\StreamNetwork_folder\\'
# name that you want to assign to the raster file of stream network travel distances
StrmNet_TravelDist_name = 'StreanNet_travelDistance'

# in order to calculate the along-stream-network travel distance, the function also needs
# the size of cells (in meters) in the Flow Accumulation raster. For a 1/3 arc second DEM, 
# the size of the cells is 10 m

cell_size = 10 # meters

SNA.StreamNet_f(flowDir_dir, flowDir_file, 
                flowAcc_dir, flowAcc_file, 
                StrmNet_dir, StrmNet_name,
                StrmNet_TravelDist_dir, StrmNet_TravelDist_name, 
                cell_size, 
                headwaters_dir = snapped_headwaters_dir, 
                headwaters_file = snapped_headwaters_file)


# name of the output StreamNetwork raster file
StrmNet_file = StrmNet_name + '.tif'

# name of the output Stream Network Travel Distances raster file
StrmNet_TravelDist_file = StrmNet_TravelDist_name + '.tif'

#%%=========GENERATE MAP OF WEIGHTS FOR THE CONNECTIVITY INDEX==========
# In this tutorial we show how to define a raster map of weights based on Manning's 
# roughness coefficients n associated with different land-use/land-cover (LULC) types. 
# See Tutorial 2 for an example on how to produce a raster map of weights as functions 
# of cell Curve Number values instead. 

# The function used below to associate Manning's coefficients to different LULC categories 
# needs a conversion table. As an example, here we define a conversion table for NLCD
# LULC categories (as per the work by Dell'Aira and Meier, 2024). Note that you will need to 
# adapt the conversion table to your own LULC dataset. 

LULCs = [11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 44, 45, 52, 71, 81, 82, 90, 95]
ns = [.035, .05, .12, .1, .07, .02, .1, .4, .4, .55, .6, .8, .4, .3, .3, .35, .5, .5]

LULC_col = 'LULC'
Manning_col = 'n'
convtab = pd.DataFrame({LULC_col: LULCs, Manning_col: ns})

# location of LULC raster map of the basin
LULC_dir = 'C:\\Users\\User1\\LULCmap_folder\\'
# file name of LULC raster map of the basin
LULC_file = 'LULCmap.img'

# location where you want to save the output raster map of Manning coefficients
ManningMap_dir = 'C:\\Users\\User1\\Manning_folder\\'
# name that you want to assign to the raster map of Manning coefficients
ManningMap_name = 'ManningMap'


UTLS.LULC_to_Manning(LULC_dir, LULC_file, 
                     ManningMap_dir, ManningMap_name, 
                     convtab, LULC_col, Manning_col)

# file name of the output map of Manning's roughness coefficients
ManningMap_file = ManningMap_name + '.tif'

# Now we can generate the map of weights W for the basin, using the formulation
# W(n) = 1-n (Hooke et al., 2021)
# The function will generate a W map resampled to match the resolution of the
# Flow Accumulation raster. 

# location where you want to save the output W raster map
W_dir = 'C:\\Users\\User1\\WeightsMap_folder\\'
# name that you want to assign to the output W raster map 
W_name = 'W_n'

# set Wflag = 1 to calculate W(n)=1-n.
# See Tutorial 2 for an example with W(CN). 
# Also note that you can optionally set your own custom function inside   
# "Wraster_f", by following the instructions in the open-source "Wraster_f"
# script.

Wflag = 1 # for W(n)= 1-n
 

HCA.Wraster_f (ManningMap_dir, ManningMap_file, 
               flowAcc_dir, flowAcc_file, 
               W_dir, W_name, Wflag)

# file name of the output raster of W
W_file = W_name + '.tif'

#%%===========CALCULATE RASTER MAP OF CONNECTIVITY INDICES==============
# In this section we calculate the raster of (absolute) hydrologic connectivity 
# indices HCI (Dell'Aira and Meier 2024)

# location where you want to save the HCI raster 
HCI_dir = 'C:\\Users\\User1\\HCI_folder\\'
# file name for the output HCI raster 
HCI_name = 'HCI_n'

# calculate the HCI raster for the basin
HCA.HCI_f(S_dir, S_file, 
          W_dir, W_file, 
          flowDir_dir, flowDir_file, 
          flowAcc_dir, flowAcc_file, 
          StrmNet_dir, StrmNet_file, 
          HCI_dir, HCI_name, cell_size)

#%%===========CALCULATE HCIU==============
# In this section we calculate the normalized connectivity index and HCIU for the watershed

# location where you want to save the raster of normalized connectivity indices
HCI_hat_dir = 'C:\\Users\\User1\\normHCI_folder\\'

# name that you want to assign to the raster of normalized connectivity indices for the basin
HCI_hat_name = 'normHCI_n'

# weight for the fully impervious benchmark
W_imp = 0.98 # for the n-based formulation HCIU(n)

# calculate the raster of normalized HCI
PyHCIU.normHCI_f(flowDir_dir, flowDir_file, 
                 flowAcc_dir, flowAcc_file, 
                 W_dir, W_file,
                 S_dir, S_file, 
                 StrmNet_dir, StrmNet_file, 
                 HCI_hat_dir, HCI_hat_name, 
                 cell_size, W_imp)

# file name of the output raster of normalized HCI
HCI_hat_file = HCI_hat_name + '.tif'


# calculate HCIU
HCIU = PyHCIU.HCIU_f(flowDir_dir, flowDir_file, 
                     flowAcc_dir, flowAcc_file,
                     HCI_hat_dir, HCI_hat_file, 
                     StrmNet_dir, StrmNet_file,
                     StrmNet_TravelDist_dir, StrmNet_TravelDist_file, 
                     cell_size)


print('The value of HCIU(n) for your basin is ' + str(np.round(HCIU, 3)) )

