# -*- coding: utf-8 -*-
"""
Author: Francesco Dell'Aira
Created on: June 14, 2024
Last Edit: June 22, 2024
Python version: 3.9.13

This script represents an asset associated with the research work by 
F. Dell'Aira and C. I. Meier (2024), titled: "Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urbanizing Watersheds".

Tutorial 2: here we show how to calculate the Hydrologic-Connectivity-based Index of Urbanization 
(HCIU; Dell'Aira and Meier, 2024), starting from the DEM and CN map of your watershed, 
using the suite of functions provided along with the article.
In this Tutorial 2 we will consider the HCIU(CN) formulation, and we will draw the stream network 
considering a minimum threshold for the number of upstream draining cells in the Flow Accumulation
raster. See Tutorial 1 (Tutorial1.py) for an example on how to consider the HCIU(n) formulation 
instead, and draw the stream network using custom headwater locations. 


REFERENCES
F. Dell'Aira and C. I. Meier, 2024."Beyond Total Impervious Area: 
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for 
Characterizing Urbanizing Watersheds".

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
This scrip use functions from the open-source Numpy and Pandas Python modules. 

https://numpy.org/
https://pandas.pydata.org/
"""
#           IMPORT SECTION
import HydroRasterAnalysis_tools as HRA
import StreamNetAnalysis_tools as SNA
import HydroConnectAnalysis_tools as HCA
import PyHCIU 
import numpy as np


#%%===========PRELIMINARY DEM PROCESSING==============

# 1) Generate sink-less DEM from raw DEM

# directory of the raw DEM (make sure you "close" the directory by adding "\\" at the end)
dem_dir = 'C:\\Users\\User1\\DEM_folder\\'
# directory where you want to save the output, sink-free DEM 
fill_dem_dir = 'C:\\Users\\User1\\SinkFreeDEM_folder\\'
# name of the raw DEM file 
dem_file = 'DEM.tif'
# name that you want to give to the output, sink-free DEM 
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
# In this Tutorial 2 we want to draw the stream network considering a minimum 
# threshold for the number of upstream draining cells in the Flow Accumulation 
# raster (traditional approach). See Tutorial 2 to see how to consider custom 
# headwater locations instead. 

# minimum threshold for the number of upstream draining cells. As an example, 
# here we consider a value of 25000, but you may need to adjust it depending 
# on your watershed. 
threshold = 25000


# We run a command that produces two rasters maps: 1) a raster with 0 and 1 cells,   
# where 1 indicates stream network cells, while 0, indicates other basin cells; and  
# 2) a raster where each stream network cell contains its travel distance to the outlet, 
# measured along the stream network (non-stream cells are labeled as 0s)

# location where to save the raster of stream network cells
StrmNet_dir = 'C:\\Users\\User1\\StreamNetwork_folder\\'
# name that you want to assign to the raster file of stream network cells
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
                step = threshold)


# name of the output StreamNetwork raster file
StrmNet_file = StrmNet_name + '.tif'

# name of the output Stream Network Travel Distances raster file
StrmNet_TravelDist_file = StrmNet_TravelDist_name + '.tif'

#%%=========GENERATE MAP OF WEIGHTS FOR THE CONNECTIVITY INDEX==========
# In this tutorial we show how to define a raster map of weights based on a CN 
# map for your basin. 
# See Tutorial 1 for an example on how to produce a raster map of weights as functions 
# of cell Manning coefficients instead. 

# location with the CN raster file
CNmap_dir = 'C:\\Users\\User1\\CNmap_folder\\'
# file name of the CN raster
CNmap_file = 'CNmap.img'

# We generate the map of weights W, using the formulation
# W(CN) = CN/100 (Dell'Aira and Meier, 2024)
# The function will generate a W map resampled to match the resolution of the
# Flow Accumulation raster. 

# location where you want to save the output W raster map
W_dir = 'C:\\Users\\User1\\WeightsMap_folder\\'
# name that you want to assign to the output W raster map 
W_name = 'W_CN'

# set Wflag = 2 to calculate W(CN)=CN/100.
# See Tutorial 1 for an example with W(n). 
# Also note that you can optionally set your own custom function inside   
# "Wraster_f", by following the instructions in the open-source "Wraster_f"
# script. 

Wflag = 2 # for W(CN)= CN/100
 

HCA.Wraster_f (CNmap_dir, CNmap_file, 
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
HCI_name = 'HCI_CN'

# calculate the HCI raster
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

# name that you want to assign to the raster of normalized connectivity indices
HCI_hat_name = 'normHCI_CN'

# weight for the fully impervious benchmark
W_imp = 0.99 # for the CN-based formulation HCIU(CN)

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


print('The value of HCIU(CN) for your basin is ' + str(np.round(HCIU, 3)) )

