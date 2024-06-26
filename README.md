# HCIU: An Urbanization Metric Sensitive to the Spatial Arrangement of Urbanized Sectors with Respect to Basin Topographic Structure and the Heterogeneous Mosaic of Other Land-use/Land-cover Patches.  
Accompanying code for our research article "Beyond Total Impervious Area: A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for Characterizing Urbanizing Watersheds".

    
```
Dell'Aira, F., and Meier, C. I. (submitted). Beyond Total Impervious Area:
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for Characterizing
Urbanizing Watersheds. Submitted to: Hydrol. Earth Syst. Sci.
```

This suite of functions allows to perform the analyses shown in our work for our case-study basins and to extend the application of HCIU to other case studies. 

## Content of the repository 
- `HydroRasterAnalysis_tools.py` Suite of functions to perform required DEM preprocessing operations, such as deriving the flow direction and flow accumulationn raster maps. Output raster maps are compatible with any GIS software and geoscripting tools.  
- `StreamNetAnalysis_tools.py` Module to identify the stream network of a watershed starting from its flow accumulation raster map. It allows to consider either a fixed minimum threshold for the number of upstream draingage cells, or else custom headwater locations. Output geospatial data are compatible with any GIS software and geoscripting tools.  
- `HydroConnectAnalysis_tools.py` Group of functions to calculate the hydrologic connectivity index proposed in our work, which is a modification of the formulation originally proposed by Borselli et al. (2008) and widely used in the literature (see references in our article). Output raster maps are compatible with any GIS software and geoscripting tools.  
- `PyHCIU.py` Module to calculate HCIU and the raster map of normalized connectivity indices, as defined in our work. Output raster maps are compatible with any GIS software and geoscripting tools.
- `HCIUfit_test_tools.py` Group of functions to fit regional peak-flow equations using competing urbanization metrics (HCIU and TIA in our work) and compare their performances, either considering the full dataset or else in a k-fold validation framework. 
- `HCIUplot_tools.py` Functions to plot the results. 
- `HCIUutils.py` Miscellanea of other functions that support our analysis.
- `Tutorial1.py` A script demonstrating function streamlining for the calculations of HCIU, starting from the DEM and LULC map of a watershed, as well as the knowledge of its headwater locations. In this tutorial, the Manning-based formulation HCIU(n) is considered (see our work for more details).
- `Tutorial2.py` Like Tutorial 1, but considering HCIU(CN) instead. It also illustrates how to draw the stream network using the traditional approach of setting a minimum threshold for the number of upstream cells in the flow accumulation raster (as an alternative to using headwater locations).
- `Tutorial3.py` It streamlines the procedures for fitting a regional model and comparing the explanatory power achieved with HCIU and the traditional fraction of total impervious area (TIA).

## How to run the code locally
<p>Download this repository to run our code on your machine. As we used some external open-source libraries to create our code, before running it you should first install the required external libraries in your working environment. Each module in this repository includes information about its specific dependencies (including a link to the libraries' documentation), and will automatically import the necessary libraries when executed. Also download the necessary raster maps and geospatial data from the links provided in the Code and Data availability statement in our article. <p>
Follow Tutorial 1 to calculate HCIU for the case-study basins, and Tutorial 2 for doing the same but using the alternative HCIU(CN) formulation. Then, test the predictive power of HCIU(n) and HCIU(CN) against the traditional fraction of total impervious area (TIA) by following the steps illustrated in Tutorial 3. You can also directly run Tutorial 3 using the basin information and calculated HCIU values reported in Table A1 of our article (table is also provided in this repository, in the Data folder). All tutorials include detailed information across all necessary steps. Additional documentation for each function is provided in the modules in this repository. 
</p>

## Citation
If you use any of this code in your experiments, please make sure to cite the following publication: 

Dell'Aira, F., and Meier, C. I. (submitted). Beyond Total Impervious Area:
A New Lumped Descriptor of Basin-Wide Hydrologic Connectivity for Characterizing
Urbanizing Watersheds. Submitted to: Hydrol. Earth Syst. Sci.

## License
Apache 2.0





