Script to run analysis for 'Mixed effects of protected areas on avian food webs'

Author: Lucie Thompson
Computational Ecology Lab. Swansea University. UK.

Using citizen science data, we explored the effect of protection on avian food webs across a network of 45 sites of protected and unprotected communities across Europe. 
We studied how protection affects food webs, and how this effect varies with environmental conditions and management purposed. 

Version 1.1 (added sensitivity analysis)

R code to replicate the analysis:

1.DownloadGBIF.R - Script used to download GBIF data (the exact dataset can be dowloaded via the dois in Supplementary Methods)
2.EventList.eBird.R - Filters and create events for eBird
3.FilteringWDPA.R - Filtering protected area shapefiles
4.ExtractEnvrionementalVariables.R - extracting suite of environmental variables from rasters 
5.DatasetCreation - Putting everything together
6.GetContrastInFoodWebMetrics - First analysis to compute difference in Food Web metrics
7.SensitivityAnalysis.R - Equivalent to 6. but with different method for selecting protected grid cells - based on coverage % 
8.CreatePlotsAndQuantifyDrivers.R - Create all figures and run variable selection in linear regressions
9.CreatePlots-SensitivityAnalysis.R - Create supplementary figures for sensitivity analysis

Supplementary data: (intermediate results, so as to not have to re run the whole analysis)

Supplementary data 1 - main contrast table.xlsx - results from 6.GetContrastInFoodWebMetrics.R
Supplementary data 2 - Food web metrics and environmental variables.xlsx - dataframe used as input for 6.GetContrastInFoodWebMetrics.R with per grid cell information on communities, including all environmental conditions and food web metrics
Supplementary data 3 - Environmental conditions.xlsx - environmental conditions for all european grid cells

######################################################################################################################################################################

To fully re-run the analysis, you will have to download the following files: 

Biodiversity records are accessible from the GBIF at https://www.gbif.org/. dois for the GBIF downloads are available in Supplementary Methods – Table 1. 

eBird version of May 202231 for all European countries is available at http://www.ebird.org. 

Land cover data are available from <https://land.copernicus.eu/pan-european/corine-land-cover/clc2018/fetch-land-file?hash=83684d24c50f069b613e0dc8e12529b893dc172f>

The TETRA-EU 1.0 dataset of species interactions and bio-geographical region shapefile are available at <https://doi.org/10.5061/dryad.jm63xsj7b>

The World Database on Protected Areas (WDPA) March 2022 version available at <www.protectedplanet.net>

Elevation from the Shuttle Radar Topography Mission (SRTM) available at https://www.worldclim.org/data/worldclim21.html under “Elevation” 

Remoteness data available at < https://malariaatlas.org/> under “Travel time to cities”

Human density data for 2015 came from the Gridded Population of the World dataset, version 4.11 (GPWv4.11, UN WPP-Adjusted Population Density, v4.11, 2015) available at < https://ghsl.jrc.ec.europa.eu/download.php?ds=pop >

You can skip running 1.DownloadGBIF.R as this is the original code used to download GBIF records, but only by using the doi generated at the original download 
time will you have the exact same subset of occurrence records at the time of download.

Run all .R scripts one by one from 2.EventList.eBird.R to fully replicate the analysis, OR to only re-run the analysis (no data initial data filtering, extraction etc), 
use the Supplementary data 2 - Food web metrics and environmental variables.xlsx file to run 6.GetContrastInFoodWebMetrics.R and Supplementary data 1 - main contrast table.xlsx to run 8.CreatePlotsAndQuantifyDrivers.R. 

For more information, consult main manuscript and supplementary methods or contact us at:

lucie.thompson@swansea.ac.uk or 
Miguel Lurgi: miguel.lurgi@swansea.ac.uk





