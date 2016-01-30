# seagliderExtraction

This repository contains MATLAB routines to extract data from the seaglider missions launched by the C-MORE group at the University of Hawaii. Raw data of glider dives are available on the ftp address of SOEST at the University of Hawaii (ftp://ftp.soest.hawaii.edu/pilot/).
Each extractSg\*.m routines creates a final output file \*data.mat containing three variables:
	- dived: a structure with data about each dive (date, position, wind, sea surface height and pressure, surface O2 flux)
	- sgd: a table containing data binned on depth for each variable (T, S, etc..)
	- isod: a table containing data binned on potential density at levels defined by the average density values at the dept of the sgd variable.

Routines in the main folder require the private routines for seaglider data extraction at https://github.com/whoi-glider/glider-kit  (David Nicholson) and routines to extract OpenDap data and compute gas exchange at https://github.com/whoi-glider/oce_tools (David Nicholson and Cara Manning). The GSW Oceanographic Toolbox is also required (http://www.teos-10.org/).

Wind speed is extracted from the OPeNDAP website for the NOAA/NCDC Blended 6-hourly 0.25-degree Sea Surface Winds (http://www.ncdc.noaa.gov/thredds/dodsC/oceanwinds6hr.html).
Sea level pressure is extracted from the NCEP reanalysis 2 (http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/catalog.html).
Sea level anomaly is retrieved from the near real-time data provided by AVISO through the OPeNDAP website (http://opendap.aviso.oceanobs.com/thredds/catalog.html). 

Other Matlab functions that are called in the scripts are found in the common/ subfolder:
- betasw_ZHH2009.m: computes scattering by pure seawater (by Xiaodong Zhang)
- vdist.m: computes distances on the WGS-84 ellipsoid (by Michael Kleder)
- binning.m: binning routine (by Benedetto Barone)

The calibrations/ subfolder contains one folder for each seaglider mission with two variables:
- oxy\_cal.mat: Winkler oxygen measurements for sensor calibration
- hplc\_cal.mat: Pigment concentration measurements for sensor calibration

Benedetto Barone - Jan 2016