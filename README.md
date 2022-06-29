# Argo-CbPM-NPP
Float-implementation of CbPM NPP model

## About 
information of implementation, and maybe stuff about directory setups?

## FUNCTIONS
### Main functions (to be called from script or command window):
1A  compile_argo.m - downloads and formats float data (SPROF files) for analysis<br/>
1B  extract_depth_horizons.m - finds and calculates parameter averages for mixed layer depth, optical depth, and euphotic depth<br/>
1C  extract_sat_data.m - geotemporally matches satellite data with float trajectory and downloads/formats corresponding data<br/>
2A  npp_depth.m - calculated from float profiles, produces depth-resolved and integrated NPP estimates<br/>
2B  npp_surface.m - calculated from float or satellite data, only integrated NPP estimates<br/>
    requires HYCOM and ZNO3 data from 

### Auxilary functions (primarily called by main functions in background):
1A  [BGC-Argo-Mat](https://github.com/NOAA-PMEL/BGC_Argo_Mat_Toolbox) toolbox<br/>
1B  calc_od_float.m, calc_mld.m, find_mld.m<br/>
1C  pos2dist.m<br/>
2A  get_doy.m, floatCbPM_JSLV3_6.m, AustinPetzold_CBpM2.m<br/>
2B  get_doy.m, floatCbPM_JSLV3_2_movtsdata.m, AustinPetzold_CBpM2.m<br/>
<br/>
plotting - cmocean.m colormap<br/>

## REQUIREMENTS (mostly taken from One Argo README)
At least MATLAB R2019b (or any newer release) is needed to use these functions without modifications.
An Internet connection is needed to get the latest versions of Argo index and Sprof files
Memory requirements depend on the number of floats that are simultaneously loaded into memory. ???
