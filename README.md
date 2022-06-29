# Argo-CbPM-NPP
Float-implementation of CbPM NPP model

## About 
information of implementation, and maybe stuff about directory setups?

## FUNCTIONS
### Main functions (to be called from script or command window):
1A compile_argo - downloads and formats float data (SPROF files) for analysis
1B extract_depth_horizons - finds and calculates parameter averages for mixed layer depth, optical depth, and euphotic depth
1C extract_sat_data - geotemporally matches satellite data with float trajectory and downloads/formats corresponding data
2A npp_depth - calculated from float profiles, produces depth-resolved and integrated NPP estimates
2B npp_surface - calculated from float or satellite data, only integrated NPP estimates

### Auxilary functions (primarily called by main functions in background):
1A BGC Argo toolbox <br />
1B calc_od_float(), calc_mld(), find_mld() <br />
1C pos2dist()
2A get_doy(), floatCbPM_JSLV3_6(), [AustinPetzold_CBpM2()]
2B get_doy(), floatCbPM_JSLV3_2_movtsdata(), AustinPetzold_CBpM2()

plotting - cmocean colormap

## REQUIREMENTS (mostly taken from One Argo README)
At least MATLAB R2019b (or any newer release) is needed to use these functions without modifications.
An Internet connection is needed to get the latest versions of Argo index and Sprof files
Memory requirements depend on the number of floats that are simultaneously loaded into memory. ???
