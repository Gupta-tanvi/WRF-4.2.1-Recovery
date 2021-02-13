# WRF-4.2.1-Recovery
This folder consists of the model changes and input files required for calculating recovery processes through NWP model: WRF-V4.2.1.


Following files contain the model changes:
1) module_wind_fitch.F: 
2) module_bl_mynn.F
3) module_pbl_driver.F
4) module_first_rk_step_part1.F
5) registry_EM.commmon

Following files are used for model input:
1) namelist.input
2) windturbines-ij.txt
3) wind-turbines-1.tbl

Please note that to calculate speed at 84M (turbine hub-height), use the variables U At 84M and V AT 84M. The wind-turbine-1.tbl file is over written in the code.
The input files (windturbines-ij.txt and wind-turbines-1.tbl) are sample files for WF case: inter-turbine spacing 0.5 km.


