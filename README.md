A-0-A
=========

These codes ingest A-0-A calibrated seafloor pressure ('POBS') data into Matlab. Subsequent routines apply calibrations to raw data, correct for changes in barometer pressure, and compare against output from ocean circulation models. Various tools are also included for inspecting the data for anomalies, etc.

Using on your local machine
---------------------------

Throughout, the paths are somewhat generalized (e.g., '../figures/blahblahblah.png), but in practice you'll likely get directory errors unless you manage to duplicate my unkown-to-you directory structure. So, running these codes will require you to update the paths to suit your system and organizational preferences.

Typical workflow
----------------

To replicate the analysis done for my 2025 (?) paper, run codes sequentially as:

*import_POBS.m* --> *process_POBS.m* --> *pressure_stitch.m* --> *apply_cals.m* --> *compare_w_models.m*

This will leave you with all the saved matfiles necessary to run the other bits of code

Code descriptions
----------

* /old

*Codes that were depricated as the project evolved*

* /tools <--- add to the path to ensure main codes work!

*Useful bits of code (many authored by others) that get called by the main codes*

* analyze_bar.m

*Called by process_POBS.m to determine reference pressure values from barometer data*

* analyze_cals.m

*Called by process_POBS.m to determine calibration values from APG data*

* anybaseString2num.m

*(Modified from Spahr Webb's original) Called when reading POBS encoded data*

* apply_cals.m (apply_cals_Y2.m)

*Correct pressure records for drift using pre-determined calibration series*

* barometer_stacks.m

*Plot barometric reference pressure time series for each instrument, as well as proxies calculated from temperature and ideal gas law*

* base_map.m

*Generates a simple regional map for Hikurangi upon which I can add markers for various parts of the GONDOR deployment*

* best_gauge.m (best_gauge_Y2.m)

*For cases where redundant gauges are inconsistent, makes various comparative plots to help determine which is at fault*

* compare_w_models.m

*Compares drift-corrected data against circulation models (ECCO2, GLORYS, HYCOM) and satellite altimetry*

* difference_all.m

*Makes plots for all combinations of differences between combined POBS/APG pressure network*

* find_cal.m

*Called by process_POBS.m to identify calibration intervals*

* import_POBS.m (import_POBS_Y2.m)

*Converts original encoded data to numeric and saves as matfiles*

* inspect_postcal_transient.m

*For each gauge, aligns and stacks the transients that follow the calibration intervals to assess their repeatability*

* inspect_stitch.m (inspect_stitch_Y2.m)

*Compares stitched pressure data against 1) redundant gauge, 2) tidal model, and 3) depth-matched reference station to assess the smoothness of the pressure record in the vicinity of calibration intervals (conversely, do we see static offsets after calibrations?)*

* manuscript_figs_X.m

*These codes generate the figures for my manuscript*

* model_filter_test.m

*Tool used to assess what cutoff to use for low-pass filtering of POBS and circulation model data to get a useful comparison*

* plot_tilts.m

*Determines and plots platform tilt from triaxial accelerometer data*

* pressure_stitch.m (pressure_stitch_Y2.m)

*Removes calibration intervals from time series and stitches data back together, then decimates and tidally filters*

* process_POBS.m (process_POBS_Y2.m)

*Extracts key calibration information from full data, storing it in a separate structure*

* readPOBStext.m

*(Modified from Spahr Webb's original) Called when reading POBS encoded data*

* stitch_A0A_years.m

*Starting point for merging 2022-2023 and 2023-2024 calibrated pressure data into a single, meaningful pressure record*

* trend_comps.m

*Plots of stacked pressure records before and after drift correction, as well as stem plot that demonstrates the resultant changes in linear trend*