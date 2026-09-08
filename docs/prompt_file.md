# We are updating the time series with 2026 data.
- Find the raw data files for 2026 in the data folder.
- Use the file add.cpue.R in the scripts folder as the template for QA/QC for 2026 data, and adding 2026 to the complete time series.
    - Read and understand the routines applied in add.cpue.R.
    - Append the R script needed to process the 2026 data to the end of add.cpue.R.
    - Run the updated add.cpue.R script to process the 2026 data and update the time series.
- Create an error-checking routine focused on the 2026 data, but also considering all years, and present the plan for approval before proceeding.

# QA/QC
- For the repeat sampling of Cook Bay sites in 2009, change the day for the second set (the one with the later time) to Julian day 197 from Julian day 196.
- - The two sets at Anton Larsen Bay on Julian day 195 in 2012 are accidental repeats with identical catches, remove the second one.
- The timestamp for the repeat Laminaria East sampling on day 239 in 2010 indicates that the second sampling event was a mistaken site entry - it should be Middle Cove, not Laminaria East. Change the 1215 set to reflect this. Confirm that the Middle Cove name fits with the naaming convention (is consistent with the name used in other years).
- For the 2023 Anton Larsen Bay, clarifying information is being sought from the data provider.
- For the 20214 Fox-1 problem: 589 was a bad tow and should be dropped (use for CPUE set to 'n').

# Additional QA/QC
- Also make these changes in the R script.
- For the 2023 Anton Larsen Bay discrepancy, the set that caught 3 age-0 and 5 age-1 Pacific cod should be changed to site = Eelgrass Point.
- For the re-use of station number 531 in 2023, for site = Jap-5, change the station number to 529.

# Stations / sites to consider dropping
- This section is meant to confirm that sites with little or no replication are dropped from analysis.
- Confirm that the following bays are dropped from the CPUE time series: Canton Harbor, Sanak, Chief Cove, Kujulik.
- Also confirm that the following stations were at sites that were only sampled in 2018 or 2024 and if so, drop them from the CPUE time series: station # 57, 122, 18, 546, 75, 1, 127.

# Stations to drop #2
- Drop the following sites from the CPUE time series: Ugak-1, Sand-1, Kilu-9, Kai-5, BB-5.
- Check that there are no remining sites in the time series that are sampled only once or twice.

# Stations to drop #3
- Also drop Kilu-7 and Kilu-8 from the CPUE time series. 
- Also drop BB-1 from the CPUE time series.

# Run cohort strength models
- Run cohort_strength.R to estimate year class strength for 2006-2026 for both cod and pollock (age-0 only).
- Confirm that the .csv and plots of the updated time series will be updated in the repo after the models run.
- Provide estimates of tun time remaining for each model; update these occasionally.