# Fishery management implications of demographic changes in salmon

## Summary

This repo contains code to reproduce results in Ohlberger et al. 2024 "Accounting for salmon size declines in fishery management can reduce conservation risks". We used an empirically based simulation approach to evaluate the fishery management implications of demographic changes that results from shifts in age, sex, and length compositions of mature salmon. 


## Repository Structure 
The repo contains an 'R' folder with the following sub-directories:

| Directory   | Description                                          |
| ----------- | ---------------------------------------------------- |
| `code`      | Contains code to run simulations and produce figures |
| `functions` | Contains functions needed to run the simulations     |
| `figures`   | Contains PDF copies of the main manuscript figures   |

The R folder also contains a file with scenarios (`scenarios.xlsx`) that is read when running the `run_scenarios.R` script. That script also creates a local folder called `out` to save the scenario output, which is then loaded when plotting the results figures.


## Analyses
The `code` directory contains the following files:

| File                    | Description                              |
| ----------------------- | ---------------------------------------- |
| `run_scenarios.R`       | Code to run the simulations              |
| `manuscript_figures.R`  | Code to produce the main figures         |

To reproduce the analyses, run simulations using `run_scenarios.R` and subsequently plot results figures using `manuscript_figures.R`.

## Authors 
You can contact the main author at <janohlberger@gmail.com>.
