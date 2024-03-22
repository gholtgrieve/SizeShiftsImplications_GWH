# Fishery management implications of demographic changes in salmon

This repo contains code to reproduce results in Ohlberger et al. 2024 "Accounting for salmon size declines in fishery management can reduce conservation risks". 

## Summary
We used an empirically based simulation approach (similar to a Management Strategy Evaluation) to assess the fishery management implications of demographic changes that result from shifts in age, sex, and length compositions of mature salmon. We assumed that the harvest control rule used in fishery management was based on a stock-recruit analysis that estimates the stock size that produces maximum sustainable yield, as typically done for Pacific salmon. We then evaluated the performance of the fishery management when the stock-recruit analyses accounted for observed demographic trends by converting spawner abundance to total egg mass and compared that to a traditional model using spawner abundance as a metric of stock size. Performance was assesses by calculating the expected long-term average harvest and return, as well as metrics of conservation risk, such as the probability that returns or spawner abundance would drop below certain abundance thresholds. 


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
Please contact the main author at <janohlberger@gmail.com>.
