# paper-dfr

Manuscript and code for the paper "Discrete Forecast Reconciliation".



## Reproduce the results

Requirements: install the package DiscreteRecon

```r
devtools::install_github("https://github.com/AngelPone/DiscreteRecon")
```

### Simulation

```shell
# used to save the results
mkdir simulation/results
# Section 4.1
Rscript simulation/cross-sectional/simulate.R
Rscript simulation/cross-sectional/summary.R

# Section 4.1
Rscript simulation/temporal/simulate.R
Rscript simulation/temporal/reconcile.R
Rscript simulation/temporal/summary.R
```


### M5 empirical study

Requirements: edit the path to your original M5 dataset in `experiment/M5/data.R`

```shell
# make dir to save results
mkdir M5/results

# filter the hierarchy
Rscript experiment/M5/data.R
# run the experiments and produce figures
Rscript experiment/M5/R/run.R
```

### DC crime 

```shell
Rscript experiment/dc_crime/basef.R
Rscript experiment/dc_crime/reconcile.R
Rscript experiment/dc_crime/evaluate.R
```


## Note

- run `Rscript time.R` in each directory to produce the Table~5 in Discussion.
- Please feel free to contact the author via email or submit an issue if there are any questions.
