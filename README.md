# Estimating the epidemic reproduction number from temporally aggregated incidence data: a statistical modelling approach and software tool

Rebecca K Nash<sup>1</sup>, Samir Bhatt<sup>1,2</sup>, Anne Cori<sup>1</sup>, Pierre Nouvellet<sup>1,3</sup>

<sup>1</sup> MRC Centre for Global Infectious Disease Analysis, Imperial College London, UK <br>
<sup>2</sup> Section of Epidemiology, Department of Public Health, University of Copenhagen <br>
<sup>3</sup> School of Life Sciences, University of Sussex, UK

This repository contains all code necessary to reproduce the analysis and figures
as described in Nash et al 2023.

## Software dependencies

You will need the latest version of EpiEstim to run the analysis. Please install
using the following:

```bash
install.packages('EpiEstim', repos = c('https://mrc-ide.r-universe.dev',
'https://cloud.r-project.org'))
```

## How to run

```bash
setwd("~/EM_EpiEstim_Nash2023")
```

### Data used

The 'real_data' folder contains all of the real-world data (US influenza and UK 
COVID-19) used in the main paper.

The 'supplementary_analysis' folder contains the simulated data used in the
supplementary analysis.

### Code

**Within the 'main_analysis' folder you will find all code necessary to recreate** 
**the influenza and COVID-19 figures in the main paper:**

* flu.R 
* covid_cases.R 
* covid_deaths.R 

**Within the 'supplementary_analysis' folder you will find all code necessary to** 
**recreate the analysis in the appendix.**

**Supplementary analysis for real-world data:**

* flu_appendix.R
* covid_cases_appendix.R
* covid_deaths_appendix.R
* incidence_weekday.R

**Simulation study:**

* constant_rt.R
* time_varying_rt_sudden.R
* time_varying_rt_gradual.R
* weekend_effects.R
* different_aggregations_aligned.R
* different_aggregations_misaligned.R
* mid_aggregation_variations.R
* loess_smoothing.R

### Vignettes

For other worked examples, FAQs, and more details about how to apply EpiEstim 
to temporally aggregated data, please see the following vignette available in 
the EpiEstim R package:

* https://mrc-ide.github.io/EpiEstim/articles/EpiEstim_aggregated_data.html

For a breakdown of how the EM algorithm works "under the hood", please see the 
vignette available in this repository:

* em_explanation.Rmd
