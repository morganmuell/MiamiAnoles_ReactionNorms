# Analyzing reaction norm evolution across seven species of *Anolis* lizards
Code for analyzing reaction norm evolution in a community of invasive *Anolis* lizards in Miami, FL, USA.

Citation:
Muell MR, Hall JM, Smith KV, Oaks JR, Wolak ME, and DA Warner. In press. Comparison of thermal developmental plasticity among seven recently sympatric Anolis species: insights into the evolution of reaction norms. Evolution.

## Guide to Repository:

`'datasets'` folder contains phenotype data in form of R objects, phylogeny used in analysis, and occurrence records and climate layers used to draw boundaries for climate data used in analysis.

`'scripts'` folder contains the following:
  - Environmental data script, where climate layers are trimmed down to species-specific ranges and values are pulled to create climate fixed effects for reaction norm modeling.
  - Modeling scripts, optimized for the Easley HPC cluster at Auburn University. Within each script, Bayesian mixed models extract reaction norms and test for the influence of common ancestry versus climate on reaction norm evolution among species. 
  - Two plotting scripts ("plots_ch1.R" and "plot_utils.R"). The main plotting script contains code used to generate Figs 2 and 3, and the heat map from Fig 1. plot_utils.R contains support functions for plots in main script.

`'results'` folder is divided by trait and contains output summary statistics and model diagnostics from the modeling scripts. "run1" and "run2" are results from identical runs of the same script, except with the random starting seed changed before each model run.

## GBIF Citations for Occurrence Data:

*Anolis carolinensis*:
GBIF.org (22 March 2023) GBIF Occurrence Download https://doi.org/10.15468/dl.cbyj5q

*Anolis chlorocyanus*:
GBIF.org (21 March 2023) GBIF Occurrence Download https://doi.org/10.15468/dl.nke443

*Anolis cristatellus*:
GBIF.org (21 March 2023) GBIF Occurrence Download https://doi.org/10.15468/dl.mjhf9x

*Anolis cybotes*:
GBIF.org (21 March 2023) GBIF Occurrence Download https://doi.org/10.15468/dl.uhsr58

*Anolis distichus*:
GBIF.org (21 March 2023) GBIF Occurrence Download https://doi.org/10.15468/dl.hucuww

*Anolis equestris*:
GBIF.org (22 March 2023) GBIF Occurrence Download https://doi.org/10.15468/dl.hj473r

*Anolis sagrei*:
GBIF.org (22 March 2023) GBIF Occurrence Download https://doi.org/10.15468/dl.2ku9ak

