readme
================

## Description

Data and code associated with the study “Repeatable differences in
snowshoe hare behavior: animal personality may facilitate adaption to
climate change”. *In prep*. Please see the publication for further
details on data collection, the experimental apparatus, procedure, and
statistical analyses.

This repository will be permanently archived at Zenodo, with the latest
release to be found here:
<!--[![DOI](https://zenodo.org/badge/279408355.svg)](https://zenodo.org/badge/latestdoi/279408355).-->

A report summary containing all the R script output, including plots,
can be found at
[analysis-summary.md](https://github.com/paulqsims/snowshoe-hare_personality/blob/master/R/analysis_summary.Rmd)
in `reports/`. Markdown file versions (.md) are recommended for viewing
on github and HMTL file versions (.html) are recommended for desktop
viewing.

R Markdown files can be rendered by sourcing
[render-Rmd.R](https://github.com/paulqsims/snowshoe-hare_personality/blob/master/R/render-Rmd.R)
in `R/`

**Last update**: 2020-08-19

## Authors

Diana J.R. Lafferty (<dlaffert@nmu.edu>)

Paul Q. Sims (<paul.q.sims@gmail.com>)

L. Scott Mills (<scott.mills@mso.umt.edu>)

## Directory structure

  - `data/` Data sets used in the analyses
  - `R/` Self-contained R markdown (\*.Rmd) scripts (can be run on their
    own) and custom functions used for modifying the data, running
    analyses, and generating plots
  - `figs/` High quality figures generated for publication
  - `reports/` Markdown and html summaries of code output

## Prerequisites: R version and R packages

R version 3.6.3 (2020-02-29)

| Package\_Name | Version | Function      |
| :------------ | :------ | :------------ |
| tidyverse     | 1.3.0   | Data clean up |
| broom.mixed   | 0.2.6   | Data clean up |
| broom         | 0.7.0   | Data clean up |
| glmmTMB       | 1.0.2.1 | Analyses      |
| MuMIn         | 1.43.17 | Analyses      |
| lme4          | 1.1.23  | Analyses      |
| lmerTest      | 3.1.2   | Analyses      |
| DHARMa        | 0.3.2.0 | Analyses      |
| ggplot2       | 3.3.2   | Plots         |
| ggsignif      | 0.6.0   | Plots         |
| scales        | 1.1.1   | Plots         |
| effects       | 4.1.4   | Plots         |
| rprojroot     | 1.3.2   | Miscellaneous |
| rmarkdown     | 2.3     | Miscellaneous |
| knitr         | 1.29    | Miscellaneous |

## Metadata

Reference files:

  - [data\_activity\_LaffertyEtAl\_2020.csv](https://github.com/paulqsims/snowshoe-hare_personality/blob/master/data/data_activity_LaffertyEtAl_2020.csv)
  - [data\_day\_LaffertyEtAl\_2020.csv](https://github.com/paulqsims/snowshoe-hare_personality/blob/master/data/data_day_LaffertyEtAl_2020.csv)
  - [data\_time\_LaffertyEtAl\_2020.csv](https://github.com/paulqsims/snowshoe-hare_personality/blob/master/data/data_time_LaffertyEtAl_2020.csv)

| Column        | Description                                                                 |
| :------------ | :-------------------------------------------------------------------------- |
| ID            | Individual ID                                                               |
| trial         | Trial: 1 or 2                                                               |
| day           | Date of measurement                                                         |
| t\_point      | Time point: 1-23                                                            |
| backg         | Background color for activity data: Light or dark                           |
| loc           | Open-Field location for activity data: Inner or outer                       |
| line\_cross   | Lines crossed                                                               |
| shel\_dur     | Time spent in the shelter                                                   |
| sex           | Sex: Female or male                                                         |
| mass          | Mass in g                                                                   |
| mass\_sc      | Mass scaled by 1 standard deviation and mean centered                       |
| trial\_series | ID for each individual by trial                                             |
| t\_point\_sc  | Time point scaled by 1 standard deviation and mean centered                 |
| bkgrd         | Background color for non-activity data: Light or dark                       |
| of\_loc       | Open-Field location for non-activity data: Inner or outer                   |
| shelter       | Shelter status: Inside or outside                                           |
| exit\_lat     | Time (s) to exit the shelter                                                |
| exit\_fast    | Whether or not hares exited the shelter immediately upon release: Yes or no |
| uniq\_tiles   | Number of unique tiles entered over the course of the day                   |
| tot\_u\_tiles | Total number of unique tiles available in the open-field                    |
| tile\_fails   | Number of unique tiles hare did not enter out of the total                  |

  - Missing values are input as `NA` text

## Licenses

Software and code is licensed under the MIT License, data is licensed
under CC0, and media files are under CC-BY-4.0 - see the
[license.md](https://github.com/paulqsims/snowshoe-hare_personality/blob/master/license.md)
file for details.

## Funding information

This work was supported by North Carolina State University, University
of Montana, and the National Science Foundation Division of
Environmental Biology (grant \#1354449 to L.S.M.).
