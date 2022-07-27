# Large Lakes Statistical Water Balance Model (L2SWBM)

## Description
Meta-repository for the Large Lakes Statistical Water Balance Model (L2SWBM). The L2SWBM used multiple independent data sets to obtain the prior distributions and likelihood functions, which were then assimilated by a Bayesian framework to infer a feasible range of each water balance component.
 <img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41597-020-00613-z/MediaObjects/41597_2020_613_Fig2_HTML.png?as=webp" width=100 align="left" />

## What is L2SWBM?
Water balance models are often employed to improve understanding of drivers of change in regional hydrologic cycles. Most of these models, however, are physically-based, and few employ state-of-the-art statistical methods to reconcile measurement uncertainty and bias.

Starting in 2015, NOAA Great Lakes Environmental Research Laboratory (GLERL), along with its partners at the University of Michigan Cooperative Institute for Great Lakes Research (CIGLR), began developing a water balance model under a Bayesian Markov chain Monte Carlo framework. Through this model, we generate new estimates of monthly runoff, over-lake evaporation, over-lake precipitation, and connecting channel flows for each of the Great Lakes. The new model reconciles discrepancies between model and measurement-based estimates of each component while closing the Laurentian Great Lakes water balance.

In 2017, funding from the International Joint Commission - through their International Watersheds Initiative was received to use the model in generating a new, balanced historical (1950 - 2015) record of the Laurentian Great Lakes water balance. The project will help in resolving the regional water budget across monthly and inter-annual time scales and represents an important stepping stone towards addressing a long-standing need in the Great Lakes for clear and defensible differentiation between hydrological, climatological, geological, and anthropogenic drivers behind seasonal and long-term changes in Laurentian Great Lakes water levels.

To run these models, you will need to do the following:

Visit the JAGS Sourceforge repository and download JAGS for your computer: https://sourceforge.net/projects/mcmc-jags/files/JAGS/

Download the R statistical programming environment https://cran.r-project.org/

Open R, and enter the command "install.packages(c('rjags'))", select your mirror, and install locally if necessary

## References
Gronewold, D., J. Bruxer, D. Durnford, J. Smith, A. Clites, F. Seglenieks, S. Quian, T. Hunter, V. Fortin, 2016: Hydrological drivers of record-setting water level rise on Earth’s largest lake system. Water Resources Research, 52(5), 4026-4042. doi: 10.1002/2015WR018209. https://www.glerl.noaa.gov/pubs/fulltext/2016/20160014.pdf

Smith, J.P., & Gronewold, A.D. (2017). Development and analysis of a Bayesian water balance model for large lake systems. arXiv: Applications. https://arxiv.org/abs/1710.10161

Smith, J.P., & Gronewold, A.D. (2018). Summary Report: Development of the Large Lake Statistical Water Balance Model for Constructing a New Historical Record of the Great Lakes Water Balance. FINAL report for the International Joint Commission. https://www.glerl.noaa.gov/pubs/fulltext/2018/20180021.pdf

Smith, J.P., Gronewold, A.D., L. Read and J.L. Crooks (2019) Large Lake Statistical Water Balance Model - Laurentian Great Lakes -  1 month time window - 1980 through 2015 monthly summary data and model output. Deep Blue, University of Michigan, Ann Arbor, MI. https://doi.org/10.7302/s6h1-d521, https://deepblue.lib.umich.edu/data/concern/data_sets/2514nk609

Gronewold, A. D., Smith, J. P., Read, L., & Crooks, J. L. (2020). Reconciling the water balance of large lake systems. Advances in Water Resources, 103505. https://doi.org/10.1016/j.advwatres.2020.103505

Do, H.X., J.P. Smith, L.M. Fry and A.D. Gronewold (2020). Seventy-year long record of monthly water balance estimates for Earth’slargest lake system. Scientific Data 7, 276. https://doi.org/10.1038/s41597-020-00613-z

Do, Hong X., J.P. Smith, L.M. Fry and A.D. Gronewold (2020). Monthly water balance estimates for the Laurentian Great Lakes from 1950 to 2019 (v1.1).  Deep Blue, University of Michigan, Ann Arbor, MI. https://doi.org/10.7302/tx97-nn12

Gronewold, A. D., Do, H. X., Mei, Y., & Stow, C. A. (2021). A tug-of-war within the hydrologic cycle of a continental freshwater basin. Geophysical Research Letters, 48, e2020GL090374. https://doi.org/10.1029/2020GL090374

## What is this repository
The L2SWBM modelling system encompasses many different versions used by various agencies, along with their configuration, input, and output files. And we are continually in the process of making these versions available for general community use. This repository serves as a hub for all of our publically available versions so that they can be more easily discovered. Furthermore, these versions are in varying stages of support and this page indicates the current support level for each version.

|Tag|Description|
|------|------|
| ![](static/badges/Status-Core-green.svg) | Fully supported core component of the L2SWBM modelling system |
![](static/badges/Status-community_supported-orange.svg) | Version that is used by the community and maintained by community contributions | 
![](static/badges/Status-unsupported-red.svg) | Unsupported new version that is in the process of being generalized |

## Public repositories
|Repository|Description|Support status|
|------|------|-----------|
[Coordinating Committee operational version](https://github.com/cc-hydrosub/L2SWBM) | Operational version of L2SWBM maintained by the Coordinating Committee on Great Lakes Basic Hydraulic and Hydrologic Data | ![](static/badges/Status-Core-green.svg) |
[NOAA-GLERL development version](https://github.com/NOAA-GLERL/L2SWBM)| Development version of L2SWBM maintained by the NOAA Great Lakes Environmental Research Laboratory (GLERL), along with its partners at the University of Michigan Cooperative Institute for Great Lakes Research (CIGLR) | ![](static/badges/Status-Core-green.svg) |
[SEAS-hydro research version](https://github.com/luo-yifan/SEAS_hydro_research)| Research version of L2SWBM developed by the University of Michigan  School for Environment and Sustainability (SEAS) hydrology laboratory| ![](static/badges/Status-Core-green.svg) |
[Estimates of the water balance of the Laurentian Great Lakes using L2SWBM](https://deepblue.lib.umich.edu/data/collections/5425k9888?locale=en) | University of Michigan's Deep Blue archive contains monthly estimates of the water balance of the Laurentian Great Lakes from 1950 to 2019 that were produced by L2SWBM| ![](static/badges/Status-unsupported-red.svg) |

## Disclaimer
XXXX
