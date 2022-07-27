# WRF-Hydro <img src="https://github.com/NCAR/wrf_hydro_nwm_public/blob/master/.github/images/wrf_hydro_symbol_logo_2017_09.png" width=100 align="left" />

## Description
Meta-repository for the WRF-Hydro modeling ecosystem. For additional information on WRF-Hydro please visit our official [website](https://ral.ucar.edu/projects/wrf_hydro/overview).

## What is WRF-Hydro?
The WRF-Hydro modeling system was originally designed as a model coupling framework to facilitate easier coupling between the Weather Research and Forecasting model and components of terrestrial hydrological models.  WRF-Hydro is both a stand-alone hydrological modeling architecture as well as a coupling architecture for coupling of hydrological models with atmospheric models.  WRF-Hydro is fully-parallelized to enable its usage on clusters and high performance computing systems alike.

Like the WRF model it does not attempt to prescribe a particular or singular suite of physics but, instead, is designed to be extensible to new hydrological parameterizations.  Although it was originally designed to be used within the WRF model, it has evolved over time to possess many additional attributes as follows:

The architecture is intended to significantly simplify the often laborious task of integrating, or coupling, existing and emerging hydrological models into the WRF modeling framework.  In doing so, an extensible, portable and scalable environment for hypothesis testing, sensitivity analysis, data assimilation and environmental prediction has emerged. 

The WRF-Hydro system has adopted a ‘community-based’ development processes with an open and participatory working group environment.  NCAR in collaboration with other NSF and university entities are developing a support structure for WRF-Hydro in the way of model documentation, public, online code repositories, test cases and many pre- and post-processing utilities.

## What is this repository
The WRF-Hydro modelling system encompasses many ancillary tools, and we are continually in the
process of making these tools available for general community use. This repository serves as a hub
for all of our publically available tools so that they can be more easily discovered. Furthermore,
these tools are in varying stages of support and this page indicates the current support level for
each tool set.

|Tag|Description|
|------|------|
| ![](static/badges/Status-Core-green.svg) | Fully supported core component of the WRF-Hydro modelling system |
![](static/badges/Status-community_supported-orange.svg) | Toolset that is used by the community and maintained by community contributions | 
![](static/badges/Status-unsupported-red.svg) | Unsupported new toolset that is in the process of being generalized |



## Official repositories
|Repository|Description|Support status|Build status|
|------|------|-----------|---------|
[wrf_hydro_nwm_public](https://github.com/NCAR/wrf_hydro_nwm_public) | Community repository for the core WRF-Hydro model code | ![](static/badges/Status-Core-green.svg) | [![Build Status](https://travis-ci.org/NCAR/wrf_hydro_nwm_public.svg?branch=master)](https://travis-ci.org/NCAR/wrf_hydro_nwm_public) |
[wrf_hydro_docker](https://github.com/NCAR/wrf_hydro_docker) | Docker files to support the WRF-Hydro modelling system | ![](static/badges/Status-Core-green.svg) | [![Build Status](https://travis-ci.org/NCAR/wrf_hydro_docker.svg?branch=master)](https://travis-ci.org/NCAR/wrf_hydro_docker) | 
[wrf_hydro_training](https://github.com/NCAR/wrf_hydro_training) | Materials from WRF-Hydro training workshops | ![](static/badges/Status-Core-green.svg) | None | 
[wrf_hydro_py](https://github.com/NCAR/wrf_hydro_py) | Python API for the WRF-Hydro model | ![](static/badges/Status-unsupported-red.svg) | [![Build Status](https://travis-ci.org/NCAR/wrf_hydro_py.svg?branch=master)](https://travis-ci.org/NCAR/wrf_hydro_py) | 
[rwrfhydro](https://github.com/NCAR/rwrfhydro) | A toolbox for WRF-Hydro input and output files in R | ![](static/badges/Status-community_supported-orange.svg) | None |

