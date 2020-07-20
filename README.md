# covid19-release
A repository used for modelling COVID-19 transmission in the province of Québec.

# Impact of the semaine de relâche

[Link to pre-print](https://www.medrxiv.org/) only leads to medrxiv for now.

The Canadian epidemics of COVID-19 exhibit distinct early trajectories, with Québec bearing a very high initial burden. The semaine de relâche, or March break, took place two weeks earlier in Québec as compared to the rest of Canada. This event may have played a role in the spread of severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2). We aimed to examine the role of case importation in the early transmission dynamics of SARS-CoV-2 in Québec.

Using detailed surveillance data, we developed and calibrated a deterministic SEIR-type compartmental model of SARS-CoV-2 transmission. We explored the impact of altering the number of imported cases on hospitalizations. Specifically, we investigated scenarios without case importation after March break, and as scenarios where cases were imported with the same frequency/timing  as neighboring Ontario.

## Installation

Downloading this repository.

```bash
git clone https://github.com/pop-health-mod/covid19-release
```

### A note on data availability

The local data used to parametrize and calibrate the model are currently not in the public domain. Code to reproduce the analyses is still made available to the public to ensure the greatest transparency possible regarding our methodology.

Data used for comparison with other provinces and selected countries can be retrieved from the **COVID-19 Data Repository by the Center for Systems Science and Engineering (CSSE) at Johns Hopkins University**. Since the time series are updated daily, pulling frequently is advised.

```bash
git clone https://github.com/CSSEGISandData/COVID-19.git
```

## Calibration and output

The model is first calibrated to local data on hospitalizations. We then produce a range of scenarios regarding the impact of Québec's March break.

```R
# Script for calibration
source run_model_quebec_relache_rw.R

# Script for model outputs
source relache.R
```

## License
[Apache 2.0](https://choosealicense.com/licenses/apache-2.0/)
