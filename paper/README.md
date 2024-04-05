# Reproducing paper figures

To build the models for that were used to obtain the COO predictions for Figure 1B, run lines 2-16 (note that you will need to set things up to run jobs in parallel, as explained in the README of the homepage of this repo). To plot the results after building the models, run line 18. This will create a PDF in `figures` called `grid_analysis.pdf`.

To build the models for that were used to obtain the COO predictions for Figure 1C, run lines 23-29. These will also create bash scripts containing commands for plotting the results (in `analysis/ML/robustness_scripts`). Each script will be named with a unique ID. Run these from `analysis/ML` e.g assuming the unique ID is `robustness_c9e9655b-a9b5-4462-a064-db7f69e33ec7`, run 
```sh robustness_scripts/robustness_c9e9655b-a9b5-4462-a064-db7f69e33ec7.sh``` 
from within `analysis/ML`. 


