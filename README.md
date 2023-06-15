# MVMR-CVD Online Supplementary Material

This online repository contains supplementary material for summary MR and MVMR analyses of coronary heart disease (CHD) and ischemic stroke.

## CVD_Review.csv

`CVD_Review` is a .csv file containing papers exploring causal determinants of CHD and stroke, corresponding to Figure 1 (right) and satisfying the following criteria:

1. Are accessible and written in English

2. Are not a review or meta-analysis

3. Are not animal or cell model studies

4. Use controls which do not have known pre-existing CVD conditions.

5. Satisfy Joanna Briggs Institute criteria for paper quality.

## Summary_Analysis.R

`Summary_Analysis` is a .R file with which summary MR and MVMR analyses can be replicated. The code is commented and includes the forward selection algorithm for selecting an appropriate MVMR model (see Methods).

## Incident_Analysis.R

`Incident_Analysis` is a .R file with which summary MR and MVMR analyses using summary GWAS data of incident ischemic stroke can be replicated. The code is commented and includes the forward selection algorithm for selecting an appropriate MVMR model (see Methods). The code requires files `summary_ISdat_v4_1_formatted.csv` and `summary_ISdat_v4_2_formatted` which contain a subset of GWAS summary data corresponding to the set of instruments, sourced from the incident stroke GWAS.


## summary_ISdat_v4_1_formatted.csv

`summary_ISdat_v4_1_formatted` is a .csv file, containing GWAS summary data for ischemic stroke corresponding to exposure instruments used in univariable summary MR analyses 

## summary_ISdat_v4_2_formatted.csv

`summary_ISdat_v4_2_formatted` is a .csv file, containing GWAS summary data for ischemic stroke corresponding to exposure instruments used in univariable summary MVMR analyses 

## Citation

The preprint version of this publication is available at:

## License

This project is licensed under GNU GPL v3.
