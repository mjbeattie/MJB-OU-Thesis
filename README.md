# MJB-OU-Thesis

This repository contains the files necessary to duplicate the analysis
in MJ Beattie's thesis for the University of Oklahoma, "Combining
Classification and Bayesian Methods to Better Model Drug Abuse."  The
files are contained in several directories as described below.

Thesis:  
The thesis itself in pdf and docx format.  Also included are the tables
in thesis in an Excel workbook.  The presentation from the thesis
defense is here as well.

Prep and PCA:  
These routines take the NSDUH 2002-2014 combined dataset and reduces it
for analysis.  prepdata_prePCA.R is run first, then NSDUH_pca.R is run.
The second program runs PCA analyses and produces the files for the
next stage in the study.

Cond Tree:
These routines conduct conditional inference random forest analyses on
the adult and youth datasets. cmatrix_cf3.R generates ROC curves for
the test datasets.

MCMC:
These routines run Bayesian MCMC analysis using the JAGS based runjags
code.  DBDA2E-utilities.R is a rountine required for graphing.
polyvalid.R and polyvalid_yth.R generate ROC curves and calculate
AUROC.
