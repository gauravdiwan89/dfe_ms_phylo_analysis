# Phylogenetic analysis for DFE MS (Sane et al 2020)

This repository contains the R script and files to reproduce the panels of Figure 5 (pre-Illustrator adjustments). The code also simulates the character histories for 11 enzymes in the main manuscript. 

Brief description of files in the data folder

`20161115_pruned_tree_cutoff_0_03.nwk` - Phylogenetic tree for 1093 bacteria used in this study

`phylo_full_table.csv` - Comma separated file containing the species name and NCBI accession for the bacteria in the above tree. This file also contains the presence-absence matrix of the 11 repair enzymes in the main manuscript as well as the corresponding values for scaled biases

`20180119_all_repair_enzymes_pa_matrix_evalue_cutoff.csv` - Comma separated file with the presence-absence matrix of other repair enzymes used in the Supplementary table

`XXX_enz_transitions.csv` - Files with these names contain the state of each enzyme at every node in the tree as inferred from ancestral reconstruction. The raw files of the ancestral reconstruction itself are very large and hence only the processed files have been included. The ancestral reconstruction can be done using the following code:
```{r}
library(phytools)


```
