# Phylogenetic analysis for DFE MS (Sane et al 2020)

This repository contains the R script and files to reproduce the panels of Figure 5 (pre-Illustrator adjustments). The code also simulates the character histories for 11 enzymes in the main manuscript. 

Brief description of files in the data folder

`20161115_pruned_tree_cutoff_0_03.nwk` - Phylogenetic tree for 1093 bacteria used in this study

`phylo_full_table.csv` - Comma separated file containing the species name and NCBI accession for the bacteria in the above tree. This file also contains the presence-absence matrix of the 11 repair enzymes in the main manuscript as well as the corresponding values for scaled biases

`20180119_all_repair_enzymes_pa_matrix_evalue_cutoff.csv` - Comma separated file with the presence-absence matrix of other repair enzymes used in the Supplementary table

`XXX_enz_transitions.csv` - Files with these names contain the state of each enzyme at every node in the tree as inferred from ancestral reconstruction. The raw files of the ancestral reconstruction itself are very large and hence only the processed files have been included. The ancestral reconstruction and generation of the transitions file can be done using the following code:
```{r}
library(phytools)

final_tree <- read.newick("data/20161115_pruned_tree_cutoff_0_03.nwk")
enz_pa_matrix <- read.csv("data/phylo_full_table.csv", colnames = 1)
enz_oi <- colnames(enz_pa_matrix)[6:16]

enz_simmap <- list()
for(i in 1:length(enz_oi) {
  enz_pa <- enz_pa_matrix[,i]
  rownames(enz_pa) <- enz_pa_matrix[,2]
  enz_simmap[[i]] <- make.simmap(tree = final_tree, x = enz_pa, model = "ARD", nsim = 500, pi="estimated", Q="mcmc")
}

simmap_summaries <- sapply(enz_oi, function(x) NULL)
enz_pps <- sapply(enz_oi, function(x) NULL)
enz_transitions <- sapply(enz_oi, function(x) NULL)

for(i in 1:length(enz_oi)){
  print(i)
  if(length(summary(enz_simmaps[[i]]))==3) next
  simmap_summaries[[i]] <- summary(enz_simmaps[[i]])
  enz_pps[[i]] <- rbind.data.frame(simmap_summaries[[i]]$tips, simmap_summaries[[i]]$ace)
  enz_transitions[[i]] <- as.data.frame(matrix(0, nrow=nrow(tree_oi_mod$edge), ncol = 7))
  
  for(edge_oi in 1:nrow(tree_oi_mod$edge)){
    p_node <- tree_oi_mod$edge[edge_oi, 1]
    d_node <- tree_oi_mod$edge[edge_oi, 2]
    p_states <- enz_pps[[i]][p_node,]
    d_states <- enz_pps[[i]][d_node,]
    enz_transitions[[i]][edge_oi, 1] <- paste0(p_node,"_",d_node)
    if(length(which(p_states>=0.7))>0) {
      enz_transitions[[i]][edge_oi, 2] <- as.numeric(names(p_states[which(p_states>=0.7)]))
    } else if(p_node==1094) {
      enz_transitions[[i]][edge_oi, 2] <- as.numeric(names(p_states[which.max(p_states)]))
    } else{
      same_as_before <- grep(p_node,tree_oi_mod$edge[,2])
      enz_transitions[[i]][edge_oi, 2] <- as.numeric(enz_transitions[[i]][same_as_before, 2])
    }
    if(length(which(d_states>=0.7))>0) {
      enz_transitions[[i]][edge_oi, 3] <- as.numeric(names(d_states[which(d_states>=0.7)]))
    } else {
      enz_transitions[[i]][edge_oi, 3] <- as.numeric(enz_transitions[[i]][edge_oi, 2])
    }
    enz_transitions[[i]][edge_oi, 5] <- as.numeric(p_states[which.max(p_states)])
    enz_transitions[[i]][edge_oi, 6] <- as.numeric(d_states[which.max(d_states)])
    if((enz_transitions[[i]][edge_oi, 2]!=enz_transitions[[i]][edge_oi, 3])) { 
      enz_transitions[[i]][edge_oi, 4] <- 1
      }
    if(enz_transitions[[i]][edge_oi, 6]>enz_transitions[[i]][edge_oi, 5]) {
        enz_transitions[[i]][edge_oi, 7] <- 1
      }
  }
  colnames(enz_transitions[[i]]) <- c("p_d_node","p_state","d_state","change","p_pp","d_pp","dpp>ppp")
  
  write.csv(x = enz_transitions[[i]], file = paste0("data/20171011_",enz_oi[i],"_transitions.csv"))
}
```
