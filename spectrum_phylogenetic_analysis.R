###Load packages####

library(phytools)
library(seqinr)
library(tidyverse)
library(cowplot)

###Load tree and character histories of repair enzymes####
final_tree <- read.tree("data/20161115_pruned_tree_cutoff_0_03.nwk")
enz_oi <- c("uvrA", "nfi", "mutT", "mutY", "mutL", "mutS", "mutH", "uvrD", "dnaQ", "ung", "mutM")

enz_transitions <- sapply(enz_oi,function(x) NULL)
anc_state <- sapply(enz_oi,function(x) NULL)
lost_nodes <- sapply(enz_oi, function(x) NULL)
gain_nodes <- sapply(enz_oi, function(x) NULL)
for(i in 1:length(enz_oi)) {
  if(file.exists(paste0("data/20180206_", enz_oi[i], "_transitions.csv"))) {
    enz_transitions[[enz_oi[i]]] <- read.csv(paste0("data/20180206_",enz_oi[i],"_transitions.csv"),header = T,row.names = 1)
    changed_idx <- which(enz_transitions[[enz_oi[i]]]$change==1 & enz_transitions[[enz_oi[i]]]$d_state==0)
    d_nodes <- as.numeric(sapply(enz_transitions[[enz_oi[i]]]$p_d_node, function(x) strsplit(x, split = "_")[[1]][2]))
    lost_nodes[[enz_oi[i]]] <- d_nodes[changed_idx]
    changed_idx <- which(enz_transitions[[enz_oi[i]]]$change==1 & enz_transitions[[enz_oi[i]]]$d_state==1)
    gain_nodes[[enz_oi[i]]] <- d_nodes[changed_idx]
    anc_state[[enz_oi[i]]] <- enz_transitions[[enz_oi[i]]]$p_state[grep((Ntip(final_tree)+1),enz_transitions[[enz_oi[i]]]$p_d_node)[1]]
  }
  
}

###Input relative Tv/Ts (ts_tv) and GC/AT (at_gc) biases####
ts_tv <- c(1.055, 1.022, 2.231, 2.122, 0.072, 0.050, 0.066, 0.096, 0.313, 0.356, 1.258)
names(ts_tv) <- enz_oi
ts_tv_multi <- c(
  mutL_mutY = 0.182,
  mutL_dnaQ = 0.138,
  mutM_mutY = 2.225,
  mutL_mutS = 0.066,
  mutL_mutS_mutH = 0.051
)

scaled_ts_tv <- ifelse(ts_tv<1, -((-1+ts_tv)/min(-1+ts_tv)), ((-1+ts_tv)/max(-1+ts_tv)))
scaled_ts_tv_multi <- sapply(ts_tv_multi, function(x) ifelse(x<1, -((-1+x)/min(-1+ts_tv)), ((-1+x)/max(-1+ts_tv))))

at_gc <- c(1.003, 0.957, 0.004, 1.717, 0.416, 0.375, 0.528, 0.607, 0.787, 1.494, 1.243)
names(at_gc) <- enz_oi
at_gc_multi <- c(
  mutL_mutY = 0.427,
  mutL_dnaQ = 1.134,
  mutM_mutY = 1.788,
  mutL_mutS = 0.295,
  mutL_mutS_mutH = 0.551
)

scaled_at_gc <- ifelse(at_gc<1, -((-1+at_gc)/min(-1+at_gc)), ((-1+at_gc)/max(-1+at_gc)))
scaled_at_gc_multi <- sapply(at_gc_multi, function(x) ifelse(x<1, -((-1+x)/min(-1+at_gc)), ((-1+x)/max(-1+at_gc))))

###Make phylogenetic heatmap and map enzyme transitions on the tree - Figure 5A####
orig_pa_matrix <- read.csv("data/phylo_full_table.csv", header = T, row.names = NULL)
cols_oi <- as.numeric(na.omit(match(enz_oi, colnames(orig_pa_matrix))))
new_pa_matrix <- orig_pa_matrix[ ,cols_oi]
rownames(new_pa_matrix) <- orig_pa_matrix$phylotree_tip_label


pdf(file = "phylo_heatmap_11_mut_enz.pdf",width = 20,height = 30)
phylo.heatmap(tree = final_tree, X = new_pa_matrix[,names(sort(ts_tv, decreasing = T))], fsize = c(0.1,1,0.1), colors = c(0,gray(0.6)), legend=F, split=c(0.6,0.4), lwd=0.5, mar=c(5,2,2,2))
dev.off()

enz_cols <- c("#c632b6","#00e189","#833f23","#8777ff","#54560e","#6935ad","#ffa545","#dfafdd","#ff7141","#95acff","#d4003a")

pdf(file = "enz_trans_for_edit.pdf",width = 12,height = 27)
par(mar=c(5,2,2,2))
# plotBranchbyTrait(final_tree,x = gc_anc_reorder,mode = "nodes",show.tip.label=F,palette = "gray",edge.width=1)
plot(final_tree, show.tip.label = F, edge.width = 1)
plot()
ver_adj <- rev(seq(-200,200,400/11))
hor_adj <- 0.4
anc_state <- sapply(enz_transitions, function(x) x$p_state[1])
for(i in 1:length(enz_oi)) {
  if(enz_oi[i]=="recBp") next
  if(anc_state[i]==0) {enz_pch <- 1} else enz_pch <- 19
  nodelabels(text = NULL,node = (Ntip(final_tree)+1),adj = c(hor_adj,ver_adj[i]),pch=enz_pch,lwd=2,col=enz_cols[i],cex = 3)
}
# tiplabels(tip = 1:length(final_tree$tip.label),pch = 15,col=col_by_family)
legend("bottomright",legend = enz_oi,col = enz_cols[1:11],pch=19,bty="n", cex = 2)

hor_adj <- seq(0.5,0.7,by=0.2/11)
for(j in 1:length(enz_oi)) {
  if(enz_oi[j]=="recBp") next
  where_changes <- which(enz_transitions[[enz_oi[j]]]$change==1)
  for(i in 1:length(where_changes)) {
    who_is_higher <- which.max(c(enz_transitions[[enz_oi[j]]]$p_state[where_changes[i]],enz_transitions[[enz_oi[j]]]$d_state[where_changes[i]]))
    if(who_is_higher==2) {enz_pch <- 19} else if(who_is_higher==1) {enz_pch <- 1}
    edgelabels(text = NULL,edge = where_changes[i],adj = hor_adj[j],col = enz_cols[j],pch = enz_pch,lwd=2,cex=2)
  }
}
dev.off()


###Supplementary figure - other enzyme transitions on the tree####
orig_pa_matrix <- read.csv("data/20180119_all_repair_enzymes_pa_matrix_evalue_cutoff.csv", header = T, row.names = NULL)

rem_enz <- colnames(orig_pa_matrix)[9:51]
rem_enz <- rem_enz[!rem_enz %in% enz_oi]

rem_enz_transitions <- sapply(rem_enz,function(x) NULL)
rem_anc_state <- sapply(rem_enz,function(x) NULL)
total_gains_m <- sapply(rem_enz,function(x) NULL)
total_losses_m <- sapply(rem_enz,function(x) NULL)

for(i in 1:length(rem_enz)) {
  if(rem_enz[i]=="recBp") next
  rem_enz_transitions[[rem_enz[i]]] <- read.csv(paste0("data/other_enzyme_files/20180206_",rem_enz[i],"_transitions.csv"),header = T,row.names = 1)
  rem_anc_state[[rem_enz[i]]] <- rem_enz_transitions[[rem_enz[i]]]$p_state[grep((Ntip(final_tree)+1),rem_enz_transitions[[rem_enz[i]]]$p_d_node)[1]]
  total_gains_m[[rem_enz[i]]] <- length(which(rem_enz_transitions[[rem_enz[i]]]$p_state==0&rem_enz_transitions[[rem_enz[i]]]$d_state==1))
  total_losses_m[[rem_enz[i]]] <- length(which(rem_enz_transitions[[rem_enz[i]]]$p_state==1&rem_enz_transitions[[rem_enz[i]]]$d_state==0))
}

write.csv(cbind(rem_enz, total_gains_m, total_losses_m), file = "rem_enz_gains_losses.csv", row.names = F)

enz_cols <- c("#ff4b85","#78dc69","#8967f2","#b5b900","#0181ec","#e6b500","#822a95","#007101","#e91e96","#00aa76","#ad008e","#4c7c00","#ebaaff","#807a00","#8ba1ff","#d07700","#02c2ff","#bf1f05","#006647","#ed2548","#b4cf94","#ff4764","#dec657","#654678","#ffae55","#8e2d62","#f3bb7c","#a60c26","#6b6330","#ff8bb3","#8c3913","#8c4d69","#ff844c","#a57453")

pdf(file = "rem_mut_enz_transitions.pdf",width = 12,height = 30)
par(mar=c(5,2,2,2))
# plotBranchbyTrait(test_tree,x = gc_anc_reorder,mode = "nodes",show.tip.label=F,palette = "heat.colors",edge.width=1)
plot(final_tree, show.tip.label = F, edge.width = 1)
ver_adj <- rev(seq(-200,600,800/50))
hor_adj <- 0.4
for(i in 1:length(rem_enz)) {
  if(rem_enz[i]=="recBp") next
  if(rem_anc_state[i]==0) {enz_pch <- 1} else enz_pch <- 19
  nodelabels(text = NULL,node = (Ntip(final_tree)+1),adj = c(hor_adj,ver_adj[i]),pch=enz_pch,lwd=2,col=enz_cols[i],cex = 3)
}
# tiplabels(tip = 1:length(final_tree$tip.label),pch = 15,col=col_by_family)
legend("bottomright",legend = rem_enz,col = enz_cols,pch=19,bty="n")

hor_adj <- seq(0.5,0.7,by=0.2/43)
for(j in 1:length(rem_enz)) {
  if(rem_enz[j]=="recBp") next
  where_changes <- which(rem_enz_transitions[[rem_enz[j]]]$change==1)
  for(i in 1:length(where_changes)) {
    who_is_higher <- which.max(c(rem_enz_transitions[[rem_enz[j]]]$p_state[where_changes[i]],rem_enz_transitions[[rem_enz[j]]]$d_state[where_changes[i]]))
    if(who_is_higher==2) {enz_pch <- 19} else if(who_is_higher==1) {enz_pch <- 1}
    edgelabels(text = NULL,edge = where_changes[i],adj = hor_adj[j],col = enz_cols[j],pch = enz_pch,lwd=2,cex=2)
  }
}
dev.off()



rem_enz_pa_matrix <- orig_pa_matrix[,rem_enz[rem_enz!="recBp"]]
rownames(rem_enz_pa_matrix) <- orig_pa_matrix$phylotree_tip_label

pdf(file = "phylo_heatmap_more_mut_enz_evalue_cutoff.pdf",width = 15,height = 30)
phylo.heatmap(tree = final_tree,X = rem_enz_pa_matrix,fsize = c(0.1,1,0.1),colors = c(0,gray(0.6)),legend=F,split=c(0.6,0.4),lwd=0.5,mar=c(5,2,2,2))
dev.off()

###Calculating aggregate biases at every node in every lineage of the tree####
##define lineages and normalised distances from the tree
tip_heights <- sapply(1:Ntip(final_tree), function(x) nodeheight(final_tree, x))
node_vec <- (Ntip(final_tree)+1):(Ntip(final_tree)+Nnode(final_tree))
desc_list <- sapply(node_vec, function(y) getDescendants(final_tree, node = y))

norm_node_heights2 <- vector("list", length(node_vec))
names(norm_node_heights2) <- as.vector(node_vec)
all_parent_nodes <- list()
for(i in 1:Ntip(final_tree)) {
  parent_nodes <- node_vec[which(sapply(desc_list, function(x) i %in% x))]
  all_parent_nodes[[i]] <- parent_nodes
  norm_hts <- sapply(parent_nodes, function(x) nodeheight(final_tree, x))/tip_heights[i]
  names(norm_hts) <- as.character(parent_nodes)
  for(j in 1:length(parent_nodes)) {
    norm_node_heights2[[as.character(parent_nodes[j])]] <- c(norm_node_heights2[[as.character(parent_nodes[j])]], norm_hts[as.character(parent_nodes[j])])
  }
}

avg_norm_node_heights <- sapply(norm_node_heights2, min)

##Tv/Ts
all_sums_list <- list()
tip_norm_dist <- list()
for(i in 1:length(all_parent_nodes)) {
  ##all
  all_state_mat2 <- sapply(enz_oi, function(x) 1-enz_transitions[[x]]$d_state[match(all_parent_nodes[[i]], d_nodes)])
  all_state_mat2[1, ] <- 1-unlist(anc_state)
  all_state_mat2 <- rbind(all_state_mat2, 1-new_pa_matrix[final_tree$tip.label[[i]], enz_oi])
  rownames(all_state_mat2) <- as.character(c(all_parent_nodes[[i]], i))
  all_sums2 <- vector(mode = "numeric", length = nrow(all_state_mat2))
  names(all_sums2) <- rownames(all_state_mat2)
  for(j in 1:nrow(all_state_mat2)) {
    ##double triple checks
    state_check <- c(
      mutL_mutS_mutH = 3-sum(all_state_mat2[j, c("mutL", "mutS", "mutH")]),
      # mutL_mutS_mutY = 3-sum(all_state_mat2[j, c("mutL", "mutS", "mutY")]),
      mutL_mutY = 2-sum(all_state_mat2[j, c("mutL", "mutY")]),
      mutL_mutS = 2-sum(all_state_mat2[j, c("mutL", "mutS")]),
      mutM_mutY = 2-sum(all_state_mat2[j, c("mutM", "mutY")]),
      mutL_dnaQ = 2-sum(all_state_mat2[j, c("mutL", "dnaQ")])
    )
    if(any(state_check == 0)) {
      multi_kos <- names(state_check)[state_check==0]
      # if("mutL_mutS_mutY" %in% multi_kos) {
      #   multi_kos <- multi_kos[!multi_kos %in% c("mutL_mutS", "mutL_mutY")]
      # } 
      if("mutL_mutS_mutH" %in% multi_kos) {
        multi_kos <- multi_kos[!multi_kos %in% c("mutL_mutS")]
      }
      multi_sums <- sum((scaled_ts_tv_multi[multi_kos]))
      enz_covered <- unique(unlist(sapply(multi_kos, function(x) strsplit(x, split = "_")[[1]])))
      rem_enz <- enz_oi[!enz_oi %in% enz_covered]
      rem_sums <- sum(all_state_mat2[j,rem_enz]*(scaled_ts_tv[rem_enz]))
      all_sums2[j] <- (rem_sums + multi_sums)
    } else {
      all_sums2[j] <- sum(all_state_mat2[j,]*(scaled_ts_tv))
    }
  }
  all_sums_list[[i]] <- all_sums2
  
  tip_norm_dist[[i]] <- c(avg_norm_node_heights[as.character(all_parent_nodes[[i]])],1)
  names(tip_norm_dist[[i]])[length(tip_norm_dist[[i]])] <- as.character(i)
}

##GC/AT
all_sums_list_at_gc <- list()
for(i in 1:length(all_parent_nodes)) {
  ##all
  all_state_mat3 <- sapply(enz_oi, function(x) 1-enz_transitions[[x]]$d_state[match(all_parent_nodes[[i]], d_nodes)])
  all_state_mat3[1, ] <- 1-unlist(anc_state)
  all_state_mat3 <- rbind(all_state_mat3, 1-new_pa_matrix[final_tree$tip.label[[i]], enz_oi])
  rownames(all_state_mat3) <- as.character(c(all_parent_nodes[[i]], i))
  all_sums3 <- vector(mode = "numeric", length = nrow(all_state_mat3))
  names(all_sums3) <- rownames(all_state_mat3)
  for(j in 1:nrow(all_state_mat3)) {
    ##double triple checks
    state_check <- c(
      mutL_mutS_mutH = 3-sum(all_state_mat3[j, c("mutL", "mutS", "mutH")]),
      # mutL_mutS_mutY = 3-sum(all_state_mat3[j, c("mutL", "mutS", "mutY")]),
      mutL_mutY = 2-sum(all_state_mat3[j, c("mutL", "mutY")]),
      mutL_mutS = 2-sum(all_state_mat3[j, c("mutL", "mutS")]),
      mutM_mutY = 2-sum(all_state_mat3[j, c("mutM", "mutY")]),
      mutL_dnaQ = 2-sum(all_state_mat3[j, c("mutL", "dnaQ")])
    )
    if(any(state_check == 0)) {
      multi_kos <- names(state_check)[state_check==0]
      # if("mutL_mutS_mutY" %in% multi_kos) {
      #   multi_kos <- multi_kos[!multi_kos %in% c("mutL_mutS", "mutL_mutY")]
      # } 
      if("mutL_mutS_mutH" %in% multi_kos) {
        multi_kos <- multi_kos[!multi_kos %in% c("mutL_mutS")]
      }
      multi_sums <- sum((scaled_at_gc_multi[multi_kos]))
      enz_covered <- unique(unlist(sapply(multi_kos, function(x) strsplit(x, split = "_")[[1]])))
      rem_enz <- enz_oi[!enz_oi %in% enz_covered]
      rem_sums <- sum(all_state_mat3[j,rem_enz]*(scaled_at_gc[rem_enz]))
      all_sums3[j] <- (rem_sums + multi_sums)
    } else {
      all_sums3[j] <- sum(all_state_mat3[j,]*(scaled_at_gc))
    }
  }
  all_sums_list_at_gc[[i]] <- all_sums3
  
}

##GC/AT - scaled bias values corrected for by GC content of each node

# gc_anc <- read.csv("F:/Idea_new_backup/Ortho_detection/StableTraits_output.csv", header = T, row.names = 1)
# 
# rows_to_tips <- match(gc_anc$Parameter[c(-1, -2)], final_tree$tip.label)
# node_nos  <- Ntip(final_tree) + as.numeric(gsub("n", "", gc_anc$Parameter)) + 1
# 
# rows_to_tips[is.na(rows_to_tips)] <- node_nos[!is.na(node_nos)]
# 
# gc_anc_values <- gc_anc$Median[c(-1, -2)]
# names(gc_anc_values) <- as.character(rows_to_tips)
# 
# all_sums_list_at_gc <- list()
# for(i in 1:length(all_parent_nodes)) {
#   ##all
#   all_state_mat3 <- sapply(enz_oi, function(x) 1-enz_transitions[[x]]$d_state[match(all_parent_nodes[[i]], d_nodes)])
#   all_state_mat3[1, ] <- 1-unlist(anc_state)
#   all_state_mat3 <- rbind(all_state_mat3, 1-new_pa_matrix[final_tree$tip.label[[i]], enz_oi])
#   rownames(all_state_mat3) <- as.character(c(all_parent_nodes[[i]], i))
#   all_sums3 <- vector(mode = "numeric", length = nrow(all_state_mat3))
#   names(all_sums3) <- rownames(all_state_mat3)
#   for(j in 1:nrow(all_state_mat3)) {
#     ##double triple checks
#     state_check <- c(
#       mutL_mutS_mutH = 3-sum(all_state_mat3[j, c("mutL", "mutS", "mutH")]),
#       # mutL_mutS_mutY = 3-sum(all_state_mat3[j, c("mutL", "mutS", "mutY")]),
#       mutL_mutY = 2-sum(all_state_mat3[j, c("mutL", "mutY")]),
#       mutL_mutS = 2-sum(all_state_mat3[j, c("mutL", "mutS")]),
#       mutM_mutY = 2-sum(all_state_mat3[j, c("mutM", "mutY")]),
#       mutL_dnaQ = 2-sum(all_state_mat3[j, c("mutL", "dnaQ")])
#     )
#     if(any(state_check == 0)) {
#       multi_kos <- names(state_check)[state_check==0]
#       # if("mutL_mutS_mutY" %in% multi_kos) {
#       #   multi_kos <- multi_kos[!multi_kos %in% c("mutL_mutS", "mutL_mutY")]
#       # } 
#       if("mutL_mutS_mutH" %in% multi_kos) {
#         multi_kos <- multi_kos[!multi_kos %in% c("mutL_mutS")]
#       }
#       
#       gc_oi <- gc_anc_values[rownames(all_state_mat3)[j]]
#       gc_corrected_at_gc <- ifelse(scaled_at_gc>0, scaled_at_gc*gc_oi, scaled_at_gc*(1-gc_oi))
#       gc_corrected_at_gc_multi <- ifelse(scaled_at_gc_multi>0, scaled_at_gc_multi*gc_oi, scaled_at_gc_multi*(1-gc_oi))
#       
#       multi_sums <- sum((gc_corrected_at_gc_multi[multi_kos]))
#       enz_covered <- unique(unlist(sapply(multi_kos, function(x) strsplit(x, split = "_")[[1]])))
#       rem_enz <- enz_oi[!enz_oi %in% enz_covered]
#       rem_sums <- sum(all_state_mat3[j,rem_enz]*(gc_corrected_at_gc[rem_enz]))
#       all_sums3[j] <- (rem_sums + multi_sums)
#     } else {
#       gc_oi <- gc_anc_values[rownames(all_state_mat3)[j]]
#       gc_corrected_at_gc <- ifelse(scaled_at_gc>0, scaled_at_gc*gc_oi, scaled_at_gc*(1-gc_oi))
#       gc_corrected_at_gc_multi <- ifelse(scaled_at_gc_multi>0, scaled_at_gc_multi*gc_oi, scaled_at_gc_multi*(1-gc_oi))
#       
#       all_sums3[j] <- sum(all_state_mat3[j,]*(gc_corrected_at_gc))
#     }
#   }
#   all_sums_list_at_gc[[i]] <- all_sums3
#   
# }

##determine types of events
bias_diff <- sapply(all_sums_list, diff)
ups_downs <- sapply(bias_diff, function(x) x[x != 0])

no_of_ups <- sapply(ups_downs, function(x) sum(x > 0))
no_of_downs <- sapply(ups_downs, function(x) sum(x < 0))

no_changes <- length(which(sapply(bias_diff, function(x) all(x==0))))
one_change <- length(which(sapply(ups_downs, length) == 1))
reinforced <- length(which(sapply(ups_downs, function(x) all(x < 0) | all(x > 0))))
reinforced <- reinforced - one_change - no_changes ##they are considered in reinforced because they are NULL

cons_up_event <- NULL
cons_down_event <- NULL
reinforce_to_reverse <- NULL
reverse_to_reinforce <- NULL
for(i in 1:length(ups_downs)) {
  cons_up_event[i] <- 0
  cons_down_event[i] <- 0
  reinforce_to_reverse[i] <- 0
  reverse_to_reinforce[i] <- 0
  if(length(ups_downs[[i]])>1) {
    for(j in 1:(length(ups_downs[[i]])-1)) {
      if(ups_downs[[i]][j]<0 & ups_downs[[i]][j+1]>0) {
        cons_up_event[i] <- cons_up_event[i]+1
      } else if(ups_downs[[i]][j]>0 & ups_downs[[i]][j+1]<0) {
        cons_down_event[i] <- cons_down_event[i]+1
      }
    }
  }
  if(length(ups_downs[[i]])>2) {
    for(j in 1:(length(ups_downs[[i]])-2)) {
      if(((ups_downs[[i]][j]<0 & ups_downs[[i]][j+1]<0) & ups_downs[[i]][j+2]>0) | ((ups_downs[[i]][j]>0 & ups_downs[[i]][j+1]>0) & ups_downs[[i]][j+2]<0)) {
        reinforce_to_reverse[i] <- reinforce_to_reverse[i]+1
      } else if(((ups_downs[[i]][j]<0 & ups_downs[[i]][j+1]>0) & ups_downs[[i]][j+2]>0) | ((ups_downs[[i]][j]>0 & ups_downs[[i]][j+1]<0) & ups_downs[[i]][j+2]<0)) {
        reverse_to_reinforce[i] <- reverse_to_reinforce[i]+1
      }
    }
  }
}
total_fluc <- cons_up_event+cons_down_event

##same for GC -> AT
bias_diff_at_gc <- sapply(all_sums_list_at_gc, diff)
ups_downs_at_gc <- sapply(bias_diff_at_gc, function(x) x[x != 0])

no_of_ups_at_gc <- sapply(ups_downs_at_gc, function(x) sum(x > 0))
no_of_downs_at_gc <- sapply(ups_downs_at_gc, function(x) sum(x < 0))

no_changes_at_gc <- length(which(sapply(bias_diff_at_gc, function(x) all(x==0))))
one_change_at_gc <- length(which(sapply(ups_downs_at_gc, length) == 1))
reinforced_at_gc <- length(which(sapply(ups_downs_at_gc, function(x) all(x < 0) | all(x > 0))))
reinforced_at_gc <- reinforced_at_gc - one_change_at_gc - no_changes_at_gc

cons_up_event_at_gc <- NULL
cons_down_event_at_gc <- NULL
reinforce_to_reverse_at_gc <- NULL
reverse_to_reinforce_at_gc <- NULL
for(i in 1:length(ups_downs_at_gc)) {
  cons_up_event_at_gc[i] <- 0
  cons_down_event_at_gc[i] <- 0
  reinforce_to_reverse_at_gc[i] <- 0
  reverse_to_reinforce_at_gc[i] <- 0
  if(length(ups_downs_at_gc[[i]])>1) {
    for(j in 1:(length(ups_downs_at_gc[[i]])-1)) {
      if(ups_downs_at_gc[[i]][j]<0 & ups_downs_at_gc[[i]][j+1]>0) {
        cons_up_event_at_gc[i] <- cons_up_event_at_gc[i]+1
      } else if(ups_downs_at_gc[[i]][j]>0 & ups_downs_at_gc[[i]][j+1]<0) {
        cons_down_event_at_gc[i] <- cons_down_event_at_gc[i]+1
      }
    }
  }
  if(length(ups_downs_at_gc[[i]])>2) {
    for(j in 1:(length(ups_downs_at_gc[[i]])-2)) {
      if(((ups_downs_at_gc[[i]][j]<0 & ups_downs_at_gc[[i]][j+1]<0) & ups_downs_at_gc[[i]][j+2]>0) | ((ups_downs_at_gc[[i]][j]>0 & ups_downs_at_gc[[i]][j+1]>0) & ups_downs_at_gc[[i]][j+2]<0)) {
        reinforce_to_reverse_at_gc[i] <- reinforce_to_reverse_at_gc[i]+1
      } else if(((ups_downs_at_gc[[i]][j]<0 & ups_downs_at_gc[[i]][j+1]>0) & ups_downs_at_gc[[i]][j+2]>0) | ((ups_downs_at_gc[[i]][j]>0 & ups_downs_at_gc[[i]][j+1]<0) & ups_downs_at_gc[[i]][j+2]<0)) {
        reverse_to_reinforce_at_gc[i] <- reverse_to_reinforce_at_gc[i]+1
      }
    }
  }
}
total_fluc_at_gc <- cons_up_event_at_gc+cons_down_event_at_gc

###Figure 5C####

ts_table <- tibble(
  category = c("No change", "Single change", "Reinforced", "Reversals", "Reversal after\nreinforcement", "Reinforcement\nafter reversal"),
  `events` = c(no_changes, one_change, reinforced, 0, sum(reverse_to_reinforce>0), sum(reinforce_to_reverse>0)),
  `1` = c(0, 0, 0, table(total_fluc)["1"], 0, 0),
  `2` = c(0, 0, 0, table(total_fluc)["2"], 0, 0),
  `3` = c(0, 0, 0, table(total_fluc)["3"], 0, 0),
  `4` = c(0, 0, 0, table(total_fluc)["4"], 0, 0),
  `5` = c(0, 0, 0, table(total_fluc)["5"], 0, 0),
  `6` = c(0, 0, 0, table(total_fluc)["6"], 0, 0),
  `7` = c(0, 0, 0, table(total_fluc)["7"], 0, 0),
  bias = "Tv/Ts bias"
)

gc_table <- tibble(
  category = c("No change", "Single change", "Reinforced", "Reversals", "Reversal after\nreinforcement", "Reinforcement\nafter reversal"),
  `events` = c(no_changes_at_gc, one_change_at_gc, reinforced_at_gc, 0, sum(reverse_to_reinforce_at_gc>0), sum(reinforce_to_reverse_at_gc>0)),
  `1` = c(0, 0, 0, table(total_fluc_at_gc)["1"], 0, 0),
  `2` = c(0, 0, 0, table(total_fluc_at_gc)["2"], 0, 0),
  `3` = c(0, 0, 0, table(total_fluc_at_gc)["3"], 0, 0),
  `4` = c(0, 0, 0, table(total_fluc_at_gc)["4"], 0, 0),
  `5` = c(0, 0, 0, table(total_fluc_at_gc)["5"], 0, 0),
  `6` = c(0, 0, 0, table(total_fluc_at_gc)["6"], 0, 0),
  `7` = c(0, 0, 0, table(total_fluc_at_gc)["7"], 0, 0),
  `8` = c(0, 0, 0, table(total_fluc_at_gc)["8"], 0, 0),
  `9` = c(0, 0, 0, table(total_fluc_at_gc)["9"], 0, 0),
  `10` = c(0, 0, 0, table(total_fluc_at_gc)["10"], 0, 0),
  `11` = c(0, 0, 0, table(total_fluc_at_gc)["11"], 0, 0),
  `12` = c(0, 0, 0, table(total_fluc_at_gc)["12"], 0, 0),
  `13` = c(0, 0, 0, table(total_fluc_at_gc)["13"], 0, 0),
  `14` = c(0, 0, 0, table(total_fluc_at_gc)["14"], 0, 0),
  `15` = c(0, 0, 0, table(total_fluc_at_gc)["15"], 0, 0),
  `16` = c(0, 0, 0, table(total_fluc_at_gc)["16"], 0, 0),
  # `7` = c(0, 0, 0, 0, 0, 0),
  bias = "GC/AT bias"
)

btable <- rbind(ts_table %>% reshape2::melt(), gc_table %>% reshape2::melt())

btable$category <- factor(btable$category, levels = c("No change", "Single change", "Reinforced","Reinforcement\nafter reversal", "Reversal after\nreinforcement", "Reversals"))
btable$bias <- factor(btable$bias, levels = c("Tv/Ts bias", "GC/AT bias"))

bcols <- c("#8DA0CB", rev(colorRampPalette(c("rosybrown1", "red3"))(16)))

bplot_h <- btable %>% 
  filter(! category %in% c("Reinforcement\nafter reversal", "Reversal after\nreinforcement")) %>% 
  mutate(variable = factor(variable, levels = c("events", "16", "15", "14", "13", "12", "11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1"))) %>% 
  ggplot(aes(category, value, fill = variable)) +
  geom_col(width = 0.7) +
  # scale_fill_brewer(palette = "Set2", direction = -1) +
  scale_fill_manual(values = bcols) +
  facet_wrap(~bias, nrow = 2) +
  theme_cowplot(font_size = 22) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1093)) +
  labs(
    x = "",
    y = "Number of lineages",
    fill = ""
  ) +
  # theme(axis.text.y = element_text(angle = 45, hjust = 1)) + 
  coord_flip()
bplot_h

write.csv(rbind(ts_table, gc_table), file = "20200913_matrix_for_trajectory_summary_barplot.csv", row.names = F)

###Figure 5B####
cols <- c( "gray60", "gray40",  "#8DA0CB", "cornflowerblue", "darkslateblue", colorRampPalette(c("rosybrown1", "red3"))(7))

#tv/ts
plot_cols <- vector("character", length = Ntip(final_tree))
plot_cols[which(sapply(ups_downs, function(x) (all(x < 0) | all(x > 0)) & length(x) > 1))] <- cols[3]
plot_cols[which(sapply(ups_downs, length) == 1)] <- cols[2]
plot_cols[which(sapply(bias_diff, function(x) all(x == 0)))] <- cols[1]
for(i in 1:7) {
  plot_cols[total_fluc==i] <- cols[i+5]
}
# plot_cols[reverse_to_reinforce>0] <- cols[4]
# plot_cols[reinforce_to_reverse>0] <- cols[5]


par(mar=c(5,5,2,2))
plot(0, 0, xlim = c(0,1), ylim = c(-3.2, 3.2), type = "n", xlab = "Normalised distance from root", ylab = expression(paste(Sigma, " Scaled Tv/Ts bias")), cex.axis = 1.5, cex.lab = 1.5)
for(i in 1:Ntip(final_tree)) {
  lines(tip_norm_dist[[i]], all_sums_list[[i]], col = plot_cols[i], lwd = 3)
}
for(i in which(sapply(bias_diff, function(x) all(x == 0)))) {
  lines(tip_norm_dist[[i]], all_sums_list[[i]], col = plot_cols[i], lwd = 3)
}


##gc/at
plot_cols_at_gc <- vector("character", length = Ntip(final_tree))
plot_cols_at_gc[which(sapply(ups_downs_at_gc, function(x) (all(x < 0) | all(x > 0)) & length(x) > 1))] <- cols[3]
plot_cols_at_gc[which(sapply(ups_downs_at_gc, length) == 1)] <- cols[2]
plot_cols_at_gc[which(sapply(bias_diff_at_gc, function(x) all(x == 0)))] <- cols[1]
for(i in 1:7) {
  plot_cols_at_gc[which(total_fluc_at_gc==i)] <- cols[i+5]
}
# plot_cols_at_gc[reverse_to_reinforce_at_gc>0] <- cols[4]
# plot_cols_at_gc[reinforce_to_reverse_at_gc>0] <- cols[5]

par(mar=c(5,5,2,2))
plot(0, 0, xlim = c(0,1), ylim = c(-3.2, 3.2), type = "n", xlab = "Normalised distance from root", ylab = expression(paste(Sigma, " Scaled GC/AT bias")), cex.axis = 1.5, cex.lab = 1.5)
for(i in 1:Ntip(final_tree)) {
  lines(tip_norm_dist[[i]], all_sums_list_at_gc[[i]], col = plot_cols_at_gc[i], lwd = 3)
}
for(i in which(sapply(bias_diff_at_gc, function(x) all(x == 0)))) {
  lines(tip_norm_dist[[i]], all_sums_list_at_gc[[i]], col = plot_cols_at_gc[i], lwd = 3)
}
legend("topleft", legend = c("no change", "single change", "reinforced\nbias", "", 1:7), col = c(cols[1:3], "white", cols[6:12]), lwd = 4, bty = "n", horiz = F, x.intersp = 0.5, y.intersp = 0.1, cex = 1.5, ncol = 3)
text(0, 3.2, labels = "Lineages with:", cex = 1.5, adj = 0)
text(0.6, 3.2, labels = "No. of reversals", cex = 1.5, adj = 0.5)


###Make phylogenetic heatmap for aggreagate biases of extant species - Figure 5A####
tips_bias_ts_tv <- sapply(all_sums_list, function(x) x[length(x)])
tips_bias_at_gc <- sapply(all_sums_list_at_gc, function(x) x[length(x)])

##Tv/Ts
bias_matrix <- data.frame(tips_bias_ts_tv, rep(0, Ntip(final_tree)))
rownames(bias_matrix) <- final_tree$tip.label
colnames(bias_matrix) <- c("Tv/Ts bias", "GC/AT bias")

pdf(file = "20200913_phylo_heatmap_tv_ts.pdf",width = 20,height = 30)
phylo.heatmap(tree = final_tree,X = bias_matrix,fsize = c(0.1,1.5,1),colors = rev(RColorBrewer::brewer.pal(n = 9, name = "Blues")), legend=T, split=c(0.6,0.4), lwd=0.5, mar=c(5,2,2,2))
dev.off()

##GC/AT
bias_matrix <- data.frame(rep(0, Ntip(final_tree)), tips_bias_at_gc)
rownames(bias_matrix) <- final_tree$tip.label
colnames(bias_matrix) <- c("Tv/Ts bias", "GC/AT bias")

pdf(file = "20200913_phylo_heatmap_gc_at.pdf",width = 20,height = 30)
phylo.heatmap(tree = final_tree,X = bias_matrix,fsize = c(0.1,1.5,1),colors = rev(RColorBrewer::brewer.pal(n = 9, name = "Greens")), legend=T, split=c(0.6,0.4), lwd=0.5, mar=c(5,2,2,2))
dev.off()


###Simulating character histories using actual transition rate matrices and calculating expected number of different events####

trans_mat <- sapply(enz_oi, function(x) NULL)
for(i in 1:length(enz_oi)) {
  # if(enz_oi[i] == "mutT") next
  print(enz_oi[i])
  load(paste0("20180206_", enz_oi[i],"_simmap"))
  trans_mat[[enz_oi[i]]] <- lapply(mut_enz_simmaps, function(x) x$Q)
}

all_parent_nodes_tips <- sapply(1:length(all_parent_nodes), function(x) c(all_parent_nodes[[x]], x))

##Tv/Ts
random_enz_state <- sapply(enz_oi, function(x) NULL)

no_changes_r <- matrix(NA, nrow = 500, ncol = 20)
one_change_r <- matrix(NA, nrow = 500, ncol = 20)
reinforced_r <- matrix(NA, nrow = 500, ncol = 20)

total_fluc_r <- list()

for(j in 1:500) {
  print(j)
  # if(enz_oi[i] == "mutT") next
  start.time <- Sys.time()
  
  bias_diff_r <- list()
  ups_downs_r <- list()
  
  no_of_ups_r <- list()
  no_of_downs_r <- list()
  
  cons_up_event_r <- list()
  cons_down_event_r <- list()
  
  total_fluc_r[[j]] <- list()
  
  for(i in 1:length(enz_oi)) {
    print(enz_oi[i])
    test <- sim.history(tree = final_tree, Q = trans_mat[[enz_oi[i]]][[j]], anc = as.character(anc_state[enz_oi[i]]), nsim = 20, message = F)
    random_states <- matrix(NA, nrow = (Ntip(final_tree)+final_tree$Nnode), ncol = length(test), dimnames = list(c("1094", test[[1]]$edge[,2])))
    for(jj in 1:length(test)) {
      random_states[,jj] <- c(as.numeric(test[[jj]]$node.states[1,1]), as.numeric(test[[jj]]$node.states[,2]))
    }
    random_enz_state[[enz_oi[i]]] <- random_states
  }
  
  print("Doing calculations...")
  #calc
  for(jj in 1:length(test)) {
    random_sums_list <- list()
    for(ii in 1:length(all_parent_nodes)) {
      random_state_mat <- sapply(enz_oi, function(x) 1-random_enz_state[[x]][,jj][match(as.character(all_parent_nodes[[ii]]), names(random_enz_state[[x]][,jj]))])
      random_state_mat <- rbind(random_state_mat, 1-sapply(enz_oi, function(x) random_enz_state[[x]][,jj][names(random_enz_state[[x]][,jj])==as.character(ii)]))
      rownames(random_state_mat) <- as.character(c(all_parent_nodes[[ii]], ii))
      random_sums <- vector(mode = "numeric", length = nrow(random_state_mat))
      names(random_sums) <- rownames(random_state_mat)
      for(k in 1:nrow(random_state_mat)) {
        state_check <- c(
          mutL_mutS_mutH = 3-sum(random_state_mat[k, c("mutL", "mutS", "mutH")]),
          # mutL_mutS_mutY = 3-sum(random_state_mat[k, c("mutL", "mutS", "mutY")]),
          mutL_mutY = 2-sum(random_state_mat[k, c("mutL", "mutY")]),
          mutL_mutS = 2-sum(random_state_mat[k, c("mutL", "mutS")]),
          mutM_mutY = 2-sum(random_state_mat[k, c("mutM", "mutY")]),
          mutL_dnaQ = 2-sum(random_state_mat[k, c("mutL", "dnaQ")])
        )
        if(any(state_check == 0)) {
          multi_kos <- names(state_check)[state_check==0]
          # if("mutL_mutS_mutY" %in% multi_kos) {
          #   multi_kos <- multi_kos[!multi_kos %in% c("mutL_mutS", "mutL_mutY")]
          # } 
          if("mutL_mutS_mutH" %in% multi_kos) {
            multi_kos <- multi_kos[!multi_kos %in% c("mutL_mutS")]
          }
          multi_sums <- sum((scaled_ts_tv_multi[multi_kos]))
          enz_covered <- unique(unlist(sapply(multi_kos, function(x) strsplit(x, split = "_")[[1]])))
          rem_enz <- enz_oi[!enz_oi %in% enz_covered]
          rem_sums <- sum(random_state_mat[k,rem_enz]*(scaled_ts_tv[rem_enz]))
          random_sums[k] <- (rem_sums + multi_sums)
        } else {
          random_sums[k] <- sum(random_state_mat[k,]*(scaled_ts_tv))
        }
      }
      random_sums_list[[ii]] <- random_sums
    }
    
    bias_diff_r[[jj]] <- sapply(random_sums_list, diff)
    ups_downs_r[[jj]] <- sapply(bias_diff_r[[jj]], function(x) x[x != 0])
    
    no_of_ups_r[[jj]] <- sapply(ups_downs_r[[jj]], function(x) sum(x > 0))
    no_of_downs_r[[jj]] <- sapply(ups_downs_r[[jj]], function(x) sum(x < 0))
    
    no_changes_r[j, jj] <- length(which(sapply(bias_diff_r[[jj]], function(x) all(x==0))))
    one_change_r[j, jj] <- length(which(sapply(ups_downs_r[[jj]], length) == 1))
    reinforced_r[j, jj] <- length(which(sapply(ups_downs_r[[jj]], function(x) all(x < 0) | all(x > 0))))
    reinforced_r[j, jj] <- reinforced_r[j, jj] - one_change_r[j, jj] - no_changes_r[j, jj] ##they are considered in reinforced because they are NULL
    
    cons_up_event_r[[jj]] <- vector("numeric", length(ups_downs_r[[jj]]))
    cons_down_event_r[[jj]] <- vector("numeric", length(ups_downs_r[[jj]]))
    for(l in 1:length(ups_downs_r[[jj]])) {
      cons_up_event_r[[jj]][l] <- 0
      cons_down_event_r[[jj]][l] <- 0
      if(length(ups_downs_r[[jj]][[l]])>1) {
        for(ll in 1:(length(ups_downs_r[[jj]][[l]])-1)) {
          if(ups_downs_r[[jj]][[l]][ll]<0 & ups_downs_r[[jj]][[l]][ll+1]>0) {
            cons_up_event_r[[jj]][l] <- cons_up_event_r[[jj]][l]+1
          } else if(ups_downs_r[[jj]][[l]][ll]>0 & ups_downs_r[[jj]][[l]][ll+1]<0) {
            cons_down_event_r[[jj]][l] <- cons_down_event_r[[jj]][l]+1
          }
        }
      }
    }
    total_fluc_r[[j]][[jj]] <- table(cons_up_event_r[[jj]]+cons_down_event_r[[jj]])
  }
  
  end.time <- Sys.time()
  print(end.time-start.time)
}

##GC/AT

no_changes_r_at_gc <- matrix(NA, nrow = 500, ncol = 20)
one_change_r_at_gc <- matrix(NA, nrow = 500, ncol = 20)
reinforced_r_at_gc <- matrix(NA, nrow = 500, ncol = 20)

total_fluc_r_at_gc <- list()

for(j in 1:500) {
  print(j)
  # if(enz_oi[i] == "mutT") next
  start.time <- Sys.time()
  
  bias_diff_r_at_gc <- list()
  ups_downs_r_at_gc <- list()
  
  no_of_ups_r_at_gc <- list()
  no_of_downs_r_at_gc <- list()
  
  cons_up_event_r_at_gc <- list()
  cons_down_event_r_at_gc <- list()
  
  total_fluc_r_at_gc[[j]] <- list()
  
  for(i in 1:length(enz_oi)) {
    print(enz_oi[i])
    test <- sim.history(tree = final_tree, Q = trans_mat[[enz_oi[i]]][[j]], anc = as.character(anc_state[enz_oi[i]]), nsim = 20, message = F)
    random_states <- matrix(NA, nrow = (Ntip(final_tree)+final_tree$Nnode), ncol = length(test), dimnames = list(c("1094", test[[1]]$edge[,2])))
    for(jj in 1:length(test)) {
      random_states[,jj] <- c(as.numeric(test[[jj]]$node.states[1,1]), as.numeric(test[[jj]]$node.states[,2]))
    }
    random_enz_state[[enz_oi[i]]] <- random_states
  }
  
  print("Doing calculations...")
  #calc
  for(jj in 1:length(test)) {
    random_sums_list_at_gc <- list()
    for(ii in 1:length(all_parent_nodes)) {
      random_state_mat <- sapply(enz_oi, function(x) 1-random_enz_state[[x]][,jj][match(as.character(all_parent_nodes[[ii]]), names(random_enz_state[[x]][,jj]))])
      random_state_mat <- rbind(random_state_mat, 1-sapply(enz_oi, function(x) random_enz_state[[x]][,jj][names(random_enz_state[[x]][,jj])==as.character(ii)]))
      rownames(random_state_mat) <- as.character(c(all_parent_nodes[[ii]], ii))
      random_sums <- vector(mode = "numeric", length = nrow(random_state_mat))
      names(random_sums) <- rownames(random_state_mat)
      for(k in 1:nrow(random_state_mat)) {
        state_check <- c(
          mutL_mutS_mutH = 3-sum(random_state_mat[k, c("mutL", "mutS", "mutH")]),
          # mutL_mutS_mutY = 3-sum(random_state_mat[k, c("mutL", "mutS", "mutY")]),
          mutL_mutY = 2-sum(random_state_mat[k, c("mutL", "mutY")]),
          mutL_mutS = 2-sum(random_state_mat[k, c("mutL", "mutS")]),
          mutM_mutY = 2-sum(random_state_mat[k, c("mutM", "mutY")]),
          mutL_dnaQ = 2-sum(random_state_mat[k, c("mutL", "dnaQ")])
        )
        if(any(state_check == 0)) {
          multi_kos <- names(state_check)[state_check==0]
          # if("mutL_mutS_mutY" %in% multi_kos) {
          #   multi_kos <- multi_kos[!multi_kos %in% c("mutL_mutS", "mutL_mutY")]
          # } 
          if("mutL_mutS_mutH" %in% multi_kos) {
            multi_kos <- multi_kos[!multi_kos %in% c("mutL_mutS")]
          }
          multi_sums <- sum((scaled_at_gc_multi[multi_kos]))
          enz_covered <- unique(unlist(sapply(multi_kos, function(x) strsplit(x, split = "_")[[1]])))
          rem_enz <- enz_oi[!enz_oi %in% enz_covered]
          rem_sums <- sum(random_state_mat[k,rem_enz]*(scaled_at_gc[rem_enz]))
          random_sums[k] <- (rem_sums + multi_sums)
        } else {
          random_sums[k] <- sum(random_state_mat[k,]*(scaled_at_gc))
        }
      }
      random_sums_list_at_gc[[ii]] <- random_sums
    }
    
    bias_diff_r_at_gc[[jj]] <- sapply(random_sums_list_at_gc, diff)
    ups_downs_r_at_gc[[jj]] <- sapply(bias_diff_r_at_gc[[jj]], function(x) x[x != 0])
    
    no_of_ups_r_at_gc[[jj]] <- sapply(ups_downs_r_at_gc[[jj]], function(x) sum(x > 0))
    no_of_downs_r_at_gc[[jj]] <- sapply(ups_downs_r_at_gc[[jj]], function(x) sum(x < 0))
    
    no_changes_r_at_gc[j, jj] <- length(which(sapply(bias_diff_r_at_gc[[jj]], function(x) all(x==0))))
    one_change_r_at_gc[j, jj] <- length(which(sapply(ups_downs_r_at_gc[[jj]], length) == 1))
    reinforced_r_at_gc[j, jj] <- length(which(sapply(ups_downs_r_at_gc[[jj]], function(x) all(x < 0) | all(x > 0))))
    reinforced_r_at_gc[j, jj] <- reinforced_r_at_gc[j, jj] - one_change_r_at_gc[j, jj] - no_changes_r_at_gc[j, jj] ##they are considered in reinforced because they are NULL
    
    cons_up_event_r_at_gc[[jj]] <- vector("numeric", length(ups_downs_r_at_gc[[jj]]))
    cons_down_event_r_at_gc[[jj]] <- vector("numeric", length(ups_downs_r_at_gc[[jj]]))
    for(l in 1:length(ups_downs_r_at_gc[[jj]])) {
      cons_up_event_r_at_gc[[jj]][l] <- 0
      cons_down_event_r_at_gc[[jj]][l] <- 0
      if(length(ups_downs_r_at_gc[[jj]][[l]])>1) {
        for(ll in 1:(length(ups_downs_r_at_gc[[jj]][[l]])-1)) {
          if(ups_downs_r_at_gc[[jj]][[l]][ll]<0 & ups_downs_r_at_gc[[jj]][[l]][ll+1]>0) {
            cons_up_event_r_at_gc[[jj]][l] <- cons_up_event_r_at_gc[[jj]][l]+1
          } else if(ups_downs_r_at_gc[[jj]][[l]][ll]>0 & ups_downs_r_at_gc[[jj]][[l]][ll+1]<0) {
            cons_down_event_r_at_gc[[jj]][l] <- cons_down_event_r_at_gc[[jj]][l]+1
          }
        }
      }
    }
    total_fluc_r_at_gc[[j]][[jj]] <- table(cons_up_event_r_at_gc[[jj]]+cons_down_event_r_at_gc[[jj]])
  }
  
  end.time <- Sys.time()
  print(end.time-start.time)
}

###Supplementary figure - barplots for the randomisations####
full_list <- unlist(total_fluc_r)
full_table <- NULL
sd_table <- NULL
t_tests <- list()
for(i in as.character(1:10)) {
  full_table[i] <- sum(full_list[names(full_list)==i])/10000
  sd_table[i] <- sd(full_list[names(full_list)==i], na.rm = T)
  if(!i %in% as.character(8:10)) {
    t_tests[[i]] <- t.test(full_list[names(full_list)==i], mu = table(total_fluc)[i])
  }
}


full_list_at_gc <- unlist(total_fluc_r_at_gc)
full_table_at_gc <- NULL
sd_table_at_gc <- NULL
t_tests_at_gc <- list()
for(i in as.character(1:11)) {
  full_table_at_gc[i] <- sum(full_list_at_gc[names(full_list_at_gc)==i])/10000
  sd_table_at_gc[i] <- sd(full_list_at_gc[names(full_list_at_gc)==i], na.rm = T)
  if(!i %in% as.character(7:11)) {
    t_tests_at_gc[[i]] <- t.test(full_list_at_gc[names(full_list_at_gc)==i], mu = table(total_fluc_at_gc)[i])
  }
}


##random events barplot

ts_table_r <- tibble(
  category = c("No change", "Single change", "Reinforced", "Reversals"),
  `events` = c(mean(as.vector(no_changes_r)), mean(as.vector(one_change_r)), mean(as.vector(reinforced_r)), 0),
  `1` = c(0, 0, 0, full_table["1"]),
  `2` = c(0, 0, 0, full_table["2"]),
  `3` = c(0, 0, 0, full_table["3"]),
  `4` = c(0, 0, 0, full_table["4"]),
  `5` = c(0, 0, 0, full_table["5"]),
  `6` = c(0, 0, 0, full_table["6"]),
  `7` = c(0, 0, 0, full_table["7"]),
  `8` = c(0, 0, 0, full_table["8"]),
  `9` = c(0, 0, 0, full_table["9"]),
  `10` = c(0, 0, 0, full_table["10"]),
  `11` = c(0, 0, 0, full_table["11"]),
  # `12` = c(0, 0, 0, full_table["12"]),
  bias = "Tv/Ts bias"
)

gc_table_r <- tibble(
  category = c("No change", "Single change", "Reinforced", "Reversals"),
  `events` = c(mean(as.vector(no_changes_r_at_gc)), mean(as.vector(one_change_r_at_gc)), mean(as.vector(reinforced_r_at_gc)), 0),
  `1` = c(0, 0, 0, full_table_at_gc["1"]),
  `2` = c(0, 0, 0, full_table_at_gc["2"]),
  `3` = c(0, 0, 0, full_table_at_gc["3"]),
  `4` = c(0, 0, 0, full_table_at_gc["4"]),
  `5` = c(0, 0, 0, full_table_at_gc["5"]),
  `6` = c(0, 0, 0, full_table_at_gc["6"]),
  `7` = c(0, 0, 0, full_table_at_gc["7"]),
  `8` = c(0, 0, 0, full_table_at_gc["8"]),
  `9` = c(0, 0, 0, full_table_at_gc["9"]),
  `10` = c(0, 0, 0, full_table_at_gc["10"]),
  `11` = c(0, 0, 0, full_table_at_gc["11"]),
  # `12` = c(0, 0, 0, full_table_at_gc["12"]),
  bias = "GC/AT bias"
)

bcols <- c("#8DA0CB", rev(colorRampPalette(c("rosybrown1", "red3"))(11)))

ts_btable_r <- ts_table_r %>% 
  reshape2::melt() %>% 
  mutate(sd = c(sd(as.vector(no_changes_r)), sd(as.vector(one_change_r)), sd(as.vector(reinforced_r)), 0, rep(sd_table[as.character(1:11)], each = 4)), sd_posn = c(value[1:4], cumsum(value[variable!="events"]))) %>% 
  mutate(sd = if_else(value == 0, 0, sd), sd_posn = if_else(value == 0, 0, sd_posn)) %>% 
  mutate(variable = factor(variable, levels = c("events", "11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1"))) 

gc_btable_r <- gc_table_r %>% 
  reshape2::melt() %>% 
  mutate(sd = c(sd(as.vector(no_changes_r_at_gc)), sd(as.vector(one_change_r_at_gc)), sd(as.vector(reinforced_r_at_gc)), 0, rep(sd_table_at_gc[as.character(1:11)], each = 4)), sd_posn = c(value[1:4], cumsum(value[variable!="events"]))) %>% 
  mutate(sd = if_else(value == 0, 0, sd), sd_posn = if_else(value == 0, 0, sd_posn)) %>% 
  mutate(variable = factor(variable, levels = c("events", "11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1"))) 

btable_r <- rbind(ts_btable_r, gc_btable_r)
btable_r$category <- factor(btable_r$category, levels = c("No change", "Single change", "Reinforced", "Reversals"))
btable_r$bias <- factor(btable_r$bias, levels = c("Tv/Ts bias", "GC/AT bias"))

bplot_r <- btable_r %>% 
  ggplot(aes(category, value, fill = variable)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = sd_posn-sd, ymax = sd_posn+sd), color = "gray50", width = 0.1, position = "identity") +
  # scale_fill_brewer(palette = "Set2", direction = -1) +
  scale_fill_manual(values = bcols) +
  facet_wrap(~bias, nrow = 2) +
  theme_cowplot(font_size = 22) +
  scale_y_continuous(expand = c(0,0), limits = c(-2.5, 1093)) +
  labs(
    x = "",
    y = "Number of lineages",
    fill = ""
  ) +
  # theme(axis.text.y = element_text(angle = 45, hjust = 1)) + 
  coord_flip()

bplot_r