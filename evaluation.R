simple.auc = function(TPR, FPR){
	dFPR <- c(diff(FPR), 0)
	dTPR <- c(diff(TPR), 0)
	sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

roc.curve = function(ranked_nodes, true_nodes){
	matches = ranked_nodes %in% true_nodes
	data.frame(TPR=cumsum(matches)/sum(matches), FPR=cumsum(!matches)/sum(!matches), matches)
}

seeds.auc = function(df_ranking_list, df_removed, df_kept){
	auc.matrix = matrix(nrow = ncol(df_ranking_list[[1]]), ncol=length(df_ranking_list))
	rownames(auc.matrix) = colnames(df_ranking_list[[1]])
	colnames(auc.matrix) = lapname
	for(i in 1:length(df_ranking_list)){
		df_ranking = df_ranking_list[[i]]
		for(seed in colnames(df_ranking)){
			nodematch_1 = df_removed[which(df_removed[,1] == seed), 2]
			nodematch_2 = df_removed[which(df_removed[,2] == seed), 1]
			removed_linked_nodes = unique(c(as.character(nodematch_1), as.character(nodematch_2)))
			all_nodes = unique(df_ranking[, seed]) #Agarro todos los nodos considerados en el ranking.

			keptmatch_1 = df_kept[which(df_kept[,1] == seed), 2]
			keptmatch_2 = df_kept[which(df_kept[,2] == seed), 1]
			kept_linked_nodes = unique(c(as.character(keptmatch_1), as.character(keptmatch_2)))
			all_nodes = all_nodes[which(!(all_nodes %in% c(seed, kept_linked_nodes)))] #Obviamente quito los que ya se que se conectan al seed, y el seed mismo.

			seed_roc = roc.curve(all_nodes, removed_linked_nodes) #Aprovecho que estan en orden.
			auc.matrix[seed, i] = simple.auc(seed_roc$TPR, seed_roc$FPR)
		}
	}

	return(auc.matrix)
}

auc.mat = seeds.auc(ranks_list, df_removed, df_kept)
print(auc.mat)
