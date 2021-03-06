lap.zhou = function(adyacencia, lambda){
	strengths = rowSums(as.array(adyacencia))
	degs = diag(t(adyacencia) %*% adyacencia)
	lap = (diag(strengths) - adyacencia) / ((degs %o% degs) ** lambda)
	return(lap)
}

label.propagation <- function(lap, seeds, alpha=0.85, iterations=100, h=0.1){#Minimo sugiero h * iterations = 10.
	res = rep(1.0/nrow(lap), nrow(lap)) #Inicializo como todo igual y normalizado.

	for(i in 1:iterations){
		res = res + h * (alpha * (seeds - res) - (1.0 - alpha) * (lap %*% res))
	}

	return(res)
}

rank.recommendations = function(lap, seeds_ranked, verbose=T){
	df_rankings = matrix(nrow=nrow(lap), ncol=seeds_considered)
	colnames(df_rankings) = seeds_ranked[1:seeds_considered]
	for(seed in seeds_ranked[1:seeds_considered]){
		if(verbose){
			print(paste('Seed actual:', seed))
		}
		seed_vector = rep(0, nrow(lap))
		names(seed_vector) = V(gph)$name
		seed_vector[seed] = 1
		ranking = rownames(lap)[order(label.propagation(lap, seed_vector, alpha = alpha_propagate), decreasing=T)]
		df_rankings[,seed] = ranking
	}
	return(df_rankings)
}

lap_list = lapply(dist_list, lap.zhou, lambda=lambda_zhou)

nodos_removed = c(as.character(df_removed[,1]), as.character(df_removed[,2]))
nodos_removed = nodos_removed[which(nodos_removed %in% V(gph)$name)] #Me quedo solo con los que encuentre en la red que conozco.
times_removed = aggregate(rep(1, length(nodos_removed)), by=list(nodos_removed), FUN=sum)
seeds_ranked = times_removed[order(times_removed[,2], decreasing=T),1]

seeds_considered = min(seeds_considered, length(seeds_ranked)) #Si tengo menos de lo que queria, me quedo con eso.

ranks_list = lapply(lap_list, rank.recommendations, seeds_ranked=seeds_ranked[1:seeds_considered])
