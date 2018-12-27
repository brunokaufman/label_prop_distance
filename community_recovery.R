source('bipartite.v2.R')

distance.matrices = function(g){
	#Grabo la matriz de adyacencia de la red y otras cosas utiles.
	adyacencia_original = as.matrix(as_adjacency_matrix(g))
	sq_adyacencia_original = (t(adyacencia_original) %*% adyacencia_original)#Adyacencia por traspuesta.
	degs = degree(g)
	degmat = degs %o% degs
	degmat[which(sq_adyacencia_original == 0)] == 0

	#Calculo la matriz de adyacencia basada en la medida de similaridad coseno.
	sq_row_sums = rowSums(as.array(adyacencia_original) ** 2)
	adyacencia_coseno = sq_adyacencia_original/sqrt(degs %o% degs)

	#Calculo la matriz de adyacencia basada en la similaridad Pearson.
	desvio_vecinos_promedio = sq_adyacencia_original - (degmat)/length(degs)
	diags_desv = diag(as.matrix(desvio_vecinos_promedio))
	adyacencia_pearson = desvio_vecinos_promedio / sqrt(diags_desv %o% diags_desv)

	min_allowed = min(min(adyacencia_pearson), 0)
	adyacencia_pearson = (adyacencia_pearson - min_allowed) * (max(adyacencia_pearson) - min_allowed) #No me gusta tener negativos asi que reacomodo para tener entre 0 y max si llego a tener negativos.

	#Calculo la matriz de adyacencia basada en la topological overlap.
	vector_unos = rep(1, length(degs)) #Vector donde toda entrada es un uno.
	matriz_fila_ks = degs %o% as.array(vector_unos) #Matriz donde cada fila tiene el mismo ki
	matriz_min_ks = min(matriz_fila_ks, t(matriz_fila_ks)) #min(ki, kj) lo hago matricialmente.
	numerador_overlap = sq_adyacencia_original + adyacencia_original
	denominador_overlap = matriz_min_ks - adyacencia_original + 1
	adyacencia_topological_overlap = numerador_overlap / denominador_overlap

	list(adyacencia_original, adyacencia_coseno, adyacencia_pearson, adyacencia_topological_overlap)
}

dist_list = distance.matrices(g_users)
lapname = c('Original', 'Coseno', 'Pearson', 'Topoverlap')

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

rank.recommendations = function(lap, seeds, verbose=T){
	df_rankings = matrix(nrow=nrow(lap), ncol=length(seeds))
	colnames(df_rankings) = seeds
	for(seed in seeds){
		if(verbose){
			print(paste('Seed actual:', seed))
		}
		seed_vector = rep(0, nrow(lap))
		names(seed_vector) = V(g_users)$name
		seed_vector[seed] = 1
		ranking = rownames(lap)[order(label.propagation(lap, seed_vector, alpha = alpha_propagate), decreasing=T)]
		df_rankings[,seed] = ranking
	}
	rownames(df_rankings) = 1:length(V(g_users))
	colnames(df_rankings) = V(g_users)$name
	return(df_rankings)
}

lambda_zhou = 0.5
alpha_propagate = 0.3
lap_list = lapply(dist_list, lap.zhou, lambda=lambda_zhou)

ranks_list = lapply(lap_list, rank.recommendations, seeds=V(g_users))

simple.auc = function(TPR, FPR){
	dFPR <- c(diff(FPR), 0)
	dTPR <- c(diff(TPR), 0)
	sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

roc.curve = function(ranked_nodes, true_nodes){
	matches = ranked_nodes %in% true_nodes
	data.frame(TPR=cumsum(matches)/sum(matches), FPR=cumsum(!matches)/sum(!matches), matches)
}

seed.aucs = function(ranks_list){
    aucs = rep(0, ncol(ranks_list))
    for(i in 1:ncol(ranks_list)){
        in_community = V(g_users)$color == V(g_users)$color[i]
        in_community = V(g_users)$name[which(in_community)]
        roc = roc.curve(ranks_list[,i], in_community)
        aucs[i] = simple.auc(roc$TPR, roc$FPR)
    }
    return(aucs)
}

aucs_dists = lapply(ranks_list, FUN=seed.aucs)

#Ahora calculo los módulos de alguna forma.
#No puedo usar cluster_optimal() Guimera (2005) porque CRAN ya no supportea GLPK.
membership = match(V(g_users)$color, unique(V(g_users)$color)) #Uso la particion que ya conozco.
horizontal_member_mat = as.array(rep(1, length(V(g_users)))) %o% membership #La matriz tiene las memberships en direccion horizontal y las repite en cada fila.
comembership_mat = as.matrix(horizontal_member_mat == t(horizontal_member_mat)) #Si lo comparo con su transpuesta hago una matriz donde la componente (i,j) me dice si i y j pertenecen al mismo modulo.
comembership_mat[which(comembership_mat)] = 1 #Paso los booleanos a valores numericos.
comembership_mat[which(!comembership_mat)] = 0

#Calculo las cantidades que usa Guimera (2005) de forma matricial.
adyacencia = as.matrix(as_adjacency_matrix(g_users))
k_i = rowSums(adyacencia * comembership_mat) #El k_i es el grado del nodo i dentro de su mismo modulo.
k_medio_si = (comembership_mat %*% k_i) / rowSums(comembership_mat) #Calculo la media del k_i para el modulo al cual pertenece cada nodo.
k_var_si = (comembership_mat %*% ((k_i - k_medio_si) ** 2)) / rowSums(comembership_mat) #Lo mismo con la varianza.
z_i = as.numeric((k_i - k_medio_si) / sqrt(k_var_si)) #Finalmente obtengo el z-score (que tan bien conectado esta a su modulo).
z_i[which(is.nan(z_i))] = 0 #Los que tenian varianza cero por ser todos identicos los nodos en el cluster, les pongo cero z-score.

#Ahora planteo una matriz que mapee del nodo i al modulo m.
modu_map = matrix(nrow=length(unique(membership)), ncol=length(V(g_users)))
modu_map[,] = 0 #Seteo default todo a cero.
for(i in 1:length(V(g_users))){
    modu_map[membership[i],i] = 1 #Senalo a que modulo pertenece el nodo i.
}
module_links = modu_map %*% adyacencia #Si mapeo la adyacencia por los modulos tengo la cantidad de links del nodo i a cada modulo.

p_i = 1 - (colSums(module_links ** 2) / (degree(g_users) ** 2)) #Calculo el p-score de Guimera (2005).

df_topo_auc = data.frame(Z.SC=z_i, P.SC=p_i, AUC.ORI=aucs_dists[[1]], AUC.COS=aucs_dists[[2]], AUC.PEA=aucs_dists[[3]], AUC.TOP=aucs_dists[[4]])

valid_z = which(df_topo_auc[,1] != 0)
valid_p = which(df_topo_auc[,2] != 0)

df_topo_auc[3:6] = log(1 - df_topo_auc[3:6])#Hago log la escala del auc alrededor del AUC=1.

colores = c('black', 'magenta', 'blue', 'orange')

agg.plot = function(x, y, xname='', yname='', points=F, pch=1, col='black'){
	agg = aggregate(x, by=list(y), FUN=mean)
	agg = aggregate(agg[,1], by=list(agg[,2]), FUN=mean) #Hago el aggregate en el otro sentido tambien.
	if(points){
        points(agg[,1], agg[,2], pch=pch, col=col)
    }
    else{
        plot(agg[,1], agg[,2], xlab=xname, ylab=yname, pch=pch, col=col)
    }
}

png('z_auc.png')
agg.plot(df_topo_auc$Z.SC[valid_z], df_topo_auc[valid_z, 3], pch=1, col=colores[1], xname='Z-score de nodo', yname='AUC de nodo como semilla')
agg.plot(df_topo_auc$Z.SC[valid_z], df_topo_auc[valid_z, 4], pch=2, col=colores[2], points=T)
agg.plot(df_topo_auc$Z.SC[valid_z], df_topo_auc[valid_z, 5], pch=3, col=colores[3], points=T)
agg.plot(df_topo_auc$Z.SC[valid_z], df_topo_auc[valid_z, 6], pch=4, col=colores[4], points=T)
legend(0, -0.2, legend=lapname, pch=1:4, col=colores)
dev.off()

png('p_auc.png')
agg.plot(df_topo_auc$P.SC[valid_p], df_topo_auc[valid_p, 3], pch=1, col=colores[1], xname='P-score de nodo', yname='AUC de nodo como semilla')
agg.plot(df_topo_auc$P.SC[valid_p], df_topo_auc[valid_p, 4], pch=2, col=colores[2], points=T)
agg.plot(df_topo_auc$P.SC[valid_p], df_topo_auc[valid_p, 5], pch=3, col=colores[3], points=T)
agg.plot(df_topo_auc$P.SC[valid_p], df_topo_auc[valid_p, 6], pch=4, col=colores[4], points=T)
legend(0.1, -0.5, legend=lapname, pch=1:4, col=colores)
dev.off()
