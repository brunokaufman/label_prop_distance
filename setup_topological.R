library(igraph)

df_fullgraph = read.table(filename, sep='\t', quote='')

gph = graph_from_edgelist(as.matrix(df_fullgraph), directed=F)

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

dist_list = distance.matrices(gph)
lapname = c('Original', 'Coseno', 'Pearson', 'Topoverlap')

