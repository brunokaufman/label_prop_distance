#Primero cargo la red y calculo las medidas de distancia.
print('Cargando y configurando redes.')
filename = './data/yeast_Y2H.txt' #El nombre del archivo.
sep='\t' #El separador del archivo.
source('setup_topological.R')
print('')

agg.plot.loglog = function(x, y, xname='', yname=''){
	agg = aggregate(x, by=list(y), FUN=mean)
	agg = aggregate(agg[,1], by=list(agg[,2]), FUN=mean) #Hago el aggregate en el otro sentido tambien.
	plot(agg[,1], agg[,2], log='xy', xlab=xname, ylab=yname)
}

hist.plot.loglog = function(values, varname=''){
	max_value = max(values)
	min_value = min(values)
	n_breaks = 10
	log_breaks = min_value + c(0, exp((log(max_value - min_value)) * (1:n_breaks) / n_breaks))
	histo = hist(values, breaks=log_breaks, plot=F)

	#Ploteo.
	plot(histo$mids, histo$density, xlab=varname, ylab='Densidad de probabilidad', log='xy')

	#Fiteo y ploteo fit.
	yfit = log(histo$density)
	xfit = log(histo$mids)
	which_finite = which(is.finite(yfit))
	fit = lm(yfit[which_finite] ~ xfit[which_finite])
	xdraw = c(log(min_value), log(max_value))
	ydraw = as.numeric(c(fit$coefficients[2] * log(min_value) + fit$coefficients[1], fit$coefficients[2] * log(max_value) + fit$coefficients[1]))
	lines(exp(xdraw), exp(ydraw), col='red', lwd=2)

	#Escribo la relacion lineal.
	sign_intercept_char = '+'
	if(fit$coefficients[1] < 0){
		sign_intercept_char = '-'
	}
	legend(histo$mids[round(n_breaks * 0.67, 1)], max(histo$density), legend=paste('y =', round(fit$coefficients[2], 2), 'x', sign_intercept_char, round(abs(fit$coefficients[1]), 2)), lwd=2, col='red')
}

##########

#A ver como se comparan entre si las medidas de distancia.

#Histogramas de las distribuciones de strength frente a distintas medidas.
for(i in 1:length(dist_list)){
	strength = rowSums(dist_list[[i]])

	png(paste('deg_dist_', lapname[i], '.png', sep=''))
	hist.plot.loglog(strength, varname=paste('Strength', lapname[i]))
	dev.off()
}

#Paso cada distancia entre par de nodos (sin importar posicion) a un vector y hago una columna para cada forma de ver la distancia.
full_unlist_distances = unlist(dist_list)
unlisted_distances = matrix(nrow=length(dist_list[[1]]), ncol=length(dist_list))
colnames(unlisted_distances) = lapname
for(i in 1:length(dist_list)){
	start_idx = 1 + (i - 1) * length(dist_list[[i]])
	end_idx = i * length(dist_list[[i]])
	unlisted_distances[,i] = full_unlist_distances[start_idx:end_idx]
}
df_dist = data.frame(unlisted_distances)

#Ploteo correlaciones.
png(paste('comparison_', lapname[2], '_', lapname[3], '.png', sep=''))
interesting = which(df_dist[,2] != 0 & df_dist[,3] != 0)
agg.plot.loglog(df_dist[interesting,2], df_dist[interesting,3], xname=lapname[2], yname=lapname[3])
dev.off()

png(paste('comparison_', lapname[2], '_', lapname[4], '.png', sep=''))
interesting = which(df_dist[,2] != 0 & df_dist[,4] != 0)
agg.plot.loglog(df_dist[interesting,2], df_dist[interesting,4], xname=lapname[2], yname=lapname[4])
dev.off()

png(paste('comparison_', lapname[4], '_', lapname[3], '.png', sep=''))
interesting = which(df_dist[,4] != 0 & df_dist[,3] != 0)
agg.plot.loglog(df_dist[interesting,4], df_dist[interesting,3], xname=lapname[4], yname=lapname[3])
dev.off()

#Ahora calculo los módulos de alguna forma.
#No puedo usar cluster_optimal() Guimera (2005) porque CRAN ya no supportea GLPK.
modules = cluster_label_prop(gph) #Uso label propagation como es la familia de algoritmos que despues vamos a usar para recomendar.
horizontal_member_mat = as.array(rep(1, length(V(gph)))) %o% modules$membership #La matriz tiene las memberships en direccion horizontal y las repite en cada fila.
comembership_mat = as.matrix(horizontal_member_mat == t(horizontal_member_mat)) #Si lo comparo con su transpuesta hago una matriz donde la componente (i,j) me dice si i y j pertenecen al mismo modulo.
comembership_mat[which(comembership_mat)] = 1 #Paso los booleanos a valores numericos.
comembership_mat[which(!comembership_mat)] = 0

#Calculo las cantidades que usa Guimera (2005) de forma matricial.
adyacencia = as.matrix(as_adjacency_matrix(gph))
k_i = rowSums(adyacencia * comembership_mat) #El k_i es el grado del nodo i dentro de su mismo modulo.
k_medio_si = (comembership_mat %*% k_i) / rowSums(comembership_mat) #Calculo la media del k_i para el modulo al cual pertenece cada nodo.
k_var_si = (comembership_mat %*% ((k_i - k_medio_si) ** 2)) / rowSums(comembership_mat) #Lo mismo con la varianza.
z_i = as.numeric((k_i - k_medio_si) / sqrt(k_var_si)) #Finalmente obtengo el z-score (que tan bien conectado esta a su modulo).
z_i[which(is.nan(z_i))] = 0 #Los que tenian varianza cero por ser todos identicos los nodos en el cluster, les pongo cero z-score.

#Ploteo distribuciones de z-score.
png('z_score_positivos.png')
hist.plot.loglog(z_i[which(z_i > 0)], varname='Magnitud de z-score positivo')
dev.off()

png('z_score_negativos.png')
hist.plot.loglog(abs(z_i[which(z_i < 0)]), varname='Magnitud de z-score negativo')
dev.off()

#Ahora planteo una matriz que mapee del nodo i al modulo m.
modu_map = matrix(nrow=length(unique(modules$membership)), ncol=length(V(gph)))
modu_map[,] = 0 #Seteo default todo a cero.
for(i in 1:length(V(gph))){
    modu_map[modules$membership[i],i] = 1 #Senalo a que modulo pertenece el nodo i.
}
module_links = modu_map %*% adyacencia #Si mapeo la adyacencia por los modulos tengo la cantidad de links del nodo i a cada modulo.

p_i = 1 - colSums(module_links ** 2) / (degree(gph) ** 2) #Calculo el p-score de Guimera (2005).
png('p_score_nocero.png')
hist(p_i[which(p_i > 0)], density=T, xlab='p-score')
dev.off()

#Ahora quiero calcular una matriz de distancias en el espacio p-z.
#Diferencias unidimensionales
z_distance_mat = outer(z_i, z_i, FUN='-')
p_distance_mat = outer(p_i, p_i, FUN='-')

#Ahora la distancia Euclideana en dos dimensiones.
zp_distance_mat = sqrt(z_distance_mat ** 2 + unname(p_distance_mat) ** 2)

i = 3
x = as.numeric(abs(z_distance_mat))
y = as.numeric(dist_list[[i]])
nonzero = (x != 0 & y != 0)
x = x[nonzero]
y = y[nonzero]

agg.plot.loglog(x, y)
