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
modulos = cluster_optimal(gph) #Lo hago optimal como en el paper de Guimera (2005).
