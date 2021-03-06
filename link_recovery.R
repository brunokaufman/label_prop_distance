#Primero cargo la red y quito alguna cantidad de links para luego probar recuperacion.
print('Cargando y configurando red.')
filename = './data/yeast_Y2H.txt' #El nombre del archivo.
sep='\t' #El separador del archivo.
prop_remove = 0.4 #La proporcion de links que quiero sacar.
source('setup_remove.R')
print('')

alphas = seq(from=0.1, to=0.9, by=0.1)
mat_means = matrix(nrow=length(alphas), ncol=length(lapname))
colnames(mat_means) = lapname
rownames(mat_means) = alphas
for(alpha in alphas){
	#Despues propago labels.
	print('Propagando labels.')
	lambda_zhou = 0.5 #Lambda del Laplaciano.
	alpha_propagate = alpha #Constante de damping de la propagacion.
	seeds_considered = 10 #Cantidad de seeds a evaluar.
	source('label_propagation.R')
	print('')

	#Despues evaluo por el AUC de las ROC.
	print('Evaluando resultados.')
	source('evaluation.R')

	means = as.numeric(gsub(as.character(summary(auc.mat)[4,]), pattern='Mean| |:', replacement=''))
	mat_means[which(rownames(mat_means) == alpha), ] = means
	print('')
}

#Grafico las tendencias encontradas para los AUCs recuperados para las distintas nociones de distancia en funcion del restart alpha de la propagacion.
png('grafico_aucs_alpha.png')
plot(rownames(mat_means), mat_means[,lapname[1]], xlim=c(0,1), ylim=c(0,1), xlab='Alpha restart', ylab='AUC de links recuperados', col='black', type='l', lwd=2)
lines(rownames(mat_means), mat_means[,lapname[2]], xlim=c(0,1), ylim=c(0,1), col='red', lwd=2)
lines(rownames(mat_means), mat_means[,lapname[3]], xlim=c(0,1), ylim=c(0,1), col='green', lwd=2)
lines(rownames(mat_means), mat_means[,lapname[4]], xlim=c(0,1), ylim=c(0,1), col='blue', lwd=2)
legend(0, 0.3, legend=lapname, col=c('black', 'red', 'green', 'blue'), lwd=2)
dev.off()

##########

#A ver que onda la adyacencia de Pearson en alpha=0.1 particularmente. Igual se va a calcular todo.
print('Evaluando todo en alpha=0.1.')
lambda_zhou = 0.5 #Lambda del Laplaciano.
alpha_propagate = 0.1 #Constante de damping de la propagacion.
seeds_considered = 100 #Cantidad de seeds a evaluar.
source('label_propagation.R')
print('')

#Despues evaluo por el AUC de las ROC.
print('Evaluando resultados.')
source('evaluation.R')

png('hist_auc_original.png')
hist(auc.mat[,1], freq=T, xlim=c(0.9, 1.0), xlab='AUC recuperado originalmente', breaks = 75)
dev.off()

png('hist_auc_coseno.png')
hist(auc.mat[,2], freq=T, xlim=c(0.9, 1), xlab='AUC recuperado por coseno', breaks = 75)
dev.off()

png('hist_auc_pearson.png')
hist(auc.mat[,3], freq=T, xlim=c(0.1, 0.5), xlab='AUC recuperado por Pearson', breaks = 75)
dev.off()

png('hist_auc_topoverlap.png')
hist(auc.mat[,4], freq=T, xlim=c(0.9, 1), xlab='AUC recuperado por topological overlap', breaks = 75)
dev.off()

