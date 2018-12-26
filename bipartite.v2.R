options(formatR.arrow=TRUE,width=90)
library(igraph)
library(poweRlaw)
library(Biobase)

nObjects = 100
nUsers   = 100
#distribucion de grado power-law para ku
alpha = 2 #param exponente power law para P(ku)
kmin  = 1 #param kmin para P(ku)
#clusters de usuarios
ngu = 3     #numero de grupos de usuarios
pin = 0.85 #probabilidad cableado intra
#clusters de objetos
kShrinkMax   = 0.05   # tamanio de cluster object = pin * sum(ku) * lShrink
kShrinkMin   = 0.05   #  kShrink = unif(entre  kShrinkMin y kShrinkMax)

##################################################################################################

# in
#  nObjects:       number of objects
#  nUsers  :       number of users
#  userCommunities: vector of labels. if NULL users are radomly assigned to size blanced clusters
#  

buildBipartiteGraph=function(userCommunities=NULL,
                              userDegrees=NULL,
                              nUsers,nObjects,
                              alpha,kmin,ngu,pin,kShrinkMin,kShrinkMax){
 
 if(!is.null(userCommunities)){
  if(nUsers!=length(userCommunities)) warning(paste("nUsers:",nUsers," not equal to length(userCommunities):",length(userCommunities),"\n"))
  ulabel = userCommunities
 }else{
  #asignacion de clusters de usuarios
  ulabel=1 + order(runif(nUsers))%%ngu  #clusters de masa uniforme  
 }

 # distribucion degree power law para usuarios
 if(is.null(userDegrees)){
  ku = rpldis(nUsers,kmin,alpha)
  ku = pmin(ku,nObjects)         #nadie puede comprar mas que la cantidad original de objetos
 }else{
  if(nUsers!=length(userDegrees)) warning(paste("nUsers:",nUsers," not equal to length(userDegrees):",length(userDegrees),"\n")) 
  ku = userDegrees
 }
 
 #tamanio de comunidades
 facShrink=runif(1,min=kShrinkMin,kShrinkMax)
 sizeClusters=aggregate(ku,list(ulabel),sum)
 sizeClusters=cbind(sizeClusters,round(facShrink*sizeClusters[,2]*pin))
 colnames(sizeClusters)=c("community","userClusterDegree","objectClusterSize")
 #print(xtable(sizeClusters,caption="Size of the users communities and their associates object clusters "),include.rownames=FALSE)
 
 #itero por clusters
 uulabel=unique(ulabel)
 res=c()
 lcommunities=list()
 for(i in seq_along(uulabel)){
  iu   = which(ulabel== uulabel[i])

  #dado un cluster de usuarios me fijo que grado total tiene
  # y elijo ktot * pin objetos
  # de un universo de overlapObjectClusters*n
                     
  objsIn = sample(1:nObjects,round(facShrink*pin*sum(ku[iu])))
  objsOut= (1:nObjects)
  objsOut=objsOut[!objsOut%in%objsIn]
 
  prob=c(rep(pin/length(objsIn),length(objsIn)), rep((1-pin)/length(objsOut),length(objsOut))) 
  for(j in seq_along(iu)){
   res=rbind(res,matrix(c(rep(iu[j],ku[iu[j]]),sample(c(objsIn,objsOut),ku[iu[j]],prob=prob,replace=FALSE)),ncol=2))
  }
 
  lcommunities[[i]]              = list()
  lcommunities[[i]][["users"]]   = iu
  lcommunities[[i]][["objects"]] = objsIn
 }
 colnames(res)=c("user","object")
 
 a    = res
 a[,2]= a[,2]+nUsers
 g = graph.bipartite(type=c(rep(1,nUsers),rep(0,nObjects)), 
                      edges=as.vector(as.matrix(t(a))),directed=FALSE)
 V(g)$name = c(paste("U",1:nUsers,sep=""),paste("O",1:nObjects,sep=""))
 V(g)$color           = rgb(.8,.8,.8,.5)
 V(g)$color[1:nUsers] = rainbow(length(uulabel),alpha=0.5)[ulabel]
 V(g)$shape = "circle"
 V(g)$shape[1:nUsers] = "square"
 V(g)$label.color = rgb(0,.2,.2,.5)
 V(g)$label.cex   = .48
 V(g)$frame.color = rgb(0.5,.5,.5,.5)
 V(g)$size          =8
 V(g)$size[1:nUsers]=10
 E(g)$color = rgb(0.5,.5,.5,.6)
 
 
 return(list(graph=g,lcommunities=lcommunities,ku=ku,ulabel=ulabel))
}

plotAdjacency=function(g,lcommunities){
  nusers  =sum(V(g)$type)
  nobjects=sum(!V(g)$type)
  
  laux       = lapply(lcommunities,function(x){return(x$users)})
  names(laux)= seq_along(lcommunities)
  ulabel     = as.numeric(unlist(reverseSplit(laux)))

  laux=lapply(lcommunities,function(x){return(x$objects)})
  names(laux)=seq_along(lcommunities)
  lObjectInUCluster=reverseSplit(laux)
  
  a=get.adjacency(g,sparse=FALSE)[1:nusers,nusers+(1:nobjects)]

  #ordeno columnas
  if(FALSE){
   dobjs=dist(t(a),method="binary")
   hc=hclust(dobjs)
   oorder=hc$order
  }else{
   aux=unique(unlist(lapply(lcommunities,function(x){return(x$objects)})))
   oorder=1:nObjects
   oorder=c(aux,oorder[!oorder%in%aux])
  }
  image(a[order(ulabel),oorder],col=c("wheat","black"),axes=FALSE,xlab="users",ylab="objects",main="adjacency matrix")
  for(i in 1:ngu){           #color para usuarios
   j=range(which(sort(ulabel)==i))/length(ulabel)
   lines(j,c(0,0)/length(ulabel),col=i+1,lwd=6)
  }

  ccol=rep("1",nObjects)     #color para objetos
  aux=unlist(lapply(lObjectInUCluster,function(x){return(x[1])}))
  ccol[as.numeric(names(aux))]=aux
  ccol=as.numeric(ccol)
  ccol[as.numeric(names(aux))]=ccol[as.numeric(names(aux))]+1
  ccol=ccol[oorder]
  points(rep(0,nObjects),seq(0,1,length=nObjects),col=ccol,pch=15)
 }

#proyecto red segun Zhou
if(FALSE){
projectBip=function(g,type=TRUE){
 if(type){
  a=t(get.incidence(g))
 }else{
  a=get.incidence(g)
 }
 
 ku =apply(a,1,sum)
 ko =apply(a,2,sum)
 
 a=a[which(ku>0),which(ko>0)]

 #normalizo por fila a_i,alpha / k_i
 a1=a/ku[which(ku>0)]

 #normalizo por columna a_i,beta / k_beta
 a2=a/ko[which(ko>0)]

 w=t(a1)%*%a2

 gw=graph.adjacency(w,mode="directed",weighted=TRUE)
}
}
##################################################################################################


set.seed(12347)

a=buildBipartiteGraph(NULL,NULL,nObjects,nUsers,alpha,kmin,ngu,pin,kShrinkMin,kShrinkMax)

lcommunities=a$lcommunities
g           =a$graph

gg = induced.subgraph(g, unique(as.vector(get.edgelist(g))))

