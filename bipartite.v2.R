options(formatR.arrow=TRUE,width=90)
library(igraph)
library(poweRlaw)
#library(Biobase)

nObjects = 200
nUsers   = 200

#Distribucion de grado ley de potencia discreta para los usuarios.
alpha = 2#Exponente.
kmin  = 1#Grado minimo.

#Clusters de usuarios.
ngu = 4#Cantidad de grupos de usuarios.
pin = 0.7#Probabilidad de que un link sea interno en un grupo.

#Clusters de objetos.
kShrinkMax = 0.05#El tamano de un cluster de objetos = pin * sum(ku) * lShrink.
kShrinkMin = 0.05#kShrink es una variable aleatoria uniforme entre  kShrinkMin y kShrinkMax.

build.bipartite.graph = function(userCommunities=NULL, userDegrees=NULL, nUsers, nObjects, alpha, kmin, ngu, pin, kShrinkMin, kShrinkMax){
    if(!is.null(userCommunities)){
        if(nUsers!=length(userCommunities)){
            warning(paste("nUsers:", nUsers, " not equal to length(userCommunities):", length(userCommunities), "\n"))
        }
    ulabel = userCommunities
    }
    else{#Asignacion de clusters de usuarios.
        ulabel = 1 + order(runif(nUsers)) %% ngu#Aleatorios y del mismo tamano.
    }

    #Distribucion de grado ley de potencia para los usuarios.
    if(is.null(userDegrees)){
        ku = rpldis(nUsers,kmin,alpha)
        ku = pmin(ku,nObjects)#Nadie puede comprar mas que la cantidad original de objetos.
    }
    else{
        if(nUsers!=length(userDegrees)){
            warning(paste("nUsers:", nUsers, " not equal to length(userDegrees):", length(userDegrees), "\n"))
        }
        ku = userDegrees
    }

    #Tamano de comunidades.
    facShrink = runif(1, min=kShrinkMin, kShrinkMax)
    sizeClusters = aggregate(ku, list(ulabel), sum)
    sizeClusters = cbind(sizeClusters, round(facShrink * sizeClusters[,2] * pin))
    colnames(sizeClusters) = c("community", "userClusterDegree", "objectClusterSize")

    #Itero por clusters.
    uulabel = unique(ulabel)
    res = c()
    lcommunities = list()
    for(i in seq_along(uulabel)){
        iu = which(ulabel == uulabel[i])

        #Dado un cluster de usuarios me fijo que grado total tiene y elijo ktot * pin objetos de un universo de overlapObjectClusters * n.
        objsIn = sample(1:nObjects, round(facShrink * pin * sum(ku[iu])))
        objsOut = (1:nObjects)
        objsOut = objsOut[!objsOut %in% objsIn]

        prob=c(rep(pin/length(objsIn), length(objsIn)), rep((1-pin)/length(objsOut), length(objsOut)))
        for(j in seq_along(iu)){
            res=rbind(res, matrix(c(rep(iu[j],ku[iu[j]]), sample(c(objsIn,objsOut), ku[iu[j]], prob=prob, replace=FALSE)), ncol=2))
        }

        lcommunities[[i]] = list()
        lcommunities[[i]][["users"]] = iu
        lcommunities[[i]][["objects"]] = objsIn
        }
    colnames(res) = c("user", "object")

    a = res
    a[,2] = a[,2] + nUsers
    g = graph.bipartite(type=c(rep(1,nUsers), rep(0,nObjects)), edges=as.vector(as.matrix(t(a))), directed=FALSE)
    V(g)$name = c(paste("U", 1:nUsers, sep=""), paste("O", 1:nObjects, sep=""))
    V(g)$color = rgb(.8, .8, .8, .5)
    V(g)$color[1:nUsers] = rainbow(length(uulabel), alpha=0.5)[ulabel]
    V(g)$shape = "circle"
    V(g)$shape[1:nUsers] = "square"
    V(g)$label.color = rgb(0, .2, .2, .5)
    V(g)$label.cex = .48
    V(g)$frame.color = rgb(0.5, .5, .5, .5)
    V(g)$size = 8
    V(g)$size[1:nUsers] = 10
    E(g)$color = rgb(0.5, .5, .5, .6)


    return(list(graph=g, lcommunities=lcommunities, ku=ku, ulabel=ulabel))
}

#Proyecto red segun Zhou.
projectBip = function(g, type=TRUE){
    if(type){
        a = t(get.incidence(g))
    }
    else{
        a = get.incidence(g)
    }

    ku = apply(a, 1, sum)
    ko = apply(a, 2, sum)

    a = a[which(ku > 0), which(ko > 0)]

    #Normalizo por fila a_i,alpha / k_i
    a1 = a / ku[which(ku > 0)]

    #Normalizo por columna a_i,beta / k_beta
    a2 = a / ko[which(ko > 0)]

    w = t(a1) %*% a2

    gw = graph.adjacency(w, mode="directed", weighted=TRUE)
}


set.seed(12347)

a = build.bipartite.graph(NULL, NULL, nObjects, nUsers, alpha, kmin, ngu, pin, kShrinkMin, kShrinkMax)

lcommunities = a$lcommunities
g = a$graph

gg = induced.subgraph(g, unique(as.vector(get.edgelist(g))))
gg = delete.vertices(g, which(components(g)$membership != 1))

gg_inc = as_incidence_matrix(gg)
n_users = length(which(V(gg)$type))
g_users = graph_from_adjacency_matrix(sign(t(gg_inc) %*% gg_inc) - diag(rep(1, n_users), n_users, n_users), mode='undirected')
V(g_users)$color = V(gg)$color[which(V(gg)$type)]
