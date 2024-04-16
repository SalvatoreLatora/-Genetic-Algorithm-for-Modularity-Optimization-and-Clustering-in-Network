dati <-read.table("dati.txt", h = T)

links_1 <- as.matrix(dati)
links <- as.vector(links_1)
vertices <- sort(unique(links))

n <- length(vertices)
N <- 400L 
L <- dim(dati)[1]



# getPopolazione allows to obtain the population of random solutions and associated information on grades. 
# - n: is the number of the nodes
# - N: is the size population
# - link: is the list of links
# - seed: is the random seed

getPopolazione <- function(n, N, link, seed = 0L){
  if (missing(n)) stop("Provide the number of nodes")
  if (missing(N)) stop("Provide the population size")
  if (missing(link)) stop("Provide the links of the network")
  set.seed(seed)
  initial_pop<- array(NA, dim = c(n, 2L, N))
  degreeK <- as.vector(table(link)) # gradi dei nodi della rete
  for(i in seq_len(N)){
    initial_pop[,1L,i] <- sample(c(1L, 2L), n,replace = TRUE)
    initial_pop[, 2L,i] <- degreeK
  }
  if(N == 1) initial_pop <- matrix(initial_pop, nrow = n, ncol = 2L)
  initial_pop
}
solDeg <- getPopolazione(n, N, links)



# modularity allows to calculate the modularity for each solution of the population
# - data: network
# - N: size population
# - initial_population: the population of solutions whose modularity is to be calculated
# - child: allows to calculate the modularity for a single solution


modularity <- function(data , N = NULL, L, initial_population = NULL, child = NULL){
  if (missing(data)) stop("Provide the network")
  if(is.null(initial_population) & is.null(child)) stop("Provide only one of the two arguments : 'initial_population' or 'child'")
  if(!(is.null(initial_population)) & !(is.null(child))) stop("Provide at least one argument: 'initial_population' or 'child'")
  links_1 <- as.matrix(dati)
  if (!(is.null(initial_population))){
    if(is.null(N)) stop("Provide the population size")
    er <- matrix(NA, nrow = N, ncol = 2L)
    ar <- matrix(NA, nrow = N, ncol = 2L)
    Q <- rep(NA, N)
    LinkPART <-array(cbind(links_1, rep(0L, L)), dim = c(dim(links_1)[1L], 3L, N))
    for(j in seq_len(N)){
      for(i in seq_len(L)){
        
        node1 = LinkPART[i, 1L, j]
        comm1 <- initial_population[node1, 1L, j] #communità del nodo 1
        
        node2 = LinkPART[i, 2L, j]
        comm2 <- initial_population[node2, 1L, j]
        
        if(comm1 == comm2)  LinkPART[i, 3L, j] <- comm1
      }
    }
    for(i in seq_len(N)){
      er[i,] <- c(length(subset(LinkPART[,,i], LinkPART[, 3L, i] == 1)[,3L])/L, length(subset(LinkPART[,,i], LinkPART[, 3L,i] == 2)[,3L])/L)
      ar[i,] <- c(((1L/(2L*L))*sum(subset(initial_population[,,i],initial_population[,1L,i]==1)[,2]))^2, ((1/(2*L))*sum(subset(initial_population[,,i],initial_population[, 1L, i]==2)[,2]))^2)
      Q[i] <- sum(er[i,] - ar[i,]) 
    }
  }
  if (!(is.null(child))){
    
    LinkPART <- cbind(links_1,rep(0, L))
    for(i in seq_len(L)){
      node1 = LinkPART[i,1L]
      comm1 <- child[node1, 1L] #communità del nodo 1
      
      node2 = LinkPART[i,2L]
      comm2 <- child[node2, 1L]
      
      if(comm1 == comm2)  LinkPART[i,3] <- comm1
      
    }
    
    er <- c(length(subset(LinkPART, LinkPART[, 3] == 1)[,3])/L, length(subset(LinkPART, LinkPART[, 3] == 2)[,3])/L) 
    ar <- c((1/(2*L))*sum(subset(child,child[,1]==1)[,2]), (1/(2*L))*sum(subset(child,child[,1]==2)[,2]))
    ar <- ar*ar
    Q <- sum(er - ar)
  }
  Q
}

Q <- modularity(data = dati ,N, L, initial_population = solDeg)



# genetic_algorithm allows to implement the genetic algorithm for the optimization of modularity
# - data: network
# - n: is the number of the nodes
# - N: size population
# - Generations: number of generations
# - links: is the list of links
# - n.links: number of links
# - initial_pop: initial population of random solutions
# - initial_modularity: modularity vectori for the initial population
# - n.point_crossover: number of point crossover (by default is equal to two)
# - verbose: allows you to observe the steps of the algorithm

genetic_algorithm <- function(data, N, Generations, n.links, links, n, initial_pop, initial_modularity, n.point_crossover = 2, verbose = FALSE){
  if(missing(data)) stop("Provide the network")
  if(missing(N)) stop("Provide the population size")
  if(missing(Generations)) stop("Provide the number of generations")
  if(missing(links)) stop("Provide the links")  # DA RIVEDERE
  if(missing(n.links)) stop("Provide the number of links")
  if(missing(initial_pop)) stop("provide the initial population")
  if(missing(initial_modularity)) stop("Provide the initial modularity vector")
  
  initial_pop_order <- array(NA, dim = c(n, 2L, N))
  bestQ <- rep(NA, Generations)
  for(i in seq_len(Generations)){
    id <- order(Q, decreasing = TRUE)
    #id variabile indicatrice che indica l'ordine delle soluzioni
    
    if (i==1){
      for(d in seq_len(N)){
        initial_pop_order[,1L,d] <- solDeg[,1L,d] 
      }
    }
    else{
      for(d in seq_len(N)){
        initial_pop_order[,,d] <- initial_pop_order[,1L,id[d], drop = FALSE] 
      }
    }
    Q <- Q[id[seq_len(N)]]
    
    n.top <- N*0.80 
    top.populations <- initial_pop_order[,1L,seq_len(n.top)]
    ProbVec <- abs(Q[seq_len(n.top)])/max(Q)
    ProbVec <- 0.5 * ProbVec
    id2 <- seq(1L:n.top)
    idprova <- cbind(id2, ProbVec)
    n.coppie <- round(n.top/2)
    
    # crossover
    
    figli <- array(NA, dim = c(n, 2L, n.top))
    indicatore <- 1 
    for(w in seq_len(n.top)){
      punto_cross <- sort(sample(1:n, size = n.point_crossover))
      punto_cross <- c(1, punto_cross,n)
      coppia <- sample(idprova[,1], size = 2L, replace = FALSE, prob = idprova[,2]) 
      indice <- rep(c(1,2), length(punto_cross))
      for(h in seq_len(n.point_crossover + 1)){
        figli[punto_cross[h]:punto_cross[h+1], 1L,w] <- top.populations[punto_cross[h]:punto_cross[h+1], coppia[indice[h]]]
      }
      indicatore <- indicatore  + 1L
    }   
    
    # mutation
    
    sol.mutation <- abind::abind(figli[,,sample(1:n.top, size = N*0.10, replace = FALSE)], initial_pop_order[,,sample(n.top:N, size = N*0.05, replace = FALSE)])
    
    m <- 2L # parametro della poisson
    numNODEStoswitch <- rpois(1, m) + 1 # numero di nodi che vogliamo cambiare
    NODEStoswitch <- sort(sample(links,size = numNODEStoswitch,replace = FALSE)) # forse TRUE
    switch_1 <- function(x) ifelse(x == 1, 2, 1) 
    for(f in seq_len(dim(sol.mutation)[3L])){
      for(j in 1:length(NODEStoswitch)){
        sol.mutation[NODEStoswitch[j],1L,f] <- switch_1(sol.mutation[NODEStoswitch[j],1L,f])
      }
    }
    
    initial_pop_order <- abind::abind(figli, sol.mutation, initial_pop_order[,,seq_len(N*0.05)])
    for (g in seq_len(N)){
      initial_pop_order[,2L,g]<- as.vector(table(links))
    }
   
    Q <- modularity(data = dati ,N = N, L = n.links, initial_population = initial_pop_order)
    
    
    bestQ[i] <- max(Q)
    BestSol <- initial_pop_order[,1L,which.max(Q)]
    Best_mod <- max(Q)
    if(verbose == TRUE) cat("GENERAZIONE:", i, "\t\t", "Best_mod=", Best_mod, "\n")
    
  }
  out <- list(Best_Modularity = bestQ,
              Modularity = Best_mod,
              Best_Solution = BestSol)
  out
}



solution <- genetic_algorithm(data = dati, N = N, n = n,Generations = 50,links = links,n.links = L, initial_pop = solDeg,initial_modularity = Q, n.point_crossover = 15, verbose = TRUE)


library(ggplot2)
step <- cbind(modularity = solution$Best_Modularity, Generation = seq(1:500))
step <- as.data.frame(step)
ggplot(step, aes(x = Generation, y = modularity))+ geom_point()+ geom_line()+ labs(title = "Genetic algorithm")



library(igraph)
rete <- graph_from_data_frame(dati, directed = FALSE)
vertex_attr(rete, "community", index= V(rete)) <- solution$Best_Solution

library(RCy3)
createNetworkFromIgraph(rete,"network")


prova <- edge.betweenness.community(graph = rete, directed = FALSE)

prova$membership
  
