
#RESILIENCE AND ROBUSTNESS

#punto 1)
#scrivere il codice per rimuovere lista di nodi da un network

#percolation 

library(igraph)
library(RColorBrewer)
library(ggplot2)
library(rgl)


#creation network, change between directed and undirected network
set.seed(10)
N<-100
g<-erdos.renyi.game(n=N , p= 2*log(N)/N, directed= T )

#information about the nodes 
nodes<-V(g)

num_nodes<- 40
#to remove a list of nodes, specify the number of nodes that I want to remove and then select them in the network, using sample we do a casual selecting 
random_deleting<-sample(V(g), num_nodes)

#new graph from the deleting
new_graph<-delete_vertices(g, random_deleting)

New_nodes<-V(new_graph)


#par(mfrow=c(1,2))
#plot(g)
#plot(new_graph)

#nodes
#New_nodes

#remove list of nodes from the rank ones 

#degree

degree<-degree(g, mode="all") 

#extraction of the degree with high value 
degree_ordened<-degree[order(degree, decreasing = T)] 

new_degree<-delete_vertices(g, degree_ordened)

#graph extraction and visualization 
#par(mfrow=c(1,2))
#plot(new_degree)
#plot(g)

#point 2)
N_1<-5000
set.seed(20) #give us different state of erg and barabasi

ERG<-erdos.renyi.game(n=N_1, p=2*log(N_1)/N_1, directed = T)
#need the mean degree

degree_mean<-mean(degree(ERG)) 

#imposing same degree for ERG and BARABASI 
#degree <k>= 2m
Barabasi<-barabasi.game(n=N_1,power = 1, m= degree_mean/2 )


#30% of nodes to be removed
fraction_nodes<-0.3

#fraction_nodes_high<-0.8
#approximation of num_nodes to be removed, vcount take into account the number!!!


remove_random_ERG<-function(fraction_nodes){
  
  #removing by the fraction, between 0, 1 
  num_nodes_ERG<-round(vcount(ERG)*fraction_nodes)
  nodes_removed_ERG<-sample(V(ERG), num_nodes_ERG) #from ERG
  new_graph_ERG<-delete_vertices(ERG, nodes_removed_ERG)
  LCC_erg<-components(new_graph_ERG)
  largest_ERG<-LCC_erg$csize[which.max(LCC_erg$csize)]
  
  return (largest_ERG)
}


remove_random_ERG(0.2)






remove_random_BA<-function(fraction_nodes){
  
  num_nodes_BA<-round(vcount(Barabasi)*fraction_nodes)
  
  #randomly selection in our network
  
  nodes_removed_barabasi<-sample(vcount(Barabasi), num_nodes_BA)
  
  
  new_graph_BA<-delete_vertices(Barabasi, nodes_removed_barabasi)
  
  
  LCC_ba<-components(new_graph_BA)
  largest_BA<-LCC_ba$csize[which.max(LCC_ba$csize)]
  
  
  return(largest_BA)
}

remove_random_BA(0.2)


#how plot modifies 
#par(mfrow=c(2,2))
#plot(new_graph_ERG_h)
#plot(ERG)
#plot(new_graph_BA_h)
#plot(Barabasi)

#RANDOM
#random failures BA and ERG (1000 nodes)

num_nodes_1<-1000
random_BA<-sample(V(Barabasi), num_nodes_1)

random_ERG<-sample(V(ERG), num_nodes_1)

delete_BA<-delete_vertices(Barabasi, random_BA)
delete_ERG<-delete_vertices(ERG, random_ERG)

#evaluation LCC sizes 
LCC_random_BA<-components(delete_BA)
LCC_random_ERG<-components(delete_ERG)

largest_random_BA<-LCC_random_BA$csize[which.max(LCC_random_BA$csize)]
largest_random_ERG<-LCC_random_ERG$csize[which.max(LCC_random_ERG$csize)]
#RANK NODES


degree_BA<-degree(Barabasi, mode="all") 
degree_ERG<-degree(ERG, mode="all") 

#extraction of the degree with high value 
degree_ordened_BA<-degree_BA[order(degree_BA, decreasing = T)]
degree_ordened_ERG<-degree_ERG[order(degree_ERG, decreasing = T)]

#graph
new_degree_BA<-delete_vertices(Barabasi, degree_ordened_BA)
new_degree_ERG<-delete_vertices(ERG, degree_ordened_ERG)

#evaluation LCC sizes 
LCC_degree_BA<-components(new_degree_BA)
LCC_degree_ERG<-components(new_degree_ERG)

largest_degree_BA<-LCC_degree_BA$csize[which.max(LCC_degree_BA$csize)]
largest_degree_ERG<-LCC_degree_ERG$csize[which.max(LCC_degree_ERG$csize)]




#betweenness 
#cambia il modo di fare la selezione, c'è bisogno di una soglia 


Bet_BA<-betweenness(Barabasi, v=V(Barabasi), directed = T)
Bet_ERG<-betweenness(ERG, v=V(ERG), directed= T)

#selection of a threshold
thr_BA<-quantile(Bet_BA,0.9) #per suddivisione in parti uguali 
thr_ERG<-quantile(Bet_ERG,0.9)

#selection ID nodes
selection_bet_BA<-which(Bet_BA>thr_BA)
selection_bet_ERG<-which(Bet_ERG>thr_ERG)

delete_bet_BA<-delete_vertices(Barabasi, selection_bet_BA)
delete_bet_ERG<-delete_vertices(ERG, selection_bet_ERG)

LCC_bet_BA<-components(new_degree_BA)
LCC_bet_ERG<-components(new_degree_ERG)

largest_bet_BA<-LCC_bet_BA$csize[which.max(LCC_bet_BA$csize)]
largest_bet_ERG<-LCC_bet_ERG$csize[which.max(LCC_bet_ERG$csize)]

#visualization 
par(mfrow=c(3, 2))

#random failure
plot(delete_BA, main="random_failure_Barabasi", directed= T, vertex.label=NA, vertex.size=0.1, edge.width=0.1, edge.arrow.size=0)
plot(delete_ERG, main="random_failure_ERG", directed= T, vertex.label=NA, vertex.size=0.1, edge.width=0.1, edge.arrow.size=0)

#betweenness
plot(delete_bet_BA, main="betweennes_approach_BA", directed= T, vertex.label=NA, vertex.size=0.1, edge.width=0.1, edge.arrow.size=0)
plot(delete_bet_ERG, main="betweennes_approach_ERG", directed= T, vertex.label=NA, vertex.size=0.1, edge.width=0.1, edge.arrow.size=0)

#degree
plot(new_degree_BA, main="degree_approach_BA", directed= T, vertex.label=NA, vertex.size=0.1, edge.width=0.1, edge.arrow.size=0)
plot(new_degree_ERG, main ="degree_approach_ERG", directed= T, vertex.label=NA, vertex.size=0.1, edge.width=0.1, edge.arrow.size=0)


#comparison 
largest_random_BA
largest_random_ERG

largest_bet_BA
largest_bet_ERG

largest_degree_BA
largest_degree_ERG

#controllare la rimozione se è corretta a parità di nodi da togliere 
#per migliorare il codice potrebbe essere utile creare una funzione a cui dargli in pasto numero di nodi e grafo da togliere 

#point 3)

library(igraph)

#dataexportation 
#air_traffic_metadata<-read.csv("C:\\Users\\fiasc\\Downloads\\international_air_traffic_metadata.csv")
air_traffic_edges<-read.csv("C:\\Users\\fiasc\\Downloads\\international_air_traffic.edges")
internet_edges<-read.csv("C:\\Users\\fiasc\\Downloads\\internet_AS_20000102.edges")
celegans_edges<-read.csv("C:\\Users\\fiasc\\Downloads\\celegans_ppi.edges")

#evaluation of random_failure 

#creation of a graph, miglioralo 
celegans_graph<-make_graph(celegans_edges$X1.2, directed = T)


#random_failure celegans
number_nodes_remove<-8200

#random_failure(celegans_graph)

nodes_cel<-V(celegans_graph)

LCC_initial_cel<-components(celegans_graph)
largest_size_LCC_cel<-print(LCC_initial_cel$csize[which.max(LCC_initial_cel$csize)])

random_selection_cel<-sample(nodes_cel, number_nodes_remove)

new_graph_celegans<-delete_vertices(celegans_graph, random_selection_cel)

LCC_celegans<-components(new_graph_celegans)

largest_size_LCC_celegans<-LCC_celegans$csize[which.max(LCC_celegans$csize)]
largest_size_LCC_celegans
#random failure is not so efficency, we need to extract 8200 nodes to reduce LCC 

#AIR TRAFFIC
air_traffic_graph<-make_graph(air_traffic_edges$AAA.HOI, directed= T)
nodes_air_traffic<-V(air_traffic_graph)
LCC_initial<-components(air_traffic_graph)
largest_size_LCC_tra<-print(LCC_initial$csize[which.max(LCC_initial$csize)])

random_selection_air_traffic<-sample(nodes_air_traffic, number_nodes_remove)

new_graph_traffic<-delete_vertices(air_traffic_graph, random_selection_air_traffic)

LCC_traffic<-components(new_graph_traffic)

largest_size_LCC_traffic<-LCC_traffic$csize[which.max(LCC_traffic$csize)]
largest_size_LCC_traffic
#we need more extraction here, due the high value of connection between the system 

#internet connection random failure 
internet_graph<-make_graph(internet_edges$X1.702, directed= T)
nodes_internet<-V(internet_graph)

Lcc_internet_initial<-components(internet_graph)
largest_size_LCC_internet<-print(Lcc_internet_initial$csize[which.max(Lcc_internet_initial$csize)])


random_selection_internet<-sample(nodes_internet, number_nodes_remove)

new_graph_internet<-delete_vertices(internet_graph, random_selection_internet)

LCC_internet<-components(new_graph_internet)

largest_size_LCC_internet<-LCC_internet$csize[which.max(LCC_internet$csize)]
largest_size_LCC_internet 
#comparison for RANDOM FAILURE 

#the strongest network is air_traffic's network if we study with random_failure  


#DEGREE attack CELEGANS

degree_cel<-degree(celegans_graph, mode="all")
degree_ordened_cel_indices<-order(degree_cel, decreasing = T) #indices
new_degree_cel<-delete_vertices(celegans_graph, degree_ordened_cel_indices[1:8200])
LCC_degree_cel<-components(new_degree_cel)

largest_degree_cel<-LCC_degree_cel$csize[which.max(LCC_degree_cel$csize)]
largest_degree_cel
vcount(new_degree_cel)

#even if I've extracted 8200 indices related to the high degree of nodes the network still has the same LCC's size

#DEGREE attack TRAFFIC 

degree_traffic<-degree(air_traffic_graph, mode="all")
degree_ordened_traf_indices<-order(degree_traffic, decreasing = T)
new_degree_traffic<-delete_vertices(air_traffic_graph, degree_ordened_traf_indices[1:8200])
LCC_degree_traffic<-components(new_degree_traffic)

largest_degree_traf<-LCC_degree_traffic$csize[which.max(LCC_degree_traffic$csize)]
largest_degree_traf
vcount(new_degree_traffic)

#still work!

#DEGREE attack INTERNET 


degree_internet<-degree(internet_graph, mode="all")
degree_ordened_inter_indices<-order(degree_internet, decreasing = T)
new_degree_internet<-delete_vertices(internet_graph, degree_ordened_inter_indices[1:8200])
LCC_degree_internet<-components(new_degree_internet)

largest_degree_internet<-LCC_degree_internet$csize[which.max(LCC_degree_internet$csize)]
largest_degree_internet
vcount(new_degree_internet)


#BETWEENNESS ATTACK

#celegans
Bet_cel<-betweenness(celegans_graph, v=V(celegans_graph), directed = T)
Bet_cel
#due the presence of size 2 for LCC we get that betweenness is not important here 


#traffic
Bet_traf<-betweenness(air_traffic_graph, v=V(air_traffic_graph), directed = T)
Bet_traf
#due the same structure of LCC we get same results for AIR TRAFFIC TRANSPORTATION

#INTERNET
Bet_inter<-betweenness(internet_graph, v=V(internet_graph), directed = T)
Bet_inter
#not relevant 


comparison<-data.frame(
  #extraction of the nodes with 8200 removed nodes RANDOM
  colnames= c("celegans", "traffic", "internet"),
  initial_nodes= c(vcount(celegans_graph), vcount(air_traffic_graph), vcount(internet_graph)),
  LCC_final = c(largest_size_LCC_celegans, largest_size_LCC_traffic, largest_size_LCC_internet),
  final_nodes= c(vcount(new_graph_celegans), vcount(new_graph_traffic), vcount(new_graph_internet)),
  
  #DEGREE attack
  Lcc_degree_final=c(largest_degree_cel,largest_degree_traf, largest_degree_internet),
  final_nodes_degree=c(vcount(new_degree_cel),vcount(new_degree_traffic),vcount(new_degree_internet)),
  #BETWEENNESS OR CENTRALITY DESCRIPTOR 
  betweenness=c("not relevant", "not relevant", "not relevant")
)

#EVALUATION FOR CONFIGURATION MDOEL 
#we start from erdos.renyi, remembering configuration model is model created imposing HARD CONSTRAIN fixed degree seq 
#configuration model 
set.seed(20)
N<-1000
G<-erdos.renyi.game(n=N, 2*log(N)/N, directed= T)
degree<-degree(G)
configuration_model<-sample_degseq(degree, method="simple")
#evaluation of LCC
Lcc_configuration <-components(configuration_model)
Lcc_configuration
#plot(configuration_model, directed= T, vertex.label=NA, vertex.size=0.1, edge.width=0.1, arrow.size=0.01)
vcount(configuration_model)

#ANALYSIS 
#random failure 
nodes_selected<-600
random_conf<-sample(V(configuration_model), nodes_selected)

delete_conf<-delete_vertices(configuration_model, random_conf)


#evaluation LCC sizes 
LCC_random_conf<-components(delete_conf)
LCC_random_conf

largest_random_conf<-LCC_random_conf$csize[which.max(LCC_random_conf$csize)]
largest_random_conf

vcount(delete_conf)


#DEGREE attack

degree_conf<-degree(configuration_model, mode="all")
degree_ordened_conf_indices<-order(degree_conf, decreasing = T)
new_degree_conf<-delete_vertices(configuration_model, degree_ordened_conf_indices[1:400])
LCC_degree_conf<-components(new_degree_conf)

largest_degree_conf<-LCC_degree_conf$csize[which.max(LCC_degree_conf$csize)]
largest_degree_conf
vcount(new_degree_conf)

#BETWEENNESS 

Bet_conf<-betweenness(configuration_model, v=V(configuration_model), directed = T)
#Bet_conf

thr_conf<-quantile(Bet_conf,0.6)#elimination over median 

selection_bet_conf<-which(Bet_conf>thr_conf)

delete_bet_conf<-delete_vertices(configuration_model, selection_bet_conf)

LCC_bet_conf<-components(delete_bet_conf)
largest_bet_conf<-LCC_bet_conf$csize[which.max(LCC_bet_conf$csize)]
largest_bet_conf

comparison_conf<-data.frame( 
  colnames=c("degree failure", "random failure", "betweenness failure"),
  initial_nodes=c(rep(vcount(configuration_model), 3)),
  final_nodes=c(vcount(new_degree_conf), vcount(delete_conf), vcount(delete_bet_conf)),
  LCC=c(largest_degree_conf, largest_random_conf, largest_bet_conf)
  
  )

#visualization of graphs
par(mfrow=c(1, 3))
plot(new_degree_conf, directed= T, vertex.label=NA, vertex.size=0.1, edge.width=0.1, edge.arrow.size=0, main="configuration model degree failure")
plot(configuration_model, directed= T, vertex.label=NA, vertex.size=0.1, edge.width=0.1, edgearrow.size=0, main="configuration model random failure")
plot(delete_bet_conf, directed= T, vertex.label=NA, vertex.size=0.1, edge.width=0.1, edge.arrow.size=0, main="configuration model betweenness failure")


comparison_conf
comparison


## comparison_conf
#colnames initial_nodes final_nodes LCC
#1      degree failure          1000         600 600
#2      random failure          1000         400 400
#3 betweenness failure          1000         600 600
#> comparison
#colnames initial_nodes LCC_final final_nodes Lcc_degree_final final_nodes_degree  betweenness
#1 celegans          8240         1          40                2                 40 not relevant
#2  traffic         27050         2       18850                2              18850 not relevant
#3 internet         12571         2        4371                2               4371 not relevant






