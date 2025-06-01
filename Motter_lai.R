#motter-lai model

library(igraph)

N <- 1000
g <- erdos.renyi.game(n = N, p = 2 * log(N) / N, directed = TRUE)
alpha <- 0.1

# Funzione per calcolare la capacità iniziale per ciascun nodo
initial_capacity <- function(alpha, g) {
  initial_load <- betweenness(g, v = V(g)) / sum(betweenness(g, v = V(g)))
  capacity <- (1 + alpha) * initial_load
  names(capacity) <- V(g)  # Associa la capacità agli ID dei nodi
  return(capacity)
}

motter_lai_model <- function(g, alpha) {
  
  Lcc_size <- c()
  capacity <- initial_capacity(alpha, g)
  
  #cat("Initial capacity length:\n")
  #print(length(capacity))
  
  # Stress iniziale
  selection_node <- sample(V(g), 1)
  g <- delete_vertices(g, selection_node)
  
  while (vcount(g) > 1) {
    
    bet <- betweenness(g)
    load <- bet / sum(bet)
    names(load) <- V(g)  # Associa la load agli ID dei nodi
    
    #cat("Load length:\n")
    #rint(length(load))
    
    # Trova i nodi rimanenti e le loro capacità
    remaining_nodes <- V(g)
    remaining_capacity <- capacity[names(capacity) %in% remaining_nodes]
    remaining_load <- load[names(load) %in% remaining_nodes]
    
    # Trova i nodi sovraccarichi
    overload_nodes <- remaining_nodes[remaining_load > remaining_capacity]
    #cat("Overloaded nodes:\n")
    #print(overload_nodes)
    
    # Se non ci sono nodi sovraccarichi, esci dal ciclo
    if (length(overload_nodes) == 0) {
      Lcc <- components(g)
      Lcc_max <- max(Lcc$csize)
      Lcc_size <- c(Lcc_size, Lcc_max)
      break 
    }
    
    # Rimuovi i nodi sovraccarichi
    g <- delete_vertices(g, overload_nodes)
    
    # Ricalcola la capacità per il grafo aggiornato
    capacity <- initial_capacity(alpha, g)
  }
  
  return(Lcc_size)
}

# Esecuzione del modello
Lcc_size <- motter_lai_model(g)
cat("Final Lcc_size:\n")
print(Lcc_size)



set.seed(20)
N_1<-5000
ERG<-erdos.renyi.game(n=N_1, p=2*log(N_1)/N_1, directed = T)
m<-mean(degree(ERG))/2
g_1<-barabasi.game(n=N_1, m=m, directed = T)



#results of largest connected components
Lcc_size_erg<-motter_lai_model(ERG, alpha=0.1)
Lcc_size_erg
Lcc_size_BA<-motter_lai_model(g_1, alpha=0.1)
Lcc_size_BA
#print(Lcc_size_erg)
#print(Lcc_size_BA)

#evaluation of tolerance in the system 

#for different alphas 

ERG<-erdos.renyi.game(n=N_1, p=2*log(N_1)/N_1, directed = T)
m<-mean(degree(ERG))/2
g_1<-barabasi.game(n=N_1, m=m, directed = T)


Lcc_size_erg_0.9<-motter_lai_model(ERG, alpha=0.9)

Lcc_size_BA_0.9<-motter_lai_model(g_1, alpha=0.9)

print(Lcc_size_erg)
print(Lcc_size_BA)


Lcc_size_erg_0.5<-motter_lai_model(ERG, alpha=0.5) #modified value 
Lcc_size_BA_0.5<-motter_lai_model(g_1, alpha=0.5)
print(Lcc_size_erg)
print(Lcc_size_BA)

comparison<-data.frame( 
  colnames=c("Lcc_ERG_0.1", "Lcc_Ba_0.1", "Lcc_ERG_0.5", "Lcc_BA_0.5", "Lcc_ERG_0.9", "Lcc_Ba_0.9"),
  data=c(Lcc_size_erg, Lcc_size_BA, Lcc_size_erg_0.5, Lcc_size_Ba_0.5, Lcc_size_erg_0.9, Lcc_size_BA_0.5)
  )

#changing for different seeds

















#point 3)
#real network


#building network
# Leggi il file .edges con spazi come delimitatori
edges_data_us_power <- read.table("C:\\Users\\fiasc\\Downloads\\USpowergrid.edges", header = FALSE, sep = " ")

# Assegna nomi alle colonne
colnames(edges_data_us_power) <- c("Source", "Target", "Weight")


# Crea una matrice di archi usando solo le prime due colonne
edge_list <- as.matrix(edges_data_us_power[, c("Source", "Target")])


#graph
us_power_grid_graph <- graph_from_edgelist(edge_list, directed = TRUE)

#weights
edge_weights <- edges_data_us_power$Weight
E(us_power_grid_graph)$weight <- edge_weights

#nodes
vcount(us_power_grid_graph)

#summary(g)
cat("Numero di nodi:", vcount(g), "\n")
cat("Numero di archi:", ecount(g), "\n")


# Visualizza il grafo
plot(us_power_grid_graph, vertex.label=NA, vertex.size=0.1, arrow.size=0.001 )

#evaluation of cascade failures for power-grid in motter_lai system I'm forcing it





# Funzione per calcolare la capacità iniziale per ciascun nodo
initial_capacity <- function(alpha, g) {
  initial_load <- betweenness(g, v = V(g)) / sum(betweenness(g, v = V(g)))
  capacity <- (1 + alpha) * initial_load
  names(capacity) <- V(g)  # Associa la capacità agli ID dei nodi
  return(capacity)
}

motter_lai_model <- function(g, alpha) {
  
  Lcc_size <- c()
  capacity <- initial_capacity(alpha, g)
  
  #cat("Initial capacity length:\n")
  #print(length(capacity))
  
  # Stress iniziale
  selection_node <- sample(V(g), 1)
  g <- delete_vertices(g, selection_node)
  
  while (vcount(g) > 1) {
    
    bet <- betweenness(g)
    load <- bet / sum(bet)
    names(load) <- V(g)  # Associa la load agli ID dei nodi
    
    #cat("Load length:\n")
    #rint(length(load))
    
    # Trova i nodi rimanenti e le loro capacità
    remaining_nodes <- V(g)
    remaining_capacity <- capacity[names(capacity) %in% remaining_nodes]
    remaining_load <- load[names(load) %in% remaining_nodes]
    
    # Trova i nodi sovraccarichi
    overload_nodes <- remaining_nodes[remaining_load > remaining_capacity]
    #cat("Overloaded nodes:\n")
    #print(overload_nodes)
    
    # Se non ci sono nodi sovraccarichi, esci dal ciclo
    if (length(overload_nodes) == 0) {
      Lcc <- components(g)
      Lcc_max <- max(Lcc$csize)
      Lcc_size <- c(Lcc_size, Lcc_max)
      break 
    }
    
    # Rimuovi i nodi sovraccarichi
    g <- delete_vertices(g, overload_nodes)
    
    # Ricalcola la capacità per il grafo aggiornato
    capacity <- initial_capacity(alpha, g)
  }
  
  return(Lcc_size)
}





library(ggplot2)

# Funzione per calcolare la dimensione del LCC per un dato alpha
calculate_lcc_for_alpha <- function(alpha) {
  # Assicurati che la funzione ritorni una lista di dimensioni del LCC
  result <- motter_lai_model(alpha, us_power_grid_graph)
  # Controlla e ritorna la dimensione massima del LCC se esiste
  return(if (length(result) > 0) max(result) else 0)
}

# Gamma di valori per alpha
alpha_values <- seq(0.01, 1, by = 0.01)
lcc_sizes <- numeric(length(alpha_values))

# Calcola la dimensione del LCC per ciascun alpha
for (i in seq_along(alpha_values)) {
  alpha <- alpha_values[i]
  cat("Calcolando per alpha =", alpha, "\n")  # Messaggio di debug
  result <- calculate_lcc_for_alpha(alpha)
  lcc_sizes[i] <- result
}

# Crea un data frame per ggplot
plot_data <- data.frame(alpha = alpha_values, lcc_size = lcc_sizes)

# Scatter plot
ggplot(plot_data, aes(x = alpha, y = lcc_size)) +
  geom_point() +
  labs(title = "Dimensione del LCC al variare di alpha",
       x = expression(alpha),
       y = "Dimensione del LCC") +
  theme_minimal()






library(ggplot2)
Lcc_us_power_grid<-motter_lai_model(us_power_grid_graph, alpha=0.1)
Lcc_us_power_grid_2<-motter_lai_model(us_power_grid_graph, alpha=0.5)
Lcc_us_power_grid_3<-motter_lai_model(us_power_grid_graph, alpha=0.999)
motter_lai_model(us_power_grid_graph, alpha=0.000001)
motter_lai_model(us_power_grid_graph, alpha=0.3)
motter_lai_model(us_power_grid_graph, alpha=0.7)

#explore another point around zero !!
lcc_sizes<-c(Lcc_us_power_grid, Lcc_us_power_grid_2, Lcc_us_power_grid_3)
plot_data <- data.frame(alpha =c(0.1, 0.5, 0.999, 0.000001, 0.3, 0.7), 
                        lcc_size =c(Lcc_us_power_grid, Lcc_us_power_grid_2, Lcc_us_power_grid_3,motter_lai_model(us_power_grid_graph, alpha=0.000001), motter_lai_model(us_power_grid_graph, alpha=0.3), motter_lai_model(us_power_grid_graph, alpha=0.7)) )

# Scatter plot, inserire qualche punto in più 
ggplot(plot_data, aes(x = alpha, y = lcc_size)) +
  geom_point() +
  labs(title = "scatter plot_depending_alpha_us_power_grid",
       x = expression(alpha),
       y = "LCc_size") +
  theme_minimal()






#maximum value that can assume tolerance to reduce lcc
Lcc_us_power_grid_9<-motter_lai_model(us_power_grid_graph, alpha=0) #0.000001 (critico)

#discontinous phace transitions for LCC !!!













#case of configuration model

set.seed(20)
N<-4941
G<-erdos.renyi.game(n=N, 2*log(N)/N, directed= T)
degree<-degree(G)
configuration_model<-sample_degseq(degree, method="simple")

motter_lai_model(configuration_model, alpha=0.1)

motter_lai_model(configuration_model, alpha=0.5)

motter_lai_model(configuration_model, alpha=0.9)








#cosa interessante sarebbe valutare su uno scatter plot l'andamento della LCC al variare di alpha 


