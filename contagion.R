

#CONTAGION 

library(igraph)
library(ggplot2)
library(deSolve)
library(reshape2)
library(viridis)

#source("C:\\Users\\fiasc\\OneDrive\\Desktop\\complex network\\Materiale lab\\Code\\common.R")

#solving ODE 
#we need function to describe a system of edo 
data_nodes<- read.csv("C:\\Users\\fiasc\\OneDrive\\Desktop\\complex network\\ITA_nodes.csv", header = T, sep=",")
data_edges<- read.csv("C:\\Users\\fiasc\\OneDrive\\Desktop\\complex network\\ITA_edges.csv", header= T, sep=",")

N<-length(data_nodes$population)

S_0<-data_nodes$population
S_0[26]<-data_nodes$population[data_nodes$adm_name=="Padova"] -20


I_0<-rep(0, N)
I_0[which(data_nodes$adm_name== "Padova")]<-20
R_0<-rep(0,N)

state<-c(S_0, I_0, R_0)
state[is.na(state)]<-0

S <- state[1:N]
I <- state[(N + 1):(2 * N)]
R <- state[(2 * N + 1):(3 * N)]


#costruzione della estimated_flow_matrix
flow_matrix<-matrix(0, nrow=nrow(data_edges), ncol=nrow(data_edges))

for ( i in 1:length(data_edges$estimated_flow)){
  for ( j in 1:length(data_edges$estimated_flow)){
    if(i !=j){
      from<-data_edges$from[i]
      to<-data_edges$to[j]
      flow_matrix[from, to]<-data_edges$estimated_flow[i]
    }
  }
  
}


id_to_index <- setNames(seq_along(unique(data_nodes$adm_id)), unique(data_nodes$adm_id))

# Ottieni gli indici delle città per la matrice di flussi
city_indices <- id_to_index[city_nodes$adm_id]

# Filtra la matrice di flussi per mantenere solo i nodi delle città
flow_matrix <-flow_matrix[city_indices, city_indices]


# Funzione per calcolare e normalizzare i flussi in uscita, conviene scrivere più funzioni diverse a restituire la stessa quantità
normalize_outflows <- function(flow_matrix, S) {
  
  num_cities <- nrow(flow_matrix)
  
  # Calcola i flussi totali in uscita per ciascuna città
  total_outflows <- rowSums(flow_matrix)
  
  # Inizializza il vettore dei flussi normalizzati
  normalized_outflows <- numeric(num_cities)
  
  
  for (i in 1:num_cities) {
    if (total_outflows[i] > S[i]) {
      # Se i flussi totali superano la popolazione disponibile, ridimensiona
      scaling_factor <- S[i] / total_outflows[i]
      normalized_outflows[i] <- S[i]
      flow_matrix[i, ] <- flow_matrix[i, ] * scaling_factor
    } else {
      normalized_outflows[i] <- total_outflows[i]
    }
  }
  
  return(flow_matrix)
}


sir_model_with_flows <- function(time, state, params, flow_matrix) {
  with(as.list(c(state, params)), {
    num_cities <- length(S)
    
    # Calcola e normalizza i flussi in uscita
    normalized_flow_matrix <- normalize_outflows(flow_matrix, S)
    
    
    
    # Calcola i flussi in uscita normalizzati
    outflows_S <- rowSums(normalized_flow_matrix*S)
    outflows_I <- rowSums(normalized_flow_matrix * I)  
    outflows_R <- rowSums(normalized_flow_matrix * R)  
    
    # Calcola i flussi in entrata
    inflows_S <- colSums(normalized_flow_matrix * S)
    inflows_I <- colSums(normalized_flow_matrix * I)
    inflows_R <- colSums(normalized_flow_matrix * R)
    
    
    # Aggiornamento delle equazioni differenziali
    dS <- -beta * S * I / (S+I+R) + (1 / epsilon) * (inflows_S - outflows_S)
    dI <- beta * S * I / (S+I+R) - gamma * I + (1 / epsilon) * (inflows_I - outflows_I)
    dR <- gamma * I + (1 / epsilon) * (inflows_R - outflows_R)
    
    
    # Restituzione dei derivati
    return(list(c(dS, dI, dR)))
  })
}
  
params<-c(beta=0.5 , gamma=0.2, epsilon=1)
params_1<-c(beta=0.5 , gamma=0.2 , epsilon=10^2)
params_2<-c(beta=0.5 , gamma=0.2, epsilon=10^4)
  
Time<-seq(from= 1, to = 365, by = 1)
out<-ode(y = state, times = time, func = sir_model_with_flows, parms = params, flow_matrix=flow_matrix)
out_1<-ode(y = state, times = time, func = sir_model_with_flows, parms = params_1, flow_matrix=flow_matrix) 
out_2<-ode(y = state, times = time, func = sir_model_with_flows, parms = params_2, flow_matrix=flow_matrix)
  
  
  
  #extraction results, qua dentro si ripetono i nodi così come si erano presentati inizialmente 
I_values<-out[,112:221]
I_values_1<-out_1[,112:221]
I_values_2<-out_2[,112:221]
  
  
heatmap_matrix<-matrix(0, nrow=N, ncol=length(Time))
heatmap_matrix_1<-matrix(0, nrow=N, ncol=length(Time))
heatmap_matrix_1<-matrix(0, nrow=N, ncol=length(Time))
  
for ( i in 1:length(data_nodes$population)){
  heatmap_matrix[i,]<-I_values[,i]/data_nodes$population[i]
  heatmap_matrix_1[i,]<-I_values_1[,i]/data_nodes$population[i]
  heatmap_matrix_2[i,]<-I_values_2[,i]/data_nodes$population[i]
  }


library(ggplot2)
library(reshape2)

# Supponiamo che heatmap_matrix sia una matrice N x 365
num_weeks <- 52
N <- nrow(heatmap_matrix)

# Creare una nuova matrice per le settimane
weekly_heatmap_matrix <- matrix(0, nrow = N, ncol = num_weeks)
weekly_heatmap_matrix_1 <- matrix(0, nrow = N, ncol = num_weeks)
weekly_heatmap_matrix_2 <- matrix(0, nrow = N, ncol = num_weeks)


# Aggregare i dati per settimane
for (week in 1:num_weeks) {
  start_day <- (week - 1) * 7 + 1
  end_day <- min(week * 7, ncol(heatmap_matrix))
  
  # Prendere la media dei valori in quella settimana
  weekly_heatmap_matrix[, week] <- rowMeans(heatmap_matrix[, start_day:end_day])
  weekly_heatmap_matrix_1[, week] <- rowMeans(heatmap_matrix_1[, start_day:end_day])
  weekly_heatmap_matrix_2[, week] <- rowMeans(heatmap_matrix_2[, start_day:end_day])
  
  
  }


# Creare i nomi delle settimane e dei nodi
weeks <- paste0("Week", 1:num_weeks)
nodes <- paste0("Node", N:1)

#matrice per epsilon 1
melt_data <- melt(weekly_heatmap_matrix, varnames = c("ID_NODES", "WEEKS"))

melt_data_1 <- melt(weekly_heatmap_matrix_1, varnames = c("ID_NODES", "WEEKS"))
melt_data_2 <- melt(weekly_heatmap_matrix_2, varnames = c("ID_NODES", "WEEKS"))

negative_indeces<-melt_data[melt_data<0]

for ( i in 1:length(melt_data$value)){
  if(melt_data$value[i]<0){
    melt_data$value[i]<- - melt_data$value[i]
    melt_data_1$value[i]<-  - melt_data_1$value[i]
    melt_data_2$value[i]<-  - melt_data_2$value[i]
  }
  else {
    melt_data$value<-melt_data$value
    melt_data_1$value<-melt_data_1$value
    melt_data_2$value<-melt_data_2$value
  }
  
}


#heatmap epsilon 1
heatmap_plot <- ggplot(melt_data, aes(x = WEEKS, y = ID_NODES, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white",high = "black") +  # Personalizza i colori
  theme_minimal() +
  labs(x = "weeks", y = "nodes", fill = "Intensity")+
  ggtitle("epsilon 1")



#heatmap epsilon 1
heatmap_plot_1 <- ggplot(melt_data_1, aes(x = WEEKS, y = ID_NODES, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +  # Personalizza i colori
  theme_minimal() +
  labs(x = "weeks", y = "nodes", fill = "Intensity")+
  ggtitle("epsilon 10^2")



#heatmap epsilon 1
heatmap_plot_2<- ggplot(melt_data_2, aes(x = WEEKS, y = ID_NODES, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +  # Personalizza i colori
  theme_minimal() +
  labs(x = "weeks", y = "nodes", fill = "Intensity")+
  ggtitle("epsilon 10^4")


library(gridExtra)

grid.arrange(heatmap_plot, heatmap_plot_1, heatmap_plot_2, ncol = 3)

#da ripetere per GRAVITY e RADIATION model 



















