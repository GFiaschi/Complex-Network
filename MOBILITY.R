#MOBILITY

library(rgl)
library(ggplot2)
library(RColorBrewer)
library(igraph)
library(geosphere)
library(viridis)
#install.packages("minpack.lm")
library(minpack.lm)
#install.packages("leaflet")
library(leaflet)
library(htmlwidgets)





#estraction of data 
data_nodes<- read.csv("C:\\Users\\fiasc\\OneDrive\\Desktop\\complex network\\ITA_nodes.csv", header = T, sep=",")
data_edges<- read.csv("C:\\Users\\fiasc\\OneDrive\\Desktop\\complex network\\ITA_edges.csv", header= T, sep=",")


data_nodes$lat <- as.numeric(data_nodes$lat)
data_nodes$lon <- as.numeric(data_nodes$lon)
data_nodes <- data_nodes[order(data_nodes$adm_id),]
rownames(data_nodes) <- data_nodes$adm_id

# Get an idea of the orders of magnitude

sum(data_edges$estimated_flow)*1e-6

# Filter out self-loops and links with zero flow

data_edges <- data_edges[which(data_edges$from!=data_edges$to & data_edges$estimated_flow>0),]
sum(data_edges$estimated_flow)*1e-6


#str(data_edges)
#str(data_nodes)

g <- graph_from_data_frame(data_edges[,c("from", "to")], directed=T, vertices=sort(data_nodes$adm_id))
E(g)$weight <- data_edges$estimated_flow
V(g)$pop <- data_nodes$population
V(g)$label <- data_nodes$adm_name


# Set the layout for plotting the network

layout <- matrix(NA, nrow(data_nodes), 2)
layout[,1] <- data_nodes$lon
layout[,2] <- data_nodes$lat


# Visualize

sizes_nodes <- sqrt(strength(g, mode="all"))
sizes_edges <- E(g)$weight

plot(g, layout=layout, vertex.size=sizes_nodes/max(sizes_nodes), edge.width=sizes_edges/max(sizes_edges))

# fix the appearance

plot(g, layout=layout, vertex.label=V(g)$label, vertex.label.cex=.5, vertex.label.color="gray20", vertex.size=10*sizes_nodes/max(sizes_nodes), edge.arrow.size=0.05, edge.arrow.width=0.05, edge.width=5*sizes_edges/max(sizes_edges))


# Dont trust the length_km field, let us calculate the distances by ourselves
# install the geosphere package
data_edges$Haversine <- 1e-3*geosphere::distHaversine( cbind(data_nodes[data_edges$from,]$lon,data_nodes[data_edges$from,]$lat), cbind(data_nodes[data_edges$to,]$lon,data_nodes[data_edges$to,]$lat) )

# if length_km was correct, then the two measures should correlate very well, and they don't!
ggplot(data_edges, aes(length_km, Haversine)) + theme_bw() + geom_point()

num_colors<-100
full_palette <- viridis(num_colors)
num_intermediate_colors <- 20
start_index <- (num_colors - num_intermediate_colors) / 2
end_index <- start_index + num_intermediate_colors - 1
intermediate_palette <- full_palette[start_index:end_index]



#gravity model 

Gravity_power_law <- function(N_i, N_j, d_ij, a, b, gamma, k) {
  gravity <- k * (N_i^a * N_j^b) / (d_ij^gamma)
}


#extraction of values for power law fitting 
a_1<-1.363e-01
b_1<-4.157e-02
gamma<-1.539e+00
k_1<-1.275e+06

#recalling function
flux<-function(N_i, N_j, d_ij) {
  (k_1 * (N_i^a_1 * N_j^b_1)) / (d_ij^gamma)
}


#building vector for the flux in power law gravity model 

flux_pow_law<-c()
for ( i in 1:nrow(data_edges)){
  flux_pow_law[i] <- flux(data_edges$from[i], data_edges$to[i], data_edges$length_km[i])
}

#building vector for the flux exp law in gravity model

gravity_exp_law <- function(N_i, N_j, d_ij, a, b, k, d){
  #exponential law
  gravity_expo_model <- (k*(N_i^a *N_j^b))/exp(d_ij/d)
}
#extraction parameters for exponential case 
#nlsLM(estimated_flow~gravity_exp_law(data_edges$from, data_edges$to,length_km, a, b, k, d), data = data_edges, start=c(a=1.3, b=2, k=1.9, d=0.9), control=list(maxiter=500), trace=T )

a_exp<-6.631e-02
b_exp<-7.945e-02
k_exp<-6.276e+04
d_exp<-2.482e+01


#recalling function
flux_exp<-function(N_i, N_j, d_ij) {
  (k_1 * (N_i^a_exp * N_j^b_exp)) / exp(d_ij/d_exp)
}


#building vector for the flux 

flux_exp_law<-c()
for ( i in 1:nrow(data_edges)){
  flux_exp_law[i] <- flux(data_edges$from[i], data_edges$to[i], data_edges$length_km[i])
}


s_ij <- function(data_nodes, data_edges) {
  num_nodes <- length(data_nodes$population)
  s_vector <- c(0, length( nrow(data_edges)))
  
  for (i in 1:nrow(data_edges)) {
    r_i <- data_edges$length_km[i]
    
    #indexes 
    filtered_population <- data_nodes$population[
      !is.na(data_edges$length_km) & 
        !is.na(data_nodes$population) &
        data_edges$length_km < r_i & 
        data_edges$length_km != data_edges$length_km[i] & 
        data_edges$length_km != data_edges$length_km[j]
    ]
    s_vector[i] <- sum(filtered_population, na.rm = TRUE)
  }
  
  return(s_vector)
}


#RADIATION MODEL 
# results of s_ij
s_vector <- s_ij(data_nodes, data_edges)
#print(result_matrix)

#building of flux for radiation model 
radiation_vec <- c(0, length(nrow(data_edges)))
for ( i in 1:(nrow(data_edges))){
  N_i<- as.numeric(data_nodes$population[data_edges$from[i]])
  N_j<- as.numeric(data_nodes$population[data_edges$to[i]])
  S_i<-s_vector[i]
  radiation_vec[i] <- (N_i * N_j) / (N_i + S_i) * (N_i + N_j + S_i) 
}










par(mfrow=c(3,2))

#real data
w <- log10(E(g)$weight)
edge_colors <-  colorRampPalette(intermediate_palette)(length(w))
sizes_edges <- exp(-data_edges$Haversine/40)


plot(g, layout=layout, 
     vertex.label=NA, 
     vertex.label.cex=.5, 
     vertex.label.color="gray20", 
     vertex.color="tomato",
     vertex.size=10*sizes_nodes/max(sizes_nodes), 
     edge.arrow.size=0.05, 
     edge.arrow.width=0.05, 
     edge.width=5*sizes_edges/max(sizes_edges), 
     edge.color=edge_colors,
     main="Link width ~ exp(-d(i,j))")


w <- E(g)$weight
edge_colors <-  colorRampPalette(intermediate_palette)(length(w))
sizes_edges <- w

plot(g, layout=layout, 
     vertex.label=NA, 
     vertex.label.cex=.5, 
     vertex.label.color="gray20", 
     vertex.color="tomato",
     vertex.size=10*sizes_nodes/max(sizes_nodes), 
     edge.arrow.size=0.05, 
     edge.arrow.width=0.05, 
     edge.width=5*sizes_edges/max(sizes_edges), 
     edge.color=edge_colors,
     main="Link width ~ flow(i,j)")

#visualization of GRAVITY MODEL
w_1<-(E(g_1)$weight)
#edge_colors <-  colorRampPalette(intermediate_palette)(length(w_1))
sizes_edges<- w_1
edge_colors <-  colorRampPalette(intermediate_palette)(length(w_1))
plot(g_1, layout=layout, 
     vertex.label=NA, 
     vertex.label.cex=.5, 
     vertex.label.color="gray20", 
     vertex.color="tomato",
     vertex.size=10*sizes_nodes/max(sizes_nodes), 
     edge.arrow.size=0.05, 
     edge.arrow.width=0.05, 
     edge.width=5*sizes_edges/max(sizes_edges), 
     edge.color=edge_colors,
     main="Link width ~ d(i,j)^gamma")


w_2 <- log10(E(g_1)$weight)
edge_colors <-  colorRampPalette(intermediate_palette)(length(w))
sizes_edges <- exp(-flux_exp_law/40)


plot(g, layout=layout, 
     vertex.label=NA, 
     vertex.label.cex=.5, 
     vertex.label.color="gray20", 
     vertex.color="tomato",
     vertex.size=10*sizes_nodes/max(sizes_nodes), 
     edge.arrow.size=0.05, 
     edge.arrow.width=0.05, 
     edge.width=5*sizes_edges/max(sizes_edges), 
     edge.color=edge_colors,
     main="Link width ~ flux_exp_law")


#visualization of RADIATION MODEL 
w_3 <- log10(E(g_1)$weight)
edge_colors <-  colorRampPalette(intermediate_palette)(length(w))
sizes_edges <- radiation_vec


plot(g, layout=layout, 
     vertex.label=NA, 
     vertex.label.cex=.5, 
     vertex.label.color="gray20", 
     vertex.color="tomato",
     vertex.size=10*sizes_nodes/max(sizes_nodes), 
     edge.arrow.size=0.05, 
     edge.arrow.width=0.05, 
     edge.width=5*sizes_edges/max(sizes_edges), 
     edge.color=edge_colors,
     main="Link width ~ radiation model ")



#extraction of data difference 
difference_data_gravity_pow_law <- c(abs(data_edges$estimated_flow- flux_pow_law ))


difference_data_gravity_exp_law<-c(abs(data_edges$estimated_flow-flux_exp_law))

difference_gravity_exp_VS_power <- c(abs(flux_exp_law-flux_pow_law))

difference_data_rad<-c()
#controllare la dimensione degli oggetti
for (i in 1:nrow(data_edges)) {
  difference_data_rad<-c(abs(data_edges$estimated_flow[i] - radiation_vec[i]))
}

difference_data_rad_VS_gravity<-c()
for( i in 1:nrow(data_edges)){
  difference_data_rad_VS_gravity<-c(abs(flux_pow_law[i]-radiation_vec[i]), abs(flux_exp_law[i]-radiation_vec[i]))
}



#building a data frame with radiation model definition and difference data between realistic data and model data 


data_edges$gravity_flux<- flux_pow_law
#building data frame for gravity model (power law)

data_gravity_frame<-data.frame( 
  data= data_edges$gravity_flux,
  lon = data_nodes$lon[data_edges$from],
  lat= data_nodes$lat[data_edges$from]
  
  )

data_edges$radiation_flux<- radiation_vec

data_gravity_frame$data_exp <- flux_exp_law

palette <- colorNumeric(palette = "YlOrRd", domain = data_gravity_frame$data)

#ATTENZIONE INSERIRE WIDGET CON INFORMAZIONE A RIGUARDO DELLA POPOLAZIONE 


#maps for gravity model
maps_gravity <- leaflet(data_gravity_frame) %>%
  addTiles() %>%
  addCircles(lng = ~data_gravity_frame$lon, lat = ~data_gravity_frame$lat, weight = 1,
             radius = ~data_gravity_frame$data, 
             color = ~palette(data_gravity_frame$data),  
             fillOpacity = 0.7) %>%
  addLegend("bottomright", 
            pal = palette, 
            values = ~data_gravity_frame$data,
            title = "Gravity Model",
            labFormat = labelFormat(prefix = ""),
            opacity = 1)
maps_gravity

maps_gravity_exp <- leaflet(data_gravity_frame) %>%
  addTiles() %>%
  addCircles(lng = ~data_gravity_frame$lon, lat = ~data_gravity_frame$lat, weight = 1,
             radius = ~data_gravity_frame$data_exp, 
             color = ~palette(data_gravity_frame$data),  
             fillOpacity = 0.7) %>%
  addLegend("bottomright", 
            pal = palette, 
            values = ~data_gravity_frame$data_exp,
            title = "Gravity Model",
            labFormat = labelFormat(prefix = ""),
            opacity = 1)

maps_gravity_exp



data_radiation_frame<-data.frame( 
  data= data_edges$radiation_flux,
  lon = data_nodes$lon[data_edges$from],
  lat= data_nodes$lat[data_edges$from])
#rescaled due the high value of radiation flux
palette_1 <- colorNumeric(palette = "YlOrRd", domain = data_radiation_frame$data*10^(-8))

data_edges$radiation_flux
maps_radiation <- leaflet(data_radiation_frame) %>%
  addTiles() %>%
  addCircles(lng = ~data_radiation_frame$lon, 
             lat = ~data_radiation_frame$lat,
             weight = 1,
             radius = ~data_radiation_frame$data*10^(-8), 
             color = ~palette_1(data_radiation_frame$data*10^(-8)),  
             fillOpacity = 0.3) %>%
  addLegend("bottomright", 
            pal = palette_1, 
            values = ~data_radiation_frame$data*10^(-8),
            title = "Radiation Model",
            labFormat = labelFormat(prefix = ""),
            opacity = 1)
maps_radiation

palette_2<-colorNumeric(palette = "YlOrRd", domain = data_edges$estimated_flow)

real_data<-leaflet(data_edges) %>%
  addTiles() %>%
  addCircles(lng = ~data_gravity_frame$lon, lat = ~data_gravity_frame$lat, weight = 1,
             radius = ~data_edges$estimated_flow, 
             color = ~palette_2(data_edges$estimated_flow),  
             fillOpacity = 0.7) %>%
  addLegend("bottomright", 
            pal = palette_2, 
            values = ~data_edges$estimated_flow,
            title = "real data",
            labFormat = labelFormat(prefix = ""),
            opacity = 1)
 # addCircleMarkers: potrebbe fregare 
    #

real_data


       