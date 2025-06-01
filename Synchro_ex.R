library(rgl)
library(RColorBrewer)
library(ggplot2)
library(igraph)

#for kuramoto we need formation of time series 
source("C:\\Users\\fiasc\\OneDrive\\Desktop\\complex network\\Materiale lab\\Code\\common.R")  #in questo modo chiamo la funzione di interesse per l'analisi 


#generation time series
M<-1000


MC<-1

#Synchro Kuramoto model for ER 

#generation of ER model 
set.seed(10)
N<- 200 #number of nodes
g<-erdos.renyi.game(N, p= 2*log(N)/N, directed= F)
E(g)$weight<-1  #weight of the edge in ER model 
#clusters(g) #controllare il significato di cosa restituisca 
#V(g)$name<-NULL
#layout <- layout.star(g)
plot(g, layout=layout)

#verifica come e se cambia la dinamica, se cambio la disposizione dei nodi (TOPOLOGIA DEL SISTEMA)

#KURAMOTO MODEL FOR WS, SBM, BA 
#BA
g1<-barabasi.game(N, power=1, m=1)
#controllare come costruire il preferential attach
clusters(g1) #controllare cosa significa avere clusters formati solo da 1 
E(g1)$weight<-1
#plot(g1)

#WS
g2<-watts.strogatz.game(dim=2, size= 5, nei=2, p= 0.3)
clusters(g2)
E(g2)$weight<-1
#plot(g2)

#sbm
pm<-matrix(c(0.01, 0.09, 0.09, 0.01), nrow=2, byrow=T) #creazione stoch matrix, each entry positive and sum of row equal to 1 
g3<-sbm.game(n=N, pref.matrix = pm, block.sizes = c(100, 100))
clusters(g3)
E(g3)$weight<-1
#plot(g3)

#Analysis for ER

MultiTS<-KURAMOTO_r(g, size=M, sd.meas.noise = 0, sigma=0.1)
plot.MultiTS(MultiTS)

MultiTS<-KURAMOTO_r(g, size=M, sd.meas.noise = 0, sigma=0.5)
plot.MultiTS(MultiTS)

MultiTS<-KURAMOTO_r(g, size=M, sd.meas.noise = 0, sigma=1)
plot.MultiTS(MultiTS)


res <- data.frame()
for(m in 1:10){
  cat(paste("MC #", m, "\n"))
  for(sigma in seq(0,0.3,0.005)){
    MultiTS <- KURAMOTO_r(g, size=M, sd.meas.noise=0., sigma=sigma)
    
    x <- 0
    y <- 0
    for(i in 1:length(MultiTS)){
      x <- x + cos(MultiTS[[i]][M])
      y <- y + sin(MultiTS[[i]][M])
    }
    x <- x/N
    y <- y/N
    
    r <- sqrt(x^2 + y^2)
    
    res <- rbind(res, data.frame(mc=m, sigma=sigma, r=r)) #inserisce res per righa dentro il data frame 
  }
}

#Analysis for WS

MultiTS<-KURAMOTO_r(g2, size=M, sd.meas.noise = 0, sigma=0.1)
plot.MultiTS(MultiTS)

MultiTS<-KURAMOTO_r(g2, size=M, sd.meas.noise = 0, sigma=0.5)
plot.MultiTS(MultiTS)

MultiTS<-KURAMOTO_r(g2, size=M, sd.meas.noise = 0, sigma=1)
plot.MultiTS(MultiTS)



res_2 <- data.frame()
for(m in 1:10){
  cat(paste("MC #", m, "\n"))
  for(sigma in seq(0,0.3,0.005)){
    MultiTS <- KURAMOTO_r(g2, size=M, sd.meas.noise=0., sigma=sigma)
    
    x <- 0
    y <- 0
    for(i in 1:length(MultiTS)){
      x <- x + cos(MultiTS[[i]][M])
      y <- y + sin(MultiTS[[i]][M])
    }
    x <- x/N
    y <- y/N
    
    r <- sqrt(x^2 + y^2)
    
    res_2 <- rbind(res_2, data.frame(mc=m, sigma=sigma, r=r)) #inserisce res per righa dentro il data frame 
  }
}



#Analysis for BA

MultiTS<-KURAMOTO_r(g1, size=M, sd.meas.noise = 0, sigma=0.1)
plot.MultiTS(MultiTS)

MultiTS<-KURAMOTO_r(g1, size=M, sd.meas.noise = 0, sigma=0.5)
plot.MultiTS(MultiTS)

MultiTS<-KURAMOTO_r(g1, size=M, sd.meas.noise = 0, sigma=1)
plot.MultiTS(MultiTS)



res_1 <- data.frame()
for(m in 1:10){
  cat(paste("MC #", m, "\n"))
  for(sigma in seq(0,0.3,0.005)){
    MultiTS <- KURAMOTO_r(g1, size=M, sd.meas.noise=0., sigma=sigma)
    
    x <- 0
    y <- 0
    for(i in 1:length(MultiTS)){
      x <- x + cos(MultiTS[[i]][M])
      y <- y + sin(MultiTS[[i]][M])
    }
    x <- x/N
    y <- y/N
    
    r <- sqrt(x^2 + y^2)
    
    res_1 <- rbind(res_1, data.frame(mc=m, sigma=sigma, r=r)) #inserisce res per righa dentro il data frame 
  }
}

#Analysis for SBM 
MultiTS<-KURAMOTO_r(g3, size=M, sd.meas.noise = 0, sigma=0.1)
plot.MultiTS(MultiTS)

MultiTS<-KURAMOTO_r(g3, size=M, sd.meas.noise = 0, sigma=0.5)
plot.MultiTS(MultiTS)

MultiTS<-KURAMOTO_r(g3, size=M, sd.meas.noise = 0, sigma=1)
plot.MultiTS(MultiTS)

res_3 <- data.frame()
for(m in 1:10){
  cat(paste("MC #", m, "\n"))
  for(sigma in seq(0,0.3,0.005)){
    MultiTS <- KURAMOTO_r(g3, size=M, sd.meas.noise=0., sigma=sigma)
    
    x <- 0
    y <- 0
    for(i in 1:length(MultiTS)){
      x <- x + cos(MultiTS[[i]][M])
      y <- y + sin(MultiTS[[i]][M])
    }
    x <- x/N
    y <- y/N
    
    r <- sqrt(x^2 + y^2)
    
    res_3 <- rbind(res_3, data.frame(mc=m, sigma=sigma, r=r)) #inserisce res per righa dentro il data frame 
  }
}
#evaluation for the strongest value of res
synchro<-c(max(res$r), max(res_1$r), max(res_2$r), max(res_3$r))

if(max(synchro)== synchro[0]){
  cat("ER model win the race for the synchronization and it's maximal value is=", res$r, "\n")
}


if ( max(synchro) == synchro[1]) {
  cat("BA model win the race for the synchronization and it's maximal value is=", res_1$r, "\n")
  
}

if ( max(synchro)== synchro[2]) {
  cat("WS model win the race for the synchronization and it's maximal value is=", res_2$r, "\n")
  
}
#al limite correggere con if !!!! controllare come è la struttura del if 
else {
  
cat("SBM model win the race for the synchronization and it's maximal value is=", res_3$r, "\n")
}


#zachary network VS configuration model, dubbio come gestisco le sequenze dei degree, ho bisogno di un vettore 

g4 <- graph("zachary")
deg_g4<-degree(g4) #to compare this grafix we need insert degree inside configuration model  

#configuration model 
g5<-sample_degseq(deg_g4, method="simple")


#now for comparision of these networks we can study: clustering coeff of networks, degrees, and then plotting 

#comparision 
par(mfrow = c(1, 2))
plot(g4, main= "zachary_club")
plot(g5, main= "configuration model")

#extracting values for comparision:

#cluster coefficient 
clusters(g4)
clusters(g5)

#degree 
par(mfrow=c(1,2))
hist(degree(g4), main="degree_distribution for zachary_model")
hist(degree(g5), main="degree_distribiution for configuration model")

#path length 
path.length.hist(g4) #for zachary
path.length.hist(g5) #for configuration model 













#metodo per alleggerire 
vector<- c(g,g1,g2,g3)

for(j in vector){
  
  MultiTS<-KURAMOTO_r(j, size=M, sd.meas.noise = 0, sigma=0.1)
  plot.MultiTS(MultiTS)
  
  MultiTS<-KURAMOTO_r(j, size=M, sd.meas.noise = 0, sigma=0.5)
  plot.MultiTS(MultiTS)
  
  MultiTS<-KURAMOTO_r(j, size=M, sd.meas.noise = 0, sigma=1)
  plot.MultiTS(MultiTS)
  
res<- data.frame(NULL)  
  
  res[j] <- data.frame()
  for(m in 1:10){
    cat(paste("MC #", m, "\n"))
    for(sigma in seq(0,0.3,0.005)){
      MultiTS <- KURAMOTO_r(j, size=M, sd.meas.noise=0., sigma=sigma)
      
      x <- 0
      y <- 0
      for(i in 1:length(MultiTS)){
        x <- x + cos(MultiTS[[i]][M])
        y <- y + sin(MultiTS[[i]][M])
      }
      x <- x/N
      y <- y/N
      
      r <- sqrt(x^2 + y^2)
      
      res[j]<- rbind(res_j,data.frame(mc=m, sigma=sigma, r=r)) #inserisce res per righa dentro il data frame 
    }
  }
  
  
}







