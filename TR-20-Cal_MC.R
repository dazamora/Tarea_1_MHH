#################################################################################
## Este Programa realiza la calibracion de una hidrografa por medio Monte Carlo
## Empleando la metodologia TR-20
## David Zamora Febrero 2016
#################################################################################

#----Inputs----
CNmin<-10
CNmax<-100

lambdamin<-0
lambdamax<-1

tcmin<-0
tcmax<-40

itermax<-1000000

#----Outputs----

#----Load Packages----
library(pracma)
library(hydroGOF)

#----Load data----  
data<-list.files(paste(getwd(),"/","DATOS_AGUACERO",sep = ""))

PrepT<-read.table(paste(getwd(),"DATOS_AGUACERO",data[4],sep="/"),h=F)
Qint<-read.table(paste(getwd(),"DATOS_AGUACERO",data[1],sep="/"),h=F)
Qout<-read.table(paste(getwd(),"DATOS_AGUACERO",data[2],sep="/"),h=F)

Qfast<-unlist(read.table(paste(getwd(),"DATOS_AGUACERO",data[3],sep="/"),h=F))
  
#----Efective Precipitation----
Parameters<-matrix(0,itermax,5);

colnames(Parameters)<-c("CN","Lambda","Tc","S","NSE")

# Random generation of parameters 
Parameters[,1]<-sample(seq(CNmin,CNmax,length=itermax),itermax,replace = F)
Parameters[,2]<-sample(seq(lambdamin,lambdamax,length=itermax),itermax,replace = F)
Parameters[,3]<-sample(seq(tcmin,tcmax,length=itermax),itermax,replace = F)

A<-0.48 # Area in km2
D<-1 # Tiempo one minute

# Cumulative Precipitation

Pc<-cumsum(PrepT[,3])
MPe<-1 # Matriz de precipitacion efectiva
MUH<-1 # Matrix de hidrogramas unitarios

for(i in 1:itermax){
  S<-25400/Parameters[i,1]-254
  Parameters[i,4]<-S
    
  Pec<-matrix(0,length(Pc),ncol=2)
  
  rule1<-which(Pc<Parameters[i,2]*S) # P efective is 0
  Pec[rule1,2]<-0
  
  if(length(rule1)!=length(Pc)){
    Pec[-rule1,2]<-(Pc[-rule1]-Parameters[i,2]*S)^2/(Pc[-rule1]+(1-Parameters[i,2]*S))
  }else{
    print(i)
  }
  
  out1<-Pec[2:end(Pec)[1],2]-Pec[1:end(Pec)[1]-1,2]
  
  Pe1<-Pec[1,2]
  Pe<-c(Pe1,out1)
  
  MPe<-cbind(MPe,Pe)
  

#----Unit Hydrograph----

A2<-A*(1/1.609)^2  # Convertir unidades de area km2 a mi2
tlag<-0.6*Parameters[i,3];  # Tiempo de retraso entre la precipitacion efectiva y el tiempo pico
tpico<-tlag+0.5*D  # Tiempo al pico del hidrograma
tbase<-8*tpico/3  # Tiempo de base del hidrograma unitario
tpicoh<-tpico/60 # Tiempo al pico en horas

Qpico<-484*A2/tpicoh*(0.3048)^3/25.4 # Q pico 

TimeT<-round(tbase/D)+1
UH<-approx(c(0,tpico,tbase),c(0,Qpico,0),c(0:(TimeT-1)))
UH$y[is.na(UH$y)]<-0

Qend<-convolve(UH$y,rev(Pe),type="o")
Parameters[i,5]<-NSE(Qend[1:1440],Qfast)

}







