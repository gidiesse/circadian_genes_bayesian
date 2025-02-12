library(Rcpp)
library(RcppArmadillo)
data <- read.csv("5k_genes.csv", check.names = FALSE)
sourceCpp("MCMC.cpp")#there we load the c++ file in R
#Pay attention to have this file on the same working directory
Sys.setenv(PKG_CXX11FLAGS = "-O3 -march=native") #this makes c++ code faster 
#You must use comma notation for decimal representation for data, 
#and must have time points on column and genes on the row 
#pay attention when you load
data <- data.frame(data)
clean_data <- data
clean_data <- as.matrix(clean_data)
t_ij <- seq(from = 0, to = 48, length.out = 13)#Example: change with your values
tg <- seq(from = 0, to = 150, length.out = 1500)#Example: change with your values
tg <- tg/1500
file_path=''#write there your path
#write in the vector the matrix and vector you are interest in for your analysis
vec=c('Thetaout_tot','Lambdaout','Etaout','theta_tilde','thresholds','B') #Here you can see an example
MCMC_Circadian_Genes(clean_data,t_ij,tg,file_path,10000,vec,1000)
