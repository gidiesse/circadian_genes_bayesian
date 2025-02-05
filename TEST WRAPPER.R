library(Rcpp)
library(RcppArmadillo)
data <- read.csv("5k_genes.csv", check.names = FALSE)
sourceCpp("MCMC2112.cpp")
Sys.setenv(PKG_CXX11FLAGS = "-O3 -march=native")
update.packages(type = "source", ask = FALSE)
data <- data.frame(lapply(data, function(x) as.numeric(gsub(",", ".", x))))
clean_data <- data[-1]
clean_data <- as.matrix(clean_data)
#data <- data.frame(data)
#clean_data <- data
#clean_data <- as.matrix(clean_data)
t_ij <- seq(from = 0, to = 48, length.out = 13)
tg <- seq(from = 0, to = 150, length.out = 1500)
tg <- tg/1500
file_path='C:\\Users\\aless\\OneDrive - Politecnico di Milano\\foto\\Desktop\\Università\\Magistrale\\2° anno\\1° semestre\\Bayesian statistics\\Bayesian Statistics project\\TEST'
#scrivi il nome delle variabili che ti interessano sotto forma di stringa
#write.csv(clean_data, "C:\\Users\\aless\\OneDrive - Politecnico di Milano\\foto\\Desktop\\Università\\Magistrale\\2° anno\\1° semestre\\Bayesian statistics\\Bayesian Statistics project\\Data\\clean_data.csv", row.names = FALSE,quote=FALSE)
vec=c('Thetaout_tot','Lambdaout','Etaout','theta_tilde','thresholds','B') 
MCMC_Circadian_Genes(clean_data,t_ij,tg,file_path,1000,vec)
