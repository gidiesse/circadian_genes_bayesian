library(Rcpp)
library(RcppArmadillo)
data <- read.csv("Yreal.csv", check.names = FALSE)
sourceCpp("MCMC2112.cpp")
Sys.setenv(PKG_CXX11FLAGS = "-O3 -march=native")
#update.packages(type = "source", ask = FALSE)
data <- data.frame(lapply(data, function(x) as.numeric(gsub(",", ".", x))))
clean_data <- data[-1]
clean_data <- as.matrix(clean_data)

#write.csv(clean_data, "C:\\Users\\aless\\OneDrive - Politecnico di Milano\\foto\\Desktop\\Università\\Magistrale\\2° anno\\1° semestre\\Bayesian statistics\\Bayesian Statistics project\\Data\\clean_data.csv", row.names = FALSE,quote=FALSE)
Y<- MCMC_Circadian_Genes(clean_data)
#997.64
#471
