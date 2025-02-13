# circadian_genes_bayesian
Repository for the Bayesian Statistics project on circadian genes

# How to run

To use the code in R, you need to download the MCMC.cpp file from the repository and place it in your working directory. After that, you can make it work in R by installing the Rcpp and RcppArmadillo libraries locally on your device.
```
install.packages(Rcpp)
install.packages(Rcpparmadillo)
```
Now you have to execute the following code: 
```
sourceCpp("MCMC.cpp")
Sys.setenv(PKG_CXX11FLAGS = "-O3 -march=native")
```
After running this line of code, you can now execute the MCMC code in R using the function "MCMC_Circadian_Genes". This function requires the following inputs:

1. A matrix of gene expression data with genes in rows and time points in columns.

2. The vector of time points "tij".

3. The vector "tg".

4. The file path where you want to save the output.

5. The number of iterations.

6. A vector of strings indicating which parameters you are interested in.

7. The burn-in value.

The code will save the resulting matrix and vector needed for your post-processing in the specified file path.
For an example, refer to the Wrapper R file in the GitHub repository.

---

With MCMC_C_version, you can run the code directly in C++ using synthetic data. In lines 53 to 55, you can choose the number of iterations, burn-in and thinning. Please note that in lines 20 and from 408 to 413, you need to specify the path where you want to save your matrix for post-processing. The following line is an example for Windows users.
```
"C:\Users\username\Documents\file.csv"
```
While this is an example for Mac/Linux users
```
/home/username/documents/file.csv
```
