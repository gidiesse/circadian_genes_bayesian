# circadian_genes_bayesian
Repository for the Bayesian Statistics project on circadian genes

## Building the project
The dependencies are: 
1. Armadillo - for macOS install via [homebrew](https://formulae.brew.sh/formula/armadillo)

### Build pipeline (behind the scenes)
Then use the CMakeLists.txt to produce the correspondent Makefile, and build the project with make:
```
mkdir build
cd build 
cmake ..
make all
```
The resulting executable (called `circadian genes`) is launched with:
```
./circadian genes
```
---
### Easy condensed pipeline used to build and launch
The above instructions are condensed into the the bash script `launch_program.sh`.  
Every time the main is modified, you have to rebuild the program from scratch: to do it quickly, just launch the bash script with:
```
bash launch_program.sh
```
This script will take care of building the program and launching the executable! 

To use the code in R, you have to download from the repository the MCMC.cpp code and positioned it into your working directory folder. After you can make it works in R  installing the libraries Rcpp and Rcpparmadillo locally on your device
```
install.packages(Rcpp)
install.packages(Rcpparmadillo)
```
and then you have to have the functions of the cpp on your device you have to run the following code 
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

With MCMC_C_version, you can run the code directly in C++ using synthetic data. Please note that in lines 437 to 442, you need to specify the path where you want to save your matrix for post-processing. The following line is an example for Windows users.
```
"C:\Users\NomeUtente\Documenti\file.csv"
```
While this is an example for Mac/Linux users
```
/home/nomeutente/documenti/file.csv
```
