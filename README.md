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