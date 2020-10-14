# Wrap the CCPCA function in a package

require(Rcpp) 
require(RcppArmadillo) 
require(devtools) 

# This works: 
# This function create a skeleton package 
# it works with RcppArmadillo 
RcppArmadillo.package.skeleton(name = "CCregPCA", list = character(), 
     environment = .GlobalEnv, path = ".", force = FALSE, 
     code_files = character(), example_code = TRUE) 

#----> DO THIS ---->  copy the .cpp file manually in the directory 

# This is needed so the cpp functions will be exported 
# when loading the package 
Rcpp::compileAttributes(pkgdir="./CCregPCA", verbose = getOption("verbose"))

# build the package with build from devtools 
# install the package from the commandline with: 
build(pkg="./CCregPCA") 

# install.packages("./CCregPCA_1.0.tar.gz") 

require(CCregPCA) 
# READ THIS: check if all functions are there! (you may need to do a restart!)
lsf.str("package:CCregPCA")

