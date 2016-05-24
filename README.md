# Replication material for: 'Diffusion indexes with sparse loadings'.
## Johannes Tang Kristensen.


[Link to the paper](http://dx.doi.org/10.1080/07350015.2015.1084308)

---

This repository contains the material necessary to replicate the empirical application in: 'Diffusion indexes with sparse loadings'. 

### Required packages 

The __forecastexp__ package is required for the estimation and the __macrods__ package
provides the data:

```r
library('devtools')
install_github('johannestang/forecastexp')
install_github('johannestang/macrods')
```

In addition the following packages should be installed : _pryr_ and _xtable_, both available from CRAN.   

Note: Total runtime is approximately 3 hours when using a 20 core machine (2 E5-2680 v2 CPUs @ 2.8 GHz). 

### /app

+ __forecastapp.R__ estimates the models and produces all output. 

