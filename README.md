# moveds
Fits models that account for non-responsive, Brownian motion of individuals 
during distance sampling surveys. 

## Install 
In R, the latest release can be installed using the <code>devtools</code> package with command 
```
devtools::install_github("r-glennie/moveds", build_vignettes = TRUE)
```
The package requires you have a C compiler installed on your system. Windows users may need to install R-tools for this reason. It is assumed Linux and Mac users have a compiler installed. 

## Help 
Two vignettes are included, describing how to use the package for line and 
point transects: 
```
library(moveds)
vignette("moveds_line_transects")
vignette("moveds_point_transects")
```
