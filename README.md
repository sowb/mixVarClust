# mixVarClust
R package for clustering continuous and categorical data, using mixture models. 

# install mixVarClust
```r
devtools::install_github("https://github.com/sowb/mixVarClust")
```

# Examples
## Load the mixVarClust
```r
library(mixVarClust)

```

## Clustering continuous variables

```r
data("iris")
mod_gaussian <- groupGaussianData(iris[-5], 3, modelType = "diagonal", endIter = FALSE)
## not run
#summaryResults(mod_gaussian)
#plotResults(mod_gaussian, iris)  

```

## Clustering categorical variables 

```r
data("HairEyeColor")
mod_mult <- groupMultinomialData(HairEyeColor, 3, endIter = TRUE)
## not run
#summaryResults(mod_mult)
#plotResults(mod_mult, "HairEyeColor")

```

# Clustering dataset with mix features

```r
data("ToothGrowth")
# mix data, 1 categorical and 2 continuous variables
str(ToothGrowth)
mod_mix <-groupMixData(ToothGrowth, 2, modelType = "spherical")
## not run
#summaryResults(mod_mix)
#plotResults(mod_mix, ToothGrowth)

```
# Help/documentation

```r
help(package = mixVarClust)

```
