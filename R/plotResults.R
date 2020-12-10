#' Visualization of the obtained clusters
#' @description Plot the results of "mixClust", "gaussianClust" and "multinomialClust" objects.
#'
#' @importFrom grDevices colors
#' @importFrom graphics barplot lines par plot points
#' @importFrom stats complete.cases cov princomp prcomp
#'
#' @inheritParams summaryResults
#' @param D The dataset containing quantitative variable and/or qualitative variable.
#' @param axis 2 numeric values, to indicate the 2 axis to plot, if the dataset has more than >2 continuous variables.
#' @example mod_iris <- groupGaussianData(iris[-5], 3)
#' #plot the clustering in 2d, show the axis 2 and 4.
#' plotResults(mod_iris, iris[-5], axis = c(2,4))
#'@export
plotResults <- function(object,D, axis = c(1,2)){
    dataQuant <- D[, sapply(D, is.numeric)]
    #print(tracemem(dataQuant))
    dataQuant <- as.matrix(dataQuant)
    dataQuali <- D[, sapply(D, is.factor)]
    if(class(object)=="gaussianClust"){
        return(plot_gaussian(dataQuant, cls = object$clusters, m=object$gaussian_prms.means, s= object$gaussian_prms.mat_cov, axis))
    }else if(class(object)=="multinomialClust"){
        plot_multinomial(dataQuali, object$multinomial_parameters)
    }else if(class(object)=="mixClust"){
        plot_gaussian(dataQuant, cls = object$clusters, m=object$gaussian_prms.means, s= object$gaussian_prms.mat_cov, axis)
        plot_multinomial(dataQuali, object$multinomial_parameters)
    }
    else{
        stop("The Object must be a class of : 'mixClust', 'gaussianClust' or 'multinomialClust'")
    }
}
