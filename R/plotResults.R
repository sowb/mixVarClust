#' Visualization of the obtained clusters
#' @description Plot the results of "mixClust", "gaussianClust" and "multinomialClust" objects.
#'
#' @importFrom grDevices colors
#' @importFrom graphics barplot lines par plot points
#' @importFrom stats complete.cases cov princomp prcomp
#'
#' @inheritParams summaryResults
#' @param D The dataset containing quantitative variable and/or qualitative variable.
#'
#'@export
plotResults <- function(object,D){
    dataQuant <- D[, sapply(D, is.numeric)]
    #print(tracemem(dataQuant))
    dataQuant <- as.matrix(dataQuant)
    dataQuali <- D[, sapply(D, is.factor)]
    if(class(object)=="gaussianClust"){
        return(plot_gaussian(dataQuant, cls = object$clusters, m=object$gaussian_prms.means, s= object$gaussian_prms.mat_cov))
    }else if(class(object)=="multinomialClust"){
        return(plot_multinomial(dataQuali, object$multinomial_parameters))
    }else if(class(object)=="mixClust"){
        plot_gaussian(dataQuant, cls = object$clusters, m=object$gaussian_prms.means, s= object$gaussian_prms.mat_cov)
        plot_multinomial(dataQuali, object$multinomial_parameters)
    }
    else{
        stop("The Object must be a class of : 'mixClust', 'gaussianClust' or 'multinomialClust'")
    }
}
