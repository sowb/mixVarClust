#' Summary of the clustering results
#'
#'@description Pretty print of the clustering results.
#'@param object mixclust object. 3 possible objects : "mixClust",
#'gaussianClust", "multinomialClust".
#'@return nice print.
#'@examples data("iris")
#'mod_gaussian <- groupGaussianData(iris[-5], 3, modelType = "diagonal", endIter = FALSE)
#'summaryResults(mod_gaussian)
#'
#'@importFrom utils head tail
#'@export
summaryResults <- function(object){
    cat("\n==============================================================\n")
    cat("Number of Observations = ", object$nb_obs)
    cat("\nNumber of variables = ", object$nb_var)
    cat("\nNumber of Clusters = ", ncol(object$proportions))
    cat("\n==============================================================\n")
    #cat("Model type : Mixture model gaussian and multinomial")
    cat("\nLog-liklihood = ", object$log.likelihood)
    cat("\nBIC criterion = ", object$BIC)
    cat("\nAIC criterion = ", object$AIC)
    cat("\n==============================================================\n")
    for (i in seq_len(ncol(object$proportions))){
        cat("\nCLUSTER ", i, "\n==============================================================\n")
        cat("** Proportion = ", round(object$proportions[,i], 3))
        if(class(object)=="gaussianClust" || class(object)=="mixClust"){
            cat("\n** Gaussian Parameters---------------------------------------- ", "\n -Means :\n")
            print(round(object$gaussian_prms.means[[i]], 3))
            cat("\n -VarCov matrix :\n ")
            print(round(object$gaussian_prms.mat_cov[[i]], 3))
        }
        if(class(object)=="multinomialClust" || class(object)=="mixClust"){
            ak <- round(object$multinomial_parameters[[i]],3)
            #colnames(ak) <-  colnames(dataQuali)
            ak[is.na(ak)] <- "."
            cat("\n** Multinomial Parameters------------------------------------- ", "\n -Alpha : \n")
            # for(x in 1:nrow(ak)) cat("\t\t\n")
            print(noquote(format(ak, justify = "right")))
        }
    }
    cat("\n==============================================================\n")
    cat("* Clusters =", head(object$clusters, 11), "....", tail(object$clusters, 11))
    cat("\n==============================================================\n")
}
