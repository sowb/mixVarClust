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
    cat("\n==============================================================")
    cat("\nNumber of Observations = ", object$nb_obs)
    cat("\nNumber of Variables = ", object$nb_var)
    cat("\nNumber of Clusters = ", ncol(object$proportions))
    cat("\n==============================================================")
    #cat("Model type : Mixture model gaussian and multinomial")
    cat("\nLog-liklihood = ", object$log.likelihood)
    cat("\nBIC criterion = ", object$BIC)
    cat("\nAIC criterion = ", object$AIC)
    if(is.character(object$model)) cat("\nGaussian Model =", object$model)
    cat("\nEM algorithm converged after = ", object$nb_iter, "iterations.")
    cat("\n==============================================================")
    for (i in seq_len(ncol(object$proportions))){
        cat("\nCLUSTER ", i, "\n==============================================================")
        cat("\n+++PROPORTIONS = ", round(object$proportions[,i], 3))
        if(class(object)=="gaussianClust" || class(object)=="mixClust"){
            #cat("\nGaussian Parameters---------------------------------------- ")
            cat("\n+++MEANS :\n")
            print(round(t(object$gaussian_prms.means[[i]]), 3))
            cat("\n+++VARIANCE-COV MAT :\n ")
            print(round(object$gaussian_prms.mat_cov[[i]], 3))
        }
        if(class(object)=="multinomialClust" || class(object)=="mixClust"){
            ak <- round(object$multinomial_parameters[[i]],3)
            #colnames(ak) <-  colnames(dataQuali)
            ak[is.na(ak)] <- "."
            #cat("\nMultinomial Parameters------------------------------------- ")
            cat("\n+++PROBABILITIES alpha : \n")
            # for(x in 1:nrow(ak)) cat("\t\t\n")
            print(noquote(format(ak, justify = "right")))
        }
    }
    cat("\n==============================================================\n")
    cat("* Clusters =", head(object$clusters, 11), "....", tail(object$clusters, 11))
    cat("\n==============================================================\n")
}
