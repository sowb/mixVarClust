#' EM for mixed futures, Gaussian-Mutinomial mixture models.
#'@description Clustering of a dataset with mix variables (gaussian and multinomial) and unknown labels. The quantitative features are assumed to be mixtures of gaussian distributions. The categorical features are assumed to be
#'mixtures of multinomial distributions. The quantitative and qualitative variables are assumed to be
#'independent with respect to the clusters.
#'
#'@param D A Dataframe. The dataset to perform the clustering on.
#'@param K The desired number of clusters
#'@param modelType a character, defining a family of gaussian models. 3 models
#'with free proportions are implemented : "general", "diagonal" and "spherical".
#'@param randomInit a numeric value. The number of random initializations of the EM algorithm.
#'The default is randomInit = 10.
#'@param endIter a boolean, specifying the stopping rule for the EM algorithm.
#'If endIter = TRUE, the EM alogrithm is stopped using the defined number of
#'iterations. If FALSE (default), the EM is stopped, if the incomplete likelihood converge.
#' @param nbIter a numeric value. The number of iterations desired to be performed by the EM algorithm.
#' This parameter is used only when endIter = TRUE.
#' @examples data("ToothGrowth")
#' # mix data, 1 categorical and 2 continuous variables
#' str(ToothGrowth)
#' mod_mix <-groupMixData(ToothGrowth, 2, modelType = "spherical")
#' ## not run
#' #summaryResults(mod_mix)
#' #plotResults(mod_mix, ToothGrowth)
#' @export
groupMixData <-
    function(D,
             K,
             modelType = "diagonal",
             randomInit = 10,
             endIter = FALSE,
             nbIter = 100) {
        #browser()
        #Check if the dataset is a dataframe contain
        if (!is.data.frame(D))
            stop("The dataset must a be a data.frame!")

        # Select Continuous data
        num <- sapply(D, is.numeric)
        if (!any(num))
            stop("The dataset to cluster must have columns of 'numeric' type.\n")
        dataQuant <- D[, num]
        dataQuant <- as.matrix(dataQuant)

        # # Check for missing values, inf values
        # if (is.na(any(dataQuant)) ||
        #     is.nan(any(dataQuant)) ||
        #     is.infinite(any(dataQuant)))
        #     stop("This package does not handle NA, NaN, Inf...!")

        # Select Categorical
        fac <- sapply(D, is.factor)
        if (!any(fac))
            stop("No categorical Data! The dataset must have at least one 1 column factor.\n")
        dataQuali <- D[, fac]
        if (is.factor(dataQuali))
            dataQuali <- as.data.frame(dataQuali)

        # One encoding of categorical data
        X <- onehot(dataQuali)

        # STEP 1 : INITIALISATION
        init_L <- c()
        init_pk <- list()
        init_means <- list()
        init_matcov <- list()
        init_akjh <- list()

        for (n in 1:randomInit) {
            # Latent classes
            classes <-
                factor(sample(1:K, size = nrow(D), replace = T))
            z <- onehot(classes)[[1]]

            # The initial parameters
            # Proportions
            pk <- matrix(colMeans(z), ncol = K)

            # Gaussian parameters-----------------
            # means
            # means <- list()
            # for (k in 1:K) {
            #     means[[k]] <- as.matrix(colMeans(dataQuant[which(z[, k] == 1), ]))
            # }
            # means

            # faster version
            means <-
                lapply(1:K, function(k)
                    colMeans(as.matrix(dataQuant[which(z[, k] == 1),])))

            # Covariance matrix
            # mat_cov <- list()
            # for (k in 1:K) {
            #     mat_cov[[k]] <- as.matrix(cov(dataQuant[which(z[, k] == 1), ]), ncol=ncol(D))
            # }
            # Covariance matrix faster version
            mat_cov <-
                lapply(1:K, function(k)
                    cov(as.matrix(dataQuant[which(z[, k] == 1), ])))

            # Alpha
            alpha_kjh <-
                lapply(as.data.frame(z), compute_aik, dataQuali, X)

            # The incomplete log-Likelihood evaluation
            loglik <- matrix(0, ncol = K, nrow = nrow(dataQuant))
            for (k in 1:K) {
                loglik[, k] <-
                    pk[, k] * dmv_multinomial(alpha_kjh[[k]], dataQuali, X) * dmv_gaussian(dataQuant, m = means[[k]], s = mat_cov[[k]])
            }
            loglik <- sum(log(rowSums(loglik)))

            # Store the values
            init_L[n] <- loglik
            init_pk[[n]] <- pk
            init_means[[n]] <- means
            init_matcov[[n]] <- mat_cov
            init_akjh[[n]] <- alpha_kjh
        }

        # The minimum liklihood
        m <- which.max(init_L)

        # Get the initial parameters
        pk <- init_pk[[m]]
        means <- init_means[[m]]
        mat_cov <- init_matcov[[m]]
        alpha_kjh <- init_akjh[[m]]


        # log-loglik start
        L <- c(11, 1)

        # iterations counter
        iter_nb <- 0

        repeat {
            # The conditionnal Probability
            # MAP
            tik <- matrix(0, nrow = nrow(dataQuali), ncol = K)
            for (k in 1:K) {
                tik[, k] <-
                    pk[, k] * dmv_multinomial(alpha_kjh[[k]], dataQuali, X) * dmv_gaussian(dataQuant, m = means[[k]], s = mat_cov[[k]])
            }
            tik <- tik / rowSums(tik)

            # New proportions
            new_pk <-
                matrix(colSums(tik) / nrow(dataQuali), ncol = K)

            # New means
            new_means <-
                lapply(1:K, function(k)
                {
                    as.matrix(colSums(tik[, k] * dataQuant) / sum(tik[, k]))
                })

            # # New covariance matrix
            # m <- matrix(rep(new_means[[k]], nrow(dataQuant)),
            #             ncol = ncol(dataQuant),
            #             byrow = TRUE)
            # new_matcov <-
            #     lapply(1:K, function(k) {
            #         as.matrix(t(dataQuant - m) %*% (as.vector(tik[, k]) * (dataQuant - m))) /
            #             sum(tik[, k])
            #     })
            # New covariance matrix
            Wk <- lapply(seq_len(K), function(k) {
                m <- matrix(rep(new_means[[k]], nrow(dataQuant)),
                            ncol = ncol(dataQuant),
                            byrow = TRUE)
                Wk <-
                    as.matrix(t(dataQuant - m) %*% (as.vector(tik[, k]) * (dataQuant - m)))
                Wk
            })

            if (modelType == "general") {
                new_matcov <-
                    lapply(seq_len(K), function(k)
                        matCov(Wk[[k]], sum(tik[, k])))
            } else if (modelType == "diagonal") {
                new_matcov <-
                    lapply(seq_len(K), function(k)
                        matCov_lkbk(Wk[[k]], sum(tik[, k]), nx = ncol(dataQuant)))
            } else if (modelType == "spherical") {
                new_matcov <-
                    lapply(seq_len(K), function(k)
                        matCov_lkI(Wk[[k]], sum(tik[, k]), nx = ncol(dataQuant)))
            } else{
                stop(
                    "Choose a correct model type. 3 possible models : general, diagonal or spherical. \n"
                )
                break
            }

            # Multinomial parameters
            new_alpha_kjh <-
                lapply(as.data.frame(tik), compute_aik, dataQuali, X)

            # Parameters reinitialization
            pk <- new_pk
            means <- new_means
            mat_cov <- new_matcov
            alpha_kjh <- new_alpha_kjh

            # The incomplete log-Likelihood evaluation
            loglik <- matrix(0, ncol = K, nrow = nrow(dataQuant))
            for (k in 1:K) {
                loglik[, k] <-
                    pk[, k] * dmv_multinomial(alpha_kjh[[k]], dataQuali, X) *
                    dmv_gaussian(dataQuant, m = means[[k]], s = mat_cov[[k]])
            }
            loglik <- sum(log(rowSums(loglik)))

            # Store the log liklihood
            L <- tail(L, 2)
            L <- c(L, loglik)
            #print(loglik)

            # Stopping the algorithm with predifined number of iterations
            iter_nb = iter_nb + 1
            if ((endIter == TRUE) & (iter_nb == nbIter)) {
                # cat(
                #     "Algorithm stopped using the pre-defined number of iterations. Number of iterations = ",
                #     nbIter,
                #     ".\n"
                # )
                break
            } else if # The algorithm is stopped using a threshold for the relative
            # change of the Likelihood : when the difference between the two
            # last loglik is smaller than 1e-6.
            (abs(L[length(L)] - L[length(L) - 1]) <= 1e-6) {
                # cat(
                #     "Algorithm stopped using the convergence of the log-likelihood. Convergence after : ",
                #     iter_nb,
                #     "iterations.\n"
                # )
                break
            }
            #
        }# end repeat

        #Nb params
        nparams = K * nbprms_guassian(ncol(D), modelType) +  nbprms.multinomial(D) + K -
            1

        # BIC
        b <- criterionBIC(loglik, nparams, nrow(D))
        # B <- c(B, b)
        #print(B)
        a <- criterionAIC(loglik, nparams)
        # Get clusters
        groups <- max.col(tik)

        # Print the results
        res <- list(
            nb_iter = iter_nb,
            nb_obs = nrow(D),
            nb_var = ncol(dataQuali) + ncol(dataQuant),
            clusters = groups,
            model = modelType,
            proportions = pk,
            gaussian_prms.means = means,
            gaussian_prms.mat_cov = mat_cov,
            multinomial_parameters = alpha_kjh,
            BIC = b,
            AIC = a,
            log.likelihood = loglik
        )
        class(res) <- "mixClust"
        res

    } # end bracket function
