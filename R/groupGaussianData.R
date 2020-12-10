#'EM for gaussian mixtures.
#'@description This function performs a clustering of a univariate or multivariate continuous futures by using a gaussian mixture model.
#' @inheritParams groupMixData
#' @param D A matrix, or data.frame containing continuous variable(s)
#' @return "gaussianClust" object, which contains the clustering results
#' @examples data("iris")
#' mod_gaussian <- groupGaussianData(iris[-5], 3, modelType = "diagonal", endIter = FALSE)
#' ## not run
#' #summaryResults(mod_gaussian)
#' #plotResults(mod_gaussian, iris)
#'
#'@export
groupGaussianData <-
  function(D,
           K,
           modelType = "diagonal",
           randomInit = 10,
           endIter = FALSE,
           nbIter = 100) {
    #browser()
    #
    # Check for missing values, inf values
    # if (is.na(any(D)) ||
    #     is.nan(any(D)) ||
    #     is.infinite(any(D)))
    #   stop("This package does not handle NA, NaN, Inf...!")

    if (!is.matrix(D)) {
      D <- as.matrix(D)
    }
    # INITIALISATION
    init_L <- c()
    init_pk <- list()
    init_means <- list()
    init_matcov <- list()

    for (n in 1:randomInit) {
      # Latent classes
      classes <- factor(sample(seq_len(K), size = nrow(D), replace = T))
      z <- onehot(classes)[[1]]

      # Maximum likelihood estimate
      # Proportions
      pk <- matrix(colMeans(z), ncol = K)

      # means
      # faster version
      means <-
        lapply(seq_len(K), function(k)
          colMeans(as.matrix(D[which(z[, k] == 1),])))

      # Covariance matrix faster version
      mat_cov <-
        lapply(seq_len(K), function(k)
          cov(as.matrix(D[which(z[, k] == 1), ])))

      # Likelihood evaluation
      lik <- matrix(0,
                    ncol = K,
                    nrow = nrow(D),
                    byrow = T)
      for (k in seq_len(K)) {
        lik[, k] <-
          pk[, k] * dmv_gaussian(D, m = means[[k]], s = mat_cov[[k]])
      }
      lik <- sum(log(rowSums(lik)))

      init_L[n] <- lik
      init_pk[[n]] <- pk
      init_means[[n]] <- means
      init_matcov[[n]] <- mat_cov
    }

    # The minimum liklihood
    m <- which.max(init_L)

    # Get the initial parameters
    pk <- init_pk[[m]]
    means <- init_means[[m]]
    mat_cov <- init_matcov[[m]]

    # log-lik start
    L <- c(11, 1)

    # count iterations
    iter_nb <- 0
    B <- c(1)

    repeat {
      # Conditional probability that x belongs to group k
      tik <- matrix(0, ncol = K, nrow = nrow(D))
      for (k in seq_len(K)) {
        tik[, k] <-
          pk[, k] * dmv_gaussian(D, m = means[[k]], s = mat_cov[[k]])
      }
      tik <- tik / rowSums(tik)

      # New parameters----------------------------------
      # New proportions
      new_pk <- matrix(colSums(tik) / nrow(D), ncol = K)

      # New means
      new_means <-
        lapply(seq_len(K), function(k)
        {
          as.matrix(colSums(tik[, k] * D) / sum(tik[, k]))
        })

      # New covariance matrix
      Wk <- lapply(seq_len(K), function(k) {
        m <- matrix(rep(new_means[[k]], nrow(D)),
                    ncol = ncol(D),
                    byrow = TRUE)
        Wk <- as.matrix(t(D - m) %*% (as.vector(tik[, k]) * (D - m)))
        Wk
      })

      if (modelType == "general") {
        new_matcov <-lapply(seq_len(K), function(k) matCov(Wk[[k]], sum(tik[, k])))
      } else if (modelType == "diagonal") {
        new_matcov <-
          lapply(seq_len(K), function(k)
            matCov_lkbk(Wk[[k]], sum(tik[, k]), nx = ncol(D)))
      } else if (modelType == "spherical") {
        new_matcov <- lapply(seq_len(K), function(k) matCov_lkI(Wk[[k]], sum(tik[, k]), nx = ncol(D)))
      } else{
        stop("Choose a correct model type. 3 possible : 'general', 'diagonal' or 'spherical'. \n")
        break
      }

      # Update : Parameters reinitialization
      pk <- new_pk
      means <- new_means
      mat_cov <- new_matcov

      # Likelihood evaluation-----------------------------------
      lik <- matrix(0, nrow = nrow(tik), ncol = K)
      for (k in seq_len(K)) {
        lik[, k] <-
          pk[, k] * dmv_gaussian(D, m = means[[k]], s = mat_cov[[k]])
      }
      lik <- sum(log(rowSums(lik)))

      # Store the log liklihood
      L <- tail(L, 2)
      L <- c(L, lik)
      #print(lik)

      # Stopping the algorithm with predifined number of iterations
      iter_nb = iter_nb + 1
      if ((endIter == TRUE) & (iter_nb == nbIter)) {
        # cat(
        #   "Algorithm stopped using the pre-defined number of iterations. Number of iterations = ",
        #   nbIter,
        #   ".\n"
        # )
        break
      } else if # The algorithm is stopped using a threshold for the relative
      # change of the Likelihood : when the difference between the two
      # last lik is smaller than 1e-6.
      (abs(L[length(L)] - L[length(L) - 1]) <= 1e-10) {
       # cat("Algorithm stopped using the convergence of the log-likelihood. Convergence after : ", iter_nb, "iterations.\n")
        break
      }
      # b <- criterion.BIC(lik, nparams, nrow(D))
      # B <- c(B, b)

    }# end repeat

    #Nb params
    nparams = K * nbprms_guassian(ncol(D), modelType) + K - 1

    # BIC
    b <- criterionBIC(lik, nparams, nrow(D))
    # B <- c(B, b)
    #print(B)
    a <- criterionAIC(lik, nparams)

    # Get clusters
    groups <- max.col(tik)

    # Print the results
    res <- list(
      nb_iter = iter_nb,
      nb_obs = nrow(D),
      nb_var = ncol(D),
      model = modelType,
      proportions = pk,
      gaussian_prms.means = means,
      gaussian_prms.mat_cov = mat_cov,
      BIC = b,
      AIC = a,
      log.likelihood = lik,
      clusters = groups
    )
    class(res) <- "gaussianClust"
    res
  }
