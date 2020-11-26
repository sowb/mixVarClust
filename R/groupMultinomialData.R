#' EM for Multinomial mixtures
#'@description This function performs clustering on a dataset composed by (d >= 1) independent categorical
#' variable(s) given the latent variable.
#'
#'@inheritParams groupMixData
#'@return "multinomialClust" object, which contains the clustering results
#'@export

groupMultinomialData <-
  function(D,
           K,
           randomInit = 100,
           endIter = FALSE,
           nbIter = 10) {
    # browser()
    if (is.table(D)){
      D <- as.data.frame(D)
      D <- as.data.frame(lapply(D, rep, D$Freq))[, which(colnames(D)!= "Freq")]
    }
    # Matrix of Dummies
    X <- onehot(D)

    # Initial containers
    init_L <- c()
    init_akjh <- list()
    init_pk <- list()
    for (n in seq_len(nbIter)) {
      classes <- factor(sample(1:K, size = nrow(D), replace = T))
      z <- onehot(classes)[[1]]

      # proportions
      pk <- matrix(colMeans(z), ncol = K)

      # Alpha
      alpha_kjh <- lapply(as.data.frame(z), compute_aik, D, X)

      # Compute the incomplete log liklihhood
      loglik <- matrix(0, ncol = K, nrow = nrow(D))
      for (k in seq_len(K)) {
        loglik[, k] <- pk[, k] * dmv_multinomial(alpha_kjh[[k]], D, X)
      }
      loglik <- sum(log(rowSums(loglik)))

      # Store the values
      init_L[n] <- loglik
      init_akjh[[n]] <- alpha_kjh
      init_pk[[n]] <- pk
    }

    # The minimum liklihood
    m <- which.max(init_L)

    # Initial parameters
    alpha_kjh <- init_akjh[[m]]
    pk <- init_pk[[m]]

    ######## EM algo ############
    # log-loglik start
    L <- c(11, 1)

    # count iterations
    iter_nb <- 0

    repeat {
      #----------E - Step ---------
      # MAP
      tik <- matrix(0, nrow = nrow(D), ncol = K)
      for (k in 1:K) {
        tik[, k] <- pk[, k] * dmv_multinomial(alpha_kjh[[k]], D, X)
      }
      tik <- tik / rowSums(tik)

      #----------M - Step ---------

      # New parameters
      # alphaik
      new_alpha_kjh <- lapply(as.data.frame(tik), compute_aik, D, X)

      # New proportions
      new_pik <- matrix(colSums(tik) / nrow(D), ncol = K)

      # Replace old parameters by new ones
      alpha_kjh <- new_alpha_kjh
      pk <- new_pik

      # Evaluate log liklihood
      loglik <- matrix(0, ncol = K, nrow = nrow(D))
      for (k in seq_len(K)) {
        loglik[, k] <- pk[, k] * dmv_multinomial(alpha_kjh[[k]], D, X)
      }
      loglik <- sum(log(rowSums(loglik)))

      # Store the log liklihood
      L <- tail(L, 2)
      L <- c(L, loglik)

      # Stopping the algorithm with predifined number of iterations
      iter_nb = iter_nb + 1
      if ((endIter == TRUE) & (iter_nb == nbIter)) {
        cat(
          "Algorithm stopped using the pre-defined number of iterations. Number of iterations = ",
          nbIter,
          "."
        )
        break
      } else if # The algorithm is stopped using a threshold for the relative
      # change of the Likelihood : when the difference between the two
      # last loglik is smaller than 1e-6.
      (abs(L[length(L)] - L[length(L) - 1]) <= 1e-6) {
        cat(
          "Algorithm stopped using the convergence of the log-likelihood. Convergence after : ",
          iter_nb,
          "iterations."
        )
        break
      }
      # end repeat
    }

    # Compute bic, aic, and ICL
    # Compute the total number of parameters
    nbparams <- nbprms.multinomial(D) + K - 1

    # BIC
    b <- criterionBIC(loglik, nbparams, nrow(D))
    # B <- c(B, b)

    a <- criterionAIC(loglik, nbparams)
    # A <- c(A, a)
    # Get the clusters
    group <- max.col(tik)

    # Print the results
    res <- list(
      nb_obs = nrow(D),
      nb_var = ncol(D),
      clusters = group,
      proportions = pk,
      multinomial_parameters = alpha_kjh,
      BIC = b,
      AIC = a,
      log.likelihood = loglik
    )
    class(res) <- "multinomialClust"
    res
  }
