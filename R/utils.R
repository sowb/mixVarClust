#
dmv_gaussian <- function(X, m, s) {
    #browser()
    if (!is.matrix(X)) X <- as.matrix(X)
    if (!is.matrix(m)) m <- as.matrix(m)
    if(!is.matrix(s)) as.matrix(s)

    sig_inv <- solve(s)
    density_denominator <- (2 * pi) ^ (ncol(X) * 0.5) * (det(s)) ^ 0.5
    f <- apply(X, 1, function(x) exp(-0.5 * t(x - m) %*% sig_inv %*% (x - m)) )
    return(as.matrix(f)/density_denominator)
}
#General family
matCov <- function(w, nk){
    return(w/nk)
}
# Diagonal family
# lambda k :
matCov_lkbk <- function(w, nk, nx){
    d <- diag(w)
    d <-`diag<-`(matrix(0, nrow = length(d), ncol = length(d)), d)
    bk <- d/(det(d)^(1/nx))
    lk <- (det(d)^(1/nx))/nk
    return(lk*bk)
}
# Spherical family
matCov_lkI <- function(w, nk, nx){
    I <- diag(nrow = nx)
    lk <- sum(diag(w))/(nx * nk)
    return(lk*I)
}

# one hot encodings
onehot <- function(d) {
    # browser()
    if (!is.data.frame(d)) {
        d <- as.data.frame(d)
    }
    dummies <- list()
    for (col in 1:ncol(d)) {
        if (is.factor(d[, col])) {
            categories <- levels(d[, col])
            Z <- matrix(0, nrow = nrow(d), ncol = length(categories))
            colnames(Z) <- categories
            for (i in categories) {
                Z[which(d[, col] == i), i] <- 1
            }
            dummies[[col]] <- Z
        } else {
            stop("The column ",
                col,
                " is not a factor. All columns of the data must be a factor.\n")
        }
    }
    return(dummies)
}

# Compute the parameters of the multinomial distribution
compute_aik <- function(tik, d, X) {
    l <- lapply(d, levels)
    l <- unique(as.vector(unlist(l)))
    akjh <- matrix(NA, ncol = ncol(d), nrow = length(l))
    row.names(akjh) <- l
    for (j in 1:ncol(d)) {
        for (h in levels(d[, j]))
        {
            akjh[h, j] <- sum(X[[j]][, h] * tik) / sum(tik)
        }
    }
    return(akjh)
}

# compute the mass probability of the multinomial distribution
dmv_multinomial <- function(prob, d, x) {
    f <- matrix(0, ncol = 1, nrow = nrow(d))
    for (i in 1:nrow(d)) {
        p <- 1
        for (j in 1:ncol(d)) {
            for (h in levels(d[, j])) {
                p <- p * prob[h, j] ^ x[[j]][i, h]
            }
        }
        f[i, 1] <- p
    }
    return(f)
}

# Compute the number of parameter of the multinomial dist
nbprms.multinomial <- function(d) {
    mj<- lapply(d, function(x) length(levels(x))-1)
    sum(unlist(mj))
}

# Compute the number of parameter of the guaussian dist
nbprms_guassian <- function(ncolData, modelType) {
    if (modelType == "general"){
        nprms <- ncolData * ((ncolData + 1) / 2 + 1)
    }else if(modelType=="diagonal"){
        nprms <- ncolData * 2
    }else if(modelType == "spherical"){
        nprms = ncolData + 1
    }else{
        stop("You must choose a model between : 'general', 'diagonal', and spherical'.")
    }
}

plot_multinomial <- function(dataQuali, ak){
    K = length(ak)
    color = sample(unique(gsub('[0-9]+', '', colors())))
    par(mfrow = c(K, ncol(dataQuali)))
    for (k in seq_len(K)){
        for(j in seq_len(ncol(ak[[k]]))){
            #Sys.sleep(2)
            barplot(ak[[k]][complete.cases(ak[[k]][,j]), j], beside = T, col = color[nrow(ak[[k]])+j:nrow(ak[[k]])+j+j],
                    main = paste("Dist. by", colnames(dataQuali)[j]),
                    las = 2,
                    sub = paste("Cluster", k))
        }
    }
    par(mfrow = c(1,1))
}


plot_gaussian <- function(D, cls, m, s){
    # Do pca for high dimensional data
    if(ncol(D)>2){
        D <- princomp(D, scores = TRUE)$scores[,1:2]
        m <- lapply(unique(cls), function(x)colMeans(D[which(cls == x),]))
        s <- lapply(unique(cls), function(x)cov(D[which(cls == x),]))
        # For 1 dimensional data
    }else if(ncol(D)==1){
        D <- cbind(seq(1, nrow(D)), D)
        m <- lapply(unique(cls), function(x)colMeans(D[which(cls == x),]))
        s <- lapply(unique(cls), function(x)cov(D[which(cls == x),]))
    }
    plot(D, col=cls, pch=19)
    for(k in seq_len(length(unique(cls)))){
        addEllipse(m[[k]], s[[k]])
    }

}
