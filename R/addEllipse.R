#'Ellipse
#'@description Draw an ellipse on top of a plot.
#'@param means a numeric vector of length 2. The center of the ellipse.
#'@param varCov a numeric value (the value will be transformed in 2x2 matrix)
#'or 2x2 matrix. The direction of the ellipse.
#'@return plot of the ellipses and their centers.
#'@export
addEllipse <- function(means, varCov){
    #browser()
    if(ncol(varCov)==1){
        varCov <- matrix(c(varCov, 0,0,varCov), ncol=2)
    }
    eigVal <- eigen(varCov)$values
    eigVec <- eigen(varCov)$vectors
    angles <- seq(0,2*pi, length.out = 300)
    ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles))
    ellRot <- eigVec %*% t(ellBase)
    X <- ellRot + as.vector(means)
    lines(t(X), asp = 2, type = "l", lwd = 4, col="purple")
    points(means[1], means[2], pch = 4, lwd = 3, col="purple")
}
