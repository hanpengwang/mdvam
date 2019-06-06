## Stack variables
#------------------------------------------------------------------------------------------------------------------------------------------------
Y.var <- function(data.y, M = NULL, category, within_transform = FALSE) {

                  Ys <- as.matrix( cbind(category, data.y) )  ### we combine category with its Ys

                  unique.category <- as.matrix(unlist(unique(category)))

                  J <- length(unique.category)

                  y_stack <- 0


                  for (j in 1:J) {
                                   selected_category <- unique.category[j]
                                   n.j <- length(Ys[Ys[, 1] == selected_category, 1])
                                   y <- matrix(Ys[Ys[, 1] == selected_category, -1])
                                   y <- unlist(y)
                                   y <- matrix(y)          ### several steps , might be simplified

                                   if (within_transform == TRUE) {


                                                                   w.j <- Within.j(M, n.j)
                                                                   y <- w.j %*% y


                                                                }

                                   y_stack <- rbind(y_stack, y)

                  }

                  y_stack <- as.matrix(y_stack[-1,])
                  return (y_stack)
}


#------------------------------------------------------------------------------------------------------------------------------------------------


Zj.var <- function(M, data.x, with_intercept = FALSE){
                    if (with_intercept == TRUE) {
                                                  intercept <- rep(1, dim(data.x)[1])
                                                  z.j <- cbind(intercept, data.x)
                                                  z.j <- as.matrix(z.j)
                    }
                    else{
                                                  z.j <- as.matrix(data.x)

                    }
                    I <- Diagonal(M) # because we have the same X's for each module
                    Z.j <- as.matrix(I %x% z.j)

                    return(Z.j)
}


Z.var <- function(data.x, M, category, intercept = FALSE, within_transform = FALSE){

                  Xs <- as.matrix(cbind(category, data.x))

                  unique.category <- as.matrix(unlist(unique(category)))

                  J <- length(unique.category)

                  z <- 0
                  for (j in 1:J) {
                                    selected_category <- unique.category[j]
                                    selected_xs <- Xs[Xs[,1] == selected_category, -1]
                                    zj <- Zj.var(M, selected_xs, with_intercept = intercept)

                                    if (within_transform == TRUE) {

                                                                    n.j <- length(selected_xs[,1])
                                                                    w.j <- Within.j(M, n.j)
                                                                    zj <- w.j %*% zj

                                                                  }

                                    z <- rBind(z, zj)


                                    }
                  z <- as.matrix(z[-1,])
                  return(z)

}


##------------------------------------------------------------------------------------------------------
