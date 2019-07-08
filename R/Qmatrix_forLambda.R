
QH_ <- function(data.x, M, R, n, category ){

                Xs <- as.matrix(cbind(category, data.x))

                unique.category <- as.matrix(unlist(unique(category)))

                J <- length(unique.category)





# variables needed in this function (data.x, M, R, n )




            #---------------------------------------------
                    list.QH <- list()
                    for (j1 in 1:J){
                                  list.QHj <- list()
                                  selected_category_j1 <- unique.category[j1]
                                  selected_xs_j1 <- Xs[Xs[,1] == selected_category_j1, -1]
                                  Z1 <- Zj.var(M, selected_xs_j1, with_intercept = T, Nx = (ncol(Xs) - 1))

                                  for (j2 in 1:J){
                                                  selected_category_j2 <- unique.category[j2]
                                                  selected_xs_j2 <- Xs[Xs[,1] == selected_category_j2, -1]
                                                  Z2 <- Zj.var(M, selected_xs_j2, with_intercept = T, Nx = (ncol(Xs) - 1))
                                                  this.block <-  0 - Z1%*%tcrossprod(R,Z2)
                                                  if (j1 == j2){
                                                                this.block <- Diagonal(M * n[j1]) + this.block
                                                                }
                                                  list.QHj[[j2]] <- as.matrix(this.block %*% t(as.matrix(H.j(M,n[j2]))))
                                                }
                                  QHj <- do.call(cbind,list.QHj)
                                  list.QH[[j1]] <- QHj

                    }
                    QH <- do.call(rbind, list.QH) #row-blocks of the Z matrix
                    return(QH)
}
