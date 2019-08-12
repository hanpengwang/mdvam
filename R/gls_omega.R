omegaj.fn <- function(M, n.j, sigma.sq, lambda.j){


                          J.Mnj <- tcrossprod(iota.nj(M * n.j))
                          J.nj <- tcrossprod(iota.nj(n.j))
                          omega.j <- (lambda.j) %x% J.nj + diag(1, nrow = M * n.j) %x% sigma.sq
                          return(omega.j)


}

gls_omega <- function(M, sigma.sq, lambda.tilde, groups){

                                ### variables needed in this function (M, sigma.sq, lambda.j, groups )

                                lst.omega <- list()
                                lst.lambda <- list()

                                lambda_j <- diag(M)

                                J <- length(groups)
                                
                                browser()
                                
                                for (j in 1:J){

                                              lambdaj <- lambda_j  # depends on M
                                              lambdaj[lower.tri(lambdaj, diag=TRUE)] <- lambda.tilde[,j]
                                              lambdaj[upper.tri(lambdaj, diag=FALSE)] <- 0
                                              lambdaj <- lambdaj + t(lambdaj) - diag(diag(lambdaj)) # this is our \lambda^j matrix

                                              lst.lambda[[j]] <- lambdaj

                                              n.j <- groups[j]

                                              omega.j <- omegaj.fn(M, n.j, sigma.sq, lambdaj)

                                              lst.omega[[j]] <- as.matrix(solve(omega.j))


                                }

                                omega.inv <- bdiag(lst.omega)
                                
                                

                                return(list(omega.inv, lst.omega))

}
