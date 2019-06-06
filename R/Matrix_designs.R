### Sequence is consistant with technical appendix

iota.nj <- function(n.j){
                          return(
                                 matrix(
                                         rep(1,n.j),
                                         ncol=1
                                         )
                                 )
                } # a vector of dimension nj x 1 for each different category

#------------------------------------------------------------------------------------------------------------------------------------------------
H.j <- function(M, n.j, average = FALSE){

                        I.m <- Diagonal(M)
                        if (average == TRUE) {

                                    H.j.output <- I.m %x% t(iota.nj(n.j)) / n.j

                        }
                        else {

                                    H.j.output <- I.m %x% t(iota.nj(n.j))

                        }

                        return(H.j.output)

                    } # a M * Mn.j matrix

#------------------------------------------------------------------------------------------------------------------------------------------------





H <- function(M, groups, transpose = FALSE, average = FALSE){ # groups is a vector containing number of observation for each category

                      x <- list()
                      J <- length(groups)
                      for (j in 1:J){
                                    if (transpose == TRUE) {

                                           x[[j]] <- t(H.j(M, groups[j]))           # MN * MJ

                                    }
                                    else if (average == TRUE){

                                            x[[j]] <- H.j(M, groups[j], average = TRUE)   #MJ * MN

                                    }

                                    else {

                                           x[[j]] <- H.j(M, groups[j])              # MJ * MN

                                    }


                                      }
                      H.output <- bdiag(x)

                      return(H.output)

                  }

#------------------------------------------------------------------------------------------------------------------------------------------------

J.nj <- function(n.j, average = FALSE){


                      if (average == TRUE) {

                        j.nj.output <- tcrossprod(iota.nj(n.j)) / n.j

                      }
                      else {

                        j.nj.output <- tcrossprod(iota.nj(n.j))

                      }

                      return(j.nj.output)


                      }


#------------------------------------------------------------------------------------------------------------------------------------------------

Within.j <- function(M, n.j){

                    I.M.nj <- Diagonal(M * n.j)
                    HH <- Diagonal(M) %x% J.nj(n.j, average = TRUE)
                    W.j <- I.M.nj - HH
                    return(W.j)
}



#------------------------------------------------------------------------------------------------------------------------------------------------


Within <- function(M, groups){

                    wx <- list()
                    J <- length(groups)     ### groups = group with num of observation
                    for (j in 1:J){
                                    n.j <- groups[j]
                                    wx[[j]] <- Within.j(M, n.j)

                                    }
                    return(bdiag(wx))

}


#------------------------------------------------------------------------------------------------------------------------------------------------
### sub error vector for each category(schools)

Pjm.fnj <- function(M, j, m, n){

                      insert1 <- Matrix(0,n[j], ((m - 1) * n[j]))

                      insert2 <- Diagonal(n[j])

                      insert3 <- Matrix(0, n[j], ((M - m) * n[j]))

                      insert <- as.matrix(cBind(insert1, insert2, insert3))

                      return(insert)

}

#------------------------------------------------------------------------------------------------------------------------------------------------

Pjm.fn <- function(M, j, m, n){

                            p <- 0
                            for (i in 1 : length(n)){

                                            if (i == j) {
                                              insert1 <- Matrix(0, n[j], ((m - 1) * n[i]))

                                              insert2 <- Diagonal(n[j])

                                              insert3 <- Matrix(0, n[j], ((M - m) * n[i]))

                                              insert <- cBind(insert1, insert2, insert3)

                                            }
                                            else {
                                              insert <- Matrix(0, n[j], M * n[i])
                                            }

                                            p <- cBind(p, insert)
                                     }

                            p <- as.matrix(p[,-1]) # delete first column
                            return(p)

}

