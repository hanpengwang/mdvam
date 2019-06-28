lambda_estimation <- function(data.x, M, category, R, sigma.sq, QH, e){



      Xs <- as.matrix(cbind(category, data.x))
      n <- as.matrix(table(category))
      J <- length(n)
      unique.category <- as.matrix(unlist(unique(category)))

      ### variables needed in this function (data.x, M, category, R, sigma.sq, QH)


      b <- c()
      v <- c()
      Gamma  <- numeric(0)




      for (j in 1:J){

                      bj <- c()
                      vj <- c()
                      alpha.j <- numeric(0) # collects the columns for all l
                      Gamma.j <- numeric(0)
                      j2 <- M * sum(n[1 : j])
                      j1 <- M * sum(n[1 : (j - 1)]) + 1


                      selected_category <- unique.category[j]

                      selected_xs <- Xs[Xs[,1] == selected_category, -1]

                      Zj <- Zj.var(M, selected_xs, with_intercept = T)


                      if (j == 1) {
                                    e.j <- as.matrix(e[1 : j2,])
                                    }

                      else        {
                                    e.j <- as.matrix(e[j1 : j2,])
                                    }

                      for (p in 1:M){ # m

                                      Pp <- Pjm.fn(M, j, p, n)
                                      Ppj <- Pjm.fnj(M, j, p, n)

                                      for (q in p:M){ # m

                                                          ###Pp <- Pjm.fn(M, j, p, n) moved to the top
                                                          Pq <- Pjm.fn(M, j, q, n)
                                                          # PpQL <- Pp%*%QL
                                                          #  PqQL <- Pq%*%QL

                                                          #-------------------------------------------------------------------- b_j
                                                          # We can't have Q, so calculate the block diagonal of Q that Pq and Pp extract:

                                                          ###Ppj <- Pjm.fnj(M, j, p, n) moved to the top
                                                          Pqj <- Pjm.fnj(M, j, q, n)

                                                          ###Zj <- Zj.var(M, j, data.x, S11.names) moved to the top

                                                          Qj.tilde <- Diagonal(M * n[j]) - Zj %*% tcrossprod(R, Zj) #R%*%t(Zj)
                                                          Qj.tilde <- as.matrix(Qj.tilde)
                                                          bjpq <- sigma.sq*(sum(
                                                                                diag(Pqj %*% tcrossprod(Qj.tilde, Ppj))
                                                                                )
                                                                            )

                                                          #bjpq <- sigma.sq*(sum(diag(Pq%*%QL%*%t(Pp))))
                                                          #-------------------------------------------------------------------- v_j
                                                          ejp <- Ppj %*% e.j
                                                          ejq <- Pqj %*% e.j
                                                          vjpq <- t(ejp) %*% ejq
                                                          #-------------------------------------------------------------------- alpha_j_pq
                                                          a.p <- Pp %*% QH
                                                          a.q <- Pq %*% QH
                                                          A.pq <- t(a.p) %*% a.q

                                                          alpha.jl.pq <- numeric(0)

                                                          for (d in 0:(J - 1)){
                                                                                d1 <- (d*M)+1
                                                                                d2 <- (d*M)+M
                                                                                Apq.j <- A.pq[d1:d2, d1:d2] # splits into submatrices MxM
                                                                                alpha.jl.pq[d+1] <- as.matrix(Apq.j[p, q]) # the only non-zero element
                                                                                ###rm(Apq.j) remove this rm to save time
                                                          }

                                                          bj <- rbind(bj, as.matrix(bjpq))
                                                          vj <- rbind(vj, as.matrix(vjpq))
                                                          alpha.j <- rbind(alpha.j, alpha.jl.pq) # this has l=1,..,J columns
                                                          ###rm(Pp,Pq,ejp,ejq,a.p,a.q,A.pq,alpha.jl.pq,vjpq,bjpq,ajpq)

                                                  }
                                      }
                      #--------------------------------------------------------------------

                      b <- rbind(b, bj)
                      v <- rbind(v, vj)
                      for (l in 1:J){

                                      Gamma.j <- cbind(Gamma.j, diag(alpha.j[,l]))
                      }

                      Gamma <- rbind(Gamma, Gamma.j)
                      ###rm(Gamma.j,bj,vj,d1,d2)
      }
      
      lambda.tilde <-  solve(Gamma) %*% (v - b)

      temp <- (M*(M-1)/2) + M
      dim(lambda.tilde) <- c(temp, J)
      
      ### rearrange lambda
      
      first_row <- 1 
      row_list <- c(first_row)
      temp <- 0
      for (i in c(1 : (M - 1))) {
        
        first_row <- first_row + M - temp 
        temp <- temp + 1 
        row_list <- append(row_list, first_row)
        
      }
      
      lambda.tilde <- lambda.tilde[c(row_list, 
                                               setdiff(
                                                        c(1:nrow(lambda.tilde)), row_list 
                                                        )
                                                                              ), ]
      
      
      return(lambda.tilde)

}



