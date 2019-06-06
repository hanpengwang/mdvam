# beta_estimation <- function(x, y, M, col_category){
#
#                     wx <- Z.var(x, M, col_category, intercept = F, within_transform = T)
#                     wy <- Y.var(y, M, col_category, within_transform = T)
#
#                     return(list(wx, wy))
#
#
# }






beta_estimation <- function(Z, omega.inv){



                              var.betaGLS <-  solve(t(Z) %*% (omega.inv) %*% Z)
                              beta.gls <- var.betaGLS %*% (t(Z) %*% omega.inv %*% Y)

                              v.betaGLS <- diag(var.betaGLS)

                              sd.betaGLS <- sqrt(v.betaGLS)

                              tstat <- beta.gls/sd.betaGLS

                              #lambda.tilde.gls <- lambda.tilde
                              #sigma.sq.gls <- sigma.sq

                              return(beta.gls)





}
