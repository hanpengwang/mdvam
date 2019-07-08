sigma_estimation <- function(WX, WY, Z, X, Y, df){

                      beta.w <- ginv(crossprod(WX)) %*% crossprod(WX, WY)


                      res <- Y - X %*% beta.w
                      Wres <- WY - WX %*% beta.w
                      rss <- crossprod(res, Wres)



                      sigma.sq <- as.matrix(rss) / df
                      sigma.sq <- as.numeric(sigma.sq)


                      # lambda^j
                      beta <- solve(crossprod(Z))%*%crossprod(Z,Y)
                      e <- Y-(Z%*%beta)

                      return(list(beta, e, sigma.sq))

}
