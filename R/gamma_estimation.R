fnc.form.Aj.theta5M <- function(deriv.alphaj.s){ # fucntion create variance covariance matrix for lambda for a single school 
                                  vector_len <- length(deriv.alphaj.s)

                                  num_alpha <- (-1 + sqrt(1 + 4 * 2 * vector_len)) / 2

                                  d <- diag(deriv.alphaj.s[1:num_alpha])

                                  row_length <- dim(d)[1]

                                  first_column <- 2
                                  covs <- deriv.alphaj.s[(row_length + 1):length(deriv.alphaj.s)]
                                  covs_count <- 1

                                  for (row in c(2:row_length - 1)) {

                                                  for (column in c(first_column:row_length)) {

                                                                                d[row,column] <- covs[covs_count]
                                                                                d[column,row] <- covs[covs_count]
                                                                                covs_count <- covs_count + 1


                                                                      }
                                                  first_column <- first_column + 1

                                                }

                                  return(d)

                  }


### variables needeed in this function (lambda.tilde.mle, beta.mle, l.Omega.inv,  )


VA_estimation <- function(lambda.tilde.mle, beta.mle, l.Omega.inv, J, Y, Z, n, M){
                                
                               
                                l.lambdaj<-list()
                                for(j in 1:J){
                                              temp.lambdaj.tilde <- lambda.tilde.mle[,j]
                                              temp.lambdaj <- fnc.form.Aj.theta5M(temp.lambdaj.tilde)
                                              l.lambdaj[[j]] <- temp.lambdaj
                                }
                                
                                
                                #-----------------------------------------------------------------
                                
                                
                                
                                
                                #  Omega.inv.mle
                                
                                # beta.mle
                                
                                # var.beta.mle<- solve(t(Z)%*%Omega.inv.mle%*%Z)
                                #
                                # # tstat.mle
                                # v.beta.mle<-diag(var.beta.mle); v.beta.mle
                                # #sd.beta.mle<-sqrt(v.beta.mle)
                                # #tstat.mle<-beta.mle/sd.beta.mle
                                
                                
                                #-----------------------------------------------------------------
                                #  gamma.mj
                                
                                lst.gamma<-list()
                                
                                
                                count_begin <- 0
                                count_end <- 0
                                v.mle <- Y - Z %*% beta.mle
                                for (j in 1:J){
                                              
                                              n.j <- n[j]
                                              
                                              count_begin <- count_end + 1 
                                              nrows <- n[j] * M
                                              count_end <- count_end + nrows
                                              
                                        
                                              omega.j.inv <- l.Omega.inv[[j]]
                                              a <- t(as.matrix(l.lambdaj[[j]])) %x% t(iota.nj(n.j))
                                              
                                              gamma.j <- a %*% omega.j.inv %*% v.mle[count_begin:count_end]
                                            
                                              lst.gamma[[j]] <- gamma.j 
                                              
                                              count_begin <- count_end + 1 
                                            
                                
                                }
                                
                                
                                
                                
                                gamma <- t(do.call(cbind, lst.gamma))
                                return(gamma)
                                #-----------------------------------------------------------------
                                # match university ID with gamma.mj's
                                
                                # va.out<-cbind(idj,gamma)
                                # 
                                # #-----------------------------------------------------------------
                                # va.out$ave_cent <- (va.out$gamma.raz+va.out$gamma.lec + va.out$gamma.eng + va.out$gamma.cit + va.out$gamma.wcom)/M
                                # 
                                # #----- weights of the eta's from the standard deviations
                                # 
                                # # turn into zeros the negative entries:
                                # 
                                # 
                                # l.lambdaj <- lapply(l.lambdaj, function(x) ifelse(x < 0, 0, x) )
                                # #dim(lambda.tilde)<-c((M*(M-1)),J)
                                # 
                                # va.out$w_cent<-0
                                # for (j in 1:J)
                                # {
                                #   #va.out$w_cent[j]<- as.matrix(gamma[j,])%*% solve((as.matrix(l.lambdaj[[j]]))^(0.5)) %*% (iota.nj(M)) # if singular, solve() does not work
                                #   va.out$w_cent[j]<- as.matrix(gamma[j,])%*% ginv((as.matrix(l.lambdaj[[j]]))^(0.5)) %*% (iota.nj(M)) # if invertable, solve() and ginv() give the same exact result!!!!
                                # }
                                # 
                                # #output <- va.out[order(va.out[,names(va.out)=="cent"]),]
                                # 
                                # output <- va.out[with(va.out, order(va.out$w_cent, decreasing=TRUE)),]
                                # va.plot<- output
                                # 
                                # # Correlation tables: among univariate VA's
                                # 
                                # out.qr<-out.qr[,c(1,2,4)]; names(out.qr)[3]<-"va.qr"; names(out.qr)[2]<-"nj.qr"
                                # out.cr<-out.cr[,c(1,4)];   names(out.cr)[2]<-"va.cr"
                                # out.en<-out.en[,c(1,4)];   names(out.en)[2]<-"va.en"
                                # out.cc<-out.cc[,c(1,4)];   names(out.cc)[2]<-"va.cc"
                                # out.wrc<-out.wrc[,c(1,4)]; names(out.wrc)[2]<-"va.wrc"
                                # 
                                # 
                                # univar.va<-merge(out.qr, out.cr,by="univ_id")
                                # univar.va<-merge(univar.va, out.en,by="univ_id")
                                # univar.va<-merge(univar.va, out.cc,by="univ_id")
                                # univar.va<-merge(univar.va, out.wrc,by="univ_id")
                                # 
                                # univar.va[1:5,]
                                # univar.cor<-cor(as.matrix(univar.va[,c(3:(M+2))]),method ="spearman")
                                # univar.cor<-round(univar.cor, 3)
                                # 
                                # # Correlation tables: among multivariate and univariate VA's
                                # 
                                # names(univar.va)[1]<-"idj"
                                # all.va<-merge(univar.va, va.out,by="idj")
                                # names(all.va)
                                # allva.cor<-cor(as.matrix(all.va[,c(3:7,10:16)]), method="spearman")
                                # allva.cor<-round(allva.cor,3)
                                # 
                                # # correlation between modules
                                # module.cor<-cor(cbind(data.y$mod_razona_cuantitativo_punt,data.y$mod_lectura_critica, data.y$mod_ingles_punt,data.y$mod_comp_ciudadanas_punt, data.y$mod_comunica_escrita_punt))
                                # module.cor<-round(module.cor,3)

}




#-------------------------------------------------------------------------------
#gamma estimation for univariate-------------------------------------------------------------------------------

Univar_VA_estimation <- function(x, y, n, x2 = NULL) {
  

                  
                  # within estimator 
                  # NOTE: For model 1 the intercept, and for models 2 & 3 both the intercept and the COMP variable turn to zeros in WZ (!)
                 
                  Z.w <- x
                  WZ <- as.matrix(Within2(Z.w,n))
                  WY <- as.matrix(Within2(y,n))                                    # wy = wx*beta + wu
                  pi.w.hat<- ginv(crossprod(WZ)) %*% crossprod(WZ,WY)            
                  ###pi.w.hat
                  
                  # zeta.sq.hat
                  u.w.hat <- y- (Z.w %*% pi.w.hat); #is.matrix(u.w.hat)
                  Wu <- as.matrix(Within2(u.w.hat, n))
                  N <- nrow(Z.w)
                  J <- length(n)
                  K <- ncol(Z.w)
                  zeta.sq.hat <- crossprod(u.w.hat, Wu)/(N - J - K)      ### sigma 2
                  ###zeta.sq.hat; N; J; K
                  
                  # between estimator   
                  # NOTE: Should remove the COMP variable because one of the covariates and its COMP are identical in BZ (!)
                  
                  # define Z.b 
                  
                  Z.b <- cbind(rep(1,sum(n)), x)   ### with intercept 
                  BZ <- as.matrix(Between2(Z.b, n))
                  BY <- as.matrix(Between2(y, n))
                  pi.b.hat<- ginv(crossprod(BZ)) %*% crossprod(BZ,BY)
                  
                  
                  # tau.sq.hat
                  u.b.hat<- BY - (BZ%*%pi.b.hat)
                  n_vector <- as.vector(n)
                  K <- ncol(Z.b)
                  A <- crossprod(u.b.hat)
                  B <- zeta.sq.hat * sum(diag(crossprod(BZ,(diag(1/n_vector)%*%BZ))%*%solve(crossprod(BZ))))  ###B
                  C <- zeta.sq.hat * sum(1/n) ; ###C
                  tau.sq.hat <- (A+B-C)/(J-K)
                  
                  # GLS estimator of pi.hat (Z should include the intercept and the COMP (if acc. to the model))
                  
                  # V.hat
                  tau <- matrix(tau.sq.hat)
                  zeta <- matrix(zeta.sq.hat) 
                  
                  # create a block diagonal matrix: V
                  listV<-list()
                  for(i in c(1:length(n))){
                    temp <- J.nj(n[i]) * tau[1] + Diagonal(n[i]) * zeta[1]
                    listV[[i]] <- solve(temp)
                  }
                  V.hat.inv<-as.matrix(bdiag(listV))

                  # pi.gls.hat
                  if (!is.null(x2)) {
                     Z <- cbind(rep(1,sum(n)),x2)
                  }
                  else{Z <- cbind(rep(1,sum(n)),x)}
                  

                  
                  pi.var <- solve((crossprod(Z,V.hat.inv))%*%Z)
                  pi.sd <- sqrt(diag(pi.var))
                  pi.gls.hat <- pi.var %*% (crossprod(Z,V.hat.inv) %*% y)
                  
                  # value added
                  u.hat<- as.matrix(y-Z%*%pi.gls.hat)  
                  N <- nrow(u.hat)
                  K<- ncol(Z)
                  
                  tau.L.N <- tau[1]*t(L.n(n_vector))
                  va <- tau.L.N%*%V.hat.inv%*%u.hat   ###; va ; N;J;K
                  
                  # link to university id and name:
                  return(va)
                  
                  
  
  
}

va_multiple <- function(x, y, n, x2 = NULL) {
                df_va <- 0
                for (i in c(1 : length(y[1,]))) {
                      y_single <- as.matrix(y[,i])
                      if (!is.null(x2)) {
                        va <- Univar_VA_estimation(x, y_single, n, x2)
                      }
                      else{va <- Univar_VA_estimation(x = x, y = y_single, n = n)}
                      
                      df_va <- cbind(df_va, va)
                }

                return(df_va)
}


