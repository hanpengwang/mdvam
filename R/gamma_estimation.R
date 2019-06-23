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


VA_estimation <- function(lambda.tilde.mle, beta.mle, l.Omega.inv, J, Y, Z, n){
                                
 
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
                                
                                count_begin <- 1
                                
                                v.mle.j <- Y - Z %*% beta.mle
                                browser()
                                for (j in 1:J){
                                
                                              n.j <- n[j]
                                              count_end <- count_begin + n.j - 1 
                                        
                                              omega.j.inv <- l.Omega.inv[[j]]
                                              a <- t(as.matrix(l.lambdaj[[j]])) %x% t(iota.nj(n.j))
                                              
                                              gamma.j <- a %*% omega.j.inv %*% v.mle.j[count_begin : count_end]
                                            
                                              lst.gamma[[j]] <- gamma.j
                                              
                                              count_begin <- count_end + 1 
                                            
                                
                                }
                                
                                
                                
                                
                                gamma <- do.call(cbind, lst.gamma)
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

