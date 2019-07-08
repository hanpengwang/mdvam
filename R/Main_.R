#' Estimate value added 
#'
#'
#' @param x dependent variables
#' @param y indenpent variables, can be one or multiple for each observation
#' @param col_category this is a vector where each value corresponding each observation's id 
#' @param comp_effcet default null. a vector contains variables names which you want to create composition effects
#' @param ... for future implementation
#' 
#' @return This function returns a phrase to print, with or without an
#'    exclamation point added. As a side effect, this function also prints out
#'    the phrase.
#'
#' @examples
#' mvam(x, y, col_category)
#' mvam(x, y, col_category, c('X1','X2',..,'Xm'))
#'
#' @export

### variable user guide
# category : a vector corresponding each id for X and Y
# groups or n : a vector with number of observations for each category
# unique.category : set(category)

mvam <- function(x, y, col_category, comp_effect = NULL, ...){
                  data_contain <- data_process(x, y, col_category)

                  col_category <- data_contain[[1]]
                  y <- data_contain[[2]]
                  x <- data_contain[[3]]

                  df_x <- cbind(col_category, x)

                  idj <- as.matrix(unlist(unique(col_category))) ### schools' names
                  idj <- sort(idj[,1])
                  M <- length(y[1,])
                  n <- as.matrix(table(col_category))
                  N <- sum(n)
                  K <- dim(X)[2]-1
                  J <- length(idj)
                  df <- M*N - M * K - M * J


                  if ( !is.null(comp_effect)) {
                    comp_df <- cbind(col_category , x[,c(comp_effect)])
                    comp_col <- comp_create(comp_df, n)
                    x <- cbind(x, comp_col)

                  }
                  
                  z_stack <- Z.var(x, M, col_category, intercept = T)
                  x_stack <- Z.var(x, M, col_category, intercept = F, within_transform = F)
                  y_stack <- Y.var(y, M, col_category, within_transform = F)


                  wx <- Z.var(x, M, col_category, intercept = F, within_transform = T)
                  wy <- Y.var(y, M, col_category, within_transform = T)

                  params <- sigma_estimation(wx, wy, z_stack, x_stack, y_stack, df)

                  beta <- params[[1]]
                  e <- params[[2]]
                  sigma.sq <- params[[3]]
                  R <- solve(crossprod(z_stack))

                  QH <- QH_(x, M, R, n, col_category)

                  lambda_tilde <- lambda_estimation(x, M, col_category, R, sigma.sq, QH, e)

                  temp_omega <- gls_omega(M, sigma.sq, lambda_tilde, n)
                  omega.inv <- temp_omega[[1]]

                  list_omega <- temp_omega[[2]]


                  beta.gls <- solve(t(z_stack) %*% (omega.inv) %*% z_stack) %*% (t(z_stack) %*% omega.inv %*% y_stack)



                  Gamma_ <- VA_estimation(lambda_tilde, beta.gls, list_omega, J, y_stack, z_stack, n, M)

                  return (cbind(idj, Gamma_))





}

umvam <- function(x, y, col_category, comp_effect = NULL, ...){

                    data_contain <- data_process(x, y, col_category)
                    col_category <- data_contain[[1]]
                    y <- data_contain[[2]]
                    x <- data_contain[[3]]
                    n <- as.matrix(table(col_category))
                    idj <- as.matrix(unlist(unique(col_category)))
                    idj <- sort(idj[,1])
                    if ( !is.null(comp_effect)) {
                        comp_df <- cbind(col_category , x[,c(comp_effect)])
                        comp_col <- comp_create(comp_df, n)
                        x_comp <- cbind(x, comp_col)
                        va_ <- va_multiple(x, y, n,  x_comp)
                      }
                    
                    va_ <- va_multiple(x = x, y = y, n = n)
                    
                    va_[,1] <- idj
                    
                    return(va_)
  
}














