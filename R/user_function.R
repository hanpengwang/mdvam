#' Print "Hello world"
#'
#' This is a simple function that, by default, prints "Hello world". You can
#' customize the text to print (using the \code{to_print} argument) and add
#' an exclamation point (\code{excited = TRUE}).
#'
#' @param to_print A character string giving the text the function will print
#' @param excited Logical value specifying whether to include an exclamation
#'    point after the text
#'
#' @return This function returns a phrase to print, with or without an
#'    exclamation point added. As a side effect, this function also prints out
#'    the phrase.
#'
#' @examples
#' hello_world()
#' hello_world(excited = TRUE)
#' hello_world(to_print = "Hi world")
#'
#' @export

### variable user guide
# category : a vector corresponding each id for X and Y
# groups or n : a vector with number of observations for each category
# unique.category : set(category)

mvam <- function(x, y, col_category){

                  idj <- as.matrix(unlist(unique(col_category))) ### schools' names

                  M <- length(y[1,])
                  n <- as.matrix(table(col_category))
                  N <- sum(n)
                  K <- dim(X)[2]-1
                  J <- length(idj)
                  df <- M*N - M * K - M * J

                  z_stack <- Z.var(x, M, col_category, intercept = T)
                  x_stack <- Z.var(x, M, col_category, intercept = F, within_transform = F)
                  y_stack <- Y.var(y, M, col_category, within_transform = F)



                  wx <- Z.var(x, M, col_category, intercept = F, within_transform = T)
                  wy <- Y.var(y, M, col_category, within_transform = T)


                  w <-  Within(M, n)


                  params <- sigma_estimation(wx, wy, z_stack, x_stack, y_stack, df)

                  beta <- params[[1]]
                  e <- params[[2]]
                  sigma.sq <- params[[3]]
                  R <- solve(crossprod(z_stack))

                  QH <- QH_estimation(x, M, R, n, col_category)

                  lambda_tilde <- lambda_estimation(x, M, col_category, R, sigma.sq, QH, e)

                  omega.inv <- gls_omega(M, sigma.sq, lambda_tilde, n)

                  #return (list(lambda_tilde, omega))




                  #add



}
