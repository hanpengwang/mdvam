data_process <- function(x, y, col_category) {
                df_list <- list(col_category, y, x)
                col_rename <- c('category', 'y', 'x')
                col_names <- list(NULL, NULL, NULL)
                for (i in c(1:3)) {
                    df_rename <- df_list[[i]]
                    if (is.null(colnames(df_rename)) || length(colnames(df_rename)) != length(df_rename[1,])) 
                      {
                      colnames(df_rename) <- paste(col_rename[i], c(1 : length(df_rename[1,])))
                      
                      
                      
                    }
                    col_names[[i]] <- colnames(df_rename)
                }
                
                df <- cbind(col_category, y, x)
                df <- df[order(df[,1]),]
                df <- as.matrix(na.omit(df))
                col_category <- df[,c(col_names[[1]])]
                y <- df[,c(col_names[[2]])]
                x <- df[,c(col_names[[3]])]
                
                return(list(col_category, y, x))
                
                  
}