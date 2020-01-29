data_process <- function(x, y, col_category) {
                df_list <- list(col_category, y, x)
                col_rename <- c('category', 'y', 'x')
                col_names <- list(NULL, NULL, NULL)
                for (i in c(1:3)) {
                    df_rename <- df_list[[i]]
                    if (is.null(colnames(df_rename)) || length(colnames(df_rename)) != length(df_rename[1,])) 
                      {
                      colnames(df_rename) <- paste(col_rename[i], c(1 : length(df_rename[1,])), sep = '_')

                      
                      
                    }
                    df_list[[i]] <- df_rename
                    col_names[[i]] <- colnames(df_rename)
                }
                
               
                df <- cbind(df_list[[1]], df_list[[2]], df_list[[3]])
                df <- df[order(df[,1]),]
                df <- as.matrix(na.omit(df))
                col_category <- as.matrix(df[,c(col_names[[1]])])
                y <- as.matrix(df[,c(col_names[[2]])])
                x <- as.matrix(df[,c(col_names[[3]])])
                colnames(col_category) <- col_names[[1]] ; colnames(y) <- col_names[[2]] ; colnames(x) <- col_names[[3]] 
                
                return(list(col_category, y, x))
                
                  
}