## combat2.R - bare bones COMBAT for two datasets
##
## by Artem Sokolov

## COMBAT correction for two datasets
## X1, X2 - two genes-by-samples matrices
combat2 <- function( X1, X2 )
{
    library( sva )
    stopifnot( all( rownames(X1) == rownames(X2) ) )

    ## Create a joint set
    n1 <- ncol(X1)
    n2 <- ncol(X2)
    XX <- cbind( X1, X2 )
    
    ## Remove any genes that have no variance in either dataset
    j1 <- which( apply( X1, 1, var ) == 0 )
    j2 <- which( apply( X2, 1, var ) == 0 )
    jj <- union( j1, j2 )
    if( length(jj) > 0 )
        XX <- XX[-jj,]

    ## Apply COMBAT
    y.b <- rep( c("set1", "set2"), c( n1, n2 ) )
    Xc <- ComBat( XX, y.b, NULL )

    ## Split the data back up
    X1 <- Xc[,1:n1]
    X2 <- Xc[,(n1+1):(n1+n2)]

    list( X1, X2 )
}
