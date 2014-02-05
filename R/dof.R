dof <- function(m){
    P <- length(m$baselearner)
    ms <- m$mstop()
    N <- nrow(as.matrix(extract(m$baselearner[[1]], what="design", asmatrix=T, expand=T)))
    H <- vector("list", P)
    for(j in sort(unique(m$xselect()))){
        X <- extract(m$baselearner[[j]], what="design", asmatrix=T, expand=T)
        K <- extract(m$baselearner[[j]], what="penalty", asmatrix=T)
        lambda <- as.numeric(extract(m$baselearner[[j]]$dpp(rep(1, N)), what = "lambda"))
        H[[j]] <- X%*%solve((t(X)%*%X)+lambda*K)%*%t(X)}
    DF <- rep(0, ms)
    nu <- m$control$nu
    backpart <- diag(N) - nu*H[[m$xselect()[1]]]
    DF[1] <- sum(diag(diag(N)-backpart))
    count <- 1
    for(i in 1:ms){
        print(paste("remaining: ", ms-count+1, " steps.", sep=""))
        backpart <- backpart%*%(diag(N) - nu*H[[m$xselect()[i]]])
        DF[i] <- sum(diag(diag(N)-backpart))
        count <- count+1}
    return(DF)}
