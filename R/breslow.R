breslow <- function(m, mstop, entry, exit, trans, event){
    f <- predict(m, aggregate="cumsum")
    f <- f[, mstop]
    n <- length(entry)
    dummy <- rep(0, n)
    ef <- exp(f)
    for (j in 1:n){
        dummy[j] <- sum(as.numeric((entry < exit[j]) & (exit[j] <= exit) & (trans[j] == trans)) * ef)
    }
    bhr <- rep(0, n)
    for(i in 1:n){
        dummy[which(dummy==0.0)] <- 1.0
        hi <- which((event > 0.5) & (entry[i] < exit) & (exit <= exit[i]) & (trans[i] == trans))
        bhr[i] = sum(1/dummy[hi])
    }
    Q <- sort(unique(trans))
    A <- vector("list", length(Q))
    for(q in 1:length(Q)){
        A[[q]]$times <- exit[which(trans==Q[q])]
        A[[q]]$bhr <- bhr[which(trans==Q[q])]}
    names(A) <- paste("t", Q, sep="")
    return(A)}
