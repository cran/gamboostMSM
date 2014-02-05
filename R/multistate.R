multistate <- function(trans = trans){
    ## #########################################################
    ## negative log likelihood of a cox type multistate model ##
    ## a.k.a. partial likelihood loss function                ##
    ## #########################################################
    plloss <- function(y, f, w) {
        entry <- y[, 1]
        exit <- y[, 2]
        event <- y[, 3]
        n <- length(entry)
        if (length(f) == 1) {f <- rep(f, n)}
        ef <- exp(f)
        risk <- rep(0, n)
        for (i in 1:n){risk[i] <- sum(as.numeric((entry < exit[i]) & (exit[i] <= exit) & (trans[i] == trans)) * ef)}
        lpl <- event * (f - log(risk))
        return(lpl)
    }
    ## #################################
    ## construct new boosting family ##
    ## #################################
    Family(
    ## #########################################
    ## negative gradient of the loss function ##
    ## #########################################
    ngradient = function(y, f, w) {
        entry <- y[, 1]
        exit <- y[, 2]
        trans = trans
        event <- y[, 3]
        n <- length(entry)
        if (length(f) == 1){f <- rep(f, n)}
        dummy <- rep(0, n)
        ef <- exp(f)
        ## nenner fuer jedes j ausrechnen (in jedem element summer ueber k)
        for (j in 1:n){
            dummy[j] <- sum(as.numeric((entry < exit[j]) & (exit[j] <= exit) & (trans[j] == trans)) * ef)
        }
        zi <- rep(0, n)
        for(i in 1:n){
            dummy[which(dummy==0.0)] <- 1.0
            ## riskset im zaehler ausrechnen:
            hi <- which((event > 0.5) & (entry[i] < exit) & (exit <= exit[i]) & (trans[i] == trans))
            zi[i] = sum(1/dummy[hi])
            ## for(j in 1:n){
            ##    if((event[j]>0.5) & (entry[i]<exit[j]) & (exit[j]<=exit[i]) & (trans[i]==trans[j])){
            ##        helpobject <- dummy[j]
            ##        if(helpobject == 0.0){
            ##            helpobject <- 1.0
            ##        }
            ##        helpobject <- ef[i]/helpobject
            ##        zi[i] = zi[i] + helpobject
            ##    }
            ## }
        }
        zi <- event - (ef * zi)
        return(zi)},
    ## ################
    ## risk function ##
    ## ################
    risk = risk <- function(y, f, w) -sum(plloss(y, f, w)),
    offset = function(y, w) 0,
    name = "family for multistate models.")}
