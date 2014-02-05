toolbox <- function(tool, input){
  if(!(tool %in% c("hi.for.log.pl", "predict.f.x"))){
    stop("please select an available tool!", call. = FALSE)
  }
  ##########################################
  ## helpindex for log partial likelihood ##
  ##########################################
  if(tool == "hi.for.log.pl"){
    if(!(all(names(input) == c("possible.trans.from.and.to", "d")))){
      stop("for tool 'hi.for.log.pl', please provide input argument list containing 'possible.trans.from.and.to' and 'd'!", 
           call. = FALSE)
    }
    possible.trans.from.and.to <- input$possible.trans.from.and.to
    d <- input$d
    hi.for.log.pl <- rep(NA, nrow(d))
    for(hi in 1:(nrow(d))){
      if((d$to[hi] != 99) & (d$from[hi] != 99)){
        hi.for.log.pl[hi] <- which((possible.trans.from.and.to[, 1] == d$from[hi]) &
                                    (possible.trans.from.and.to[, 2] == d$to[hi]))
      }
    }
    output <- list(hi.for.log.pl = hi.for.log.pl)
  }
  #################################################
  ## predict (non)linear covariate effect f of x ##
  #################################################
  if(tool == "predict.f.x"){
    if(!(all(names(input) == c("x", "trans", "mod", "d", "seq.length")))){
      stop("for tool 'predict.f.x', please provide input argument list containing 'x', 'trans', 'mod', 'd', 'seq.length'!", 
           call. = FALSE)
    }
    x <- input$x
    trans <- input$trans
    mod <- input$mod
    d <- input$d
    seq.length <- input$seq.length
    f.x <- rep(0, seq.length)
    for(bl in 1:length(mod$baselearner)){
      ding <- mod$baselearner[[bl]]$get_names()
      ## (i) transition specific baselearner:
      if(length(ding) > 1.5){
        if(x == ding[1]){
          ding <- strsplit(ding[2], "_")[[1]]
          if(any(trans == ding[2:length(ding)])){
            ho <- NULL
            if(any(bl %in% unique(mod$xselect()))){
              var.names <- mod$baselearner[[bl]]$get_names()
              range.x <- range(d[, var.names[1]])
              newd <- data.frame(x = seq(range.x[1], range.x[2], length=seq.length))
              names(newd)[1] <- var.names[1]
              if(length(var.names)>1.5){
                newd$tr <- TRUE
                names(newd)[2] <- var.names[2]}
              ho <- predict(mod, which=bl, newd=newd, aggregate="cumsum")[[1]]
            }else{
              ho <- matrix(nrow=seq.length, ncol=mod$mstop(), data=0)
            }
            f.x <- f.x + ho
          }
        }
        ## (ii) baselearner parameterized over all transitions:
      }else{
        if(x == ding[1]){
          ho <- NULL
          if(any(bl %in% unique(mod$xselect()))){
            var.names <- mod$baselearner[[bl]]$get_names()
            range.x <- range(d[, var.names[1]])
            newd <- data.frame(x = seq(range.x[1], range.x[2], length=seq.length))
            names(newd)[1] <- var.names[1]
            if(length(var.names)>1.5){
              newd$tr <- TRUE
              names(newd)[2] <- var.names[2]}
            ho <- predict(mod, which=bl, newd=newd, aggregate="cumsum")[[1]]
          }else{
            ho <- matrix(nrow=seq.length, ncol=mod$mstop(), data=0)
          }
          f.x <- f.x + ho
        }
      }
    }
    range.x <- range(d[, x])
    output = list(y=f.x,
                  x=seq(range.x[1], range.x[2], length=seq.length))
  }
  ##########
  ## done ##
  ##########
  return(output)
}
