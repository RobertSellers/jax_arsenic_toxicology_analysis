standrm_harmony <- function(formula, data, fct, curveid=NULL, random=NULL, priors=standrc_priors(), ...){
  callDetail <- match.call()
  
  mf <- model.frame(formula, data)
  
  if (!is.null(curveid)){
    cid <- as.character(curveid)[3]
    idc <- as.numeric(data$individual_id)
    J <- length(unique(idc))
    stsp <- strsplit(as.character(curveid)[2], "+")[[1]]
    cin <- stsp[!stsp %in% c(" ", "+")]
    pnl <- sapply(fct$names, function(x) x %in% cin)
    curvenames <- levels(as.factor(data$individual_id))
  } else {
    cid <- idc <- NULL
    J <- 1
    pnl <- rep(FALSE, length(fct$names))
    curvenames <- NULL
  }
  
  if (!is.null(random)){
    rid <- as.character(random)[3]
    idr <- as.numeric(data$dir_id)
    K <- length(unique(idr))
    stsp <- strsplit(as.character(random)[2], "+")[[1]]
    rin <- stsp[!stsp %in% c(" ", "+")]
    pnlr <- sapply(c("b", "c", "d", "e", "f"), function(x) x %in% rin)
  } else {
    rid <- idr <- NULL
    K <- 1
    pnlr <- rep(FALSE, 5)
  }
  
  fix <- fct$fixed
  if (fct$name %in% c("W1.4", "W1.3", "W2.4", "W2.3", "LN.4", "LN.3") | attr(fct, "class") == "fp-logistic") fix <- c(fix, 0)  
  isfix <- !is.na(fix)
  
  jv <- rep(J, 5)
  jv[!isfix][!pnl] <- 1
  
  N <- nrow(mf) 
  y <- mf[,1]
  x <- mf[,2]   
  if (fct$name %in% c("LL.5", "LL.4", "LL.3", "W1.4", "W1.3", "W2.4", "W2.3", "LN.4", "LN.3", "MM.2", "MM.3") | attr(fct, "class") == "fp-logistic") x[x == 0] <- 0.5*min(x[x > 0])
  if (fct$name %in% c("LL.4", "LL.3", "L.4", "L.3", "MM.2", "MM.3") | attr(fct, "class") == "fp-logistic") fix[5] <- 0
  if (fct$name %in% c("MM.2", "MM.3")){
    fix[1] <- 0
    fct$name <- "LL.5"
  } 
  if (attr(fct, "class") == "fp-logistic"){
    p1 <- get("p1", environment(fct$fct))
    p2 <- get("p2", environment(fct$fct))    
  }
  
  if (is.null(priors$pb)) pb <- rep(0, jv[1]) else pb <- priors$pb
  if (mean(y[x == min(x)]) < mean(y[x == max(x)])){
    if (is.null(priors$pc)) pc <- rep(min(y), jv[2]) else pc <- priors$pc
    if (is.null(priors$pd)) pd <- rep(max(y), jv[3]) else pd <- priors$pd
  } else {
    if (is.null(priors$pd)) pd <- rep(min(y), jv[3]) else pd <- priors$pd
    if (is.null(priors$pc)) pc <- rep(max(y), jv[2]) else pc <- priors$pc
  }  
  if (is.null(priors$pe)) pe <- rep(median(x), jv[4]) else pe <- priors$pe
  if (is.null(priors$pf)) pf <- rep(0, jv[5]) else pf <- priors$pf
  
  stan_dat <- list(N=N, J=J, K=K, idc=idc, idr=idr, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  if (is.null(curveid)) stan_dat <- list(N=N, K=K, idr=idr, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  if (is.null(random)) stan_dat <- list(N=N, J=J, idc=idc, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  if (is.null(curveid) & is.null(random)) stan_dat <- list(N=N, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  
  # data 
  dJ <- "int<lower=0> idc[N]; int<lower=0> J;" 
  dK <- "int<lower=0> idr[N]; int<lower=0> K;" 
  
  # priors
  dJp <- c("real pb[J];", "real pc[J];", "real pd[J];", "real pe[J];", "real pf[J];")[!isfix]
  dp <- c("real pb;", "real pc;", "real pd;", "real pe;", "real pf;")[!isfix]  
  dpc <- paste(sapply(1:length(pnl), function(i) if (pnl[i]) dJp[i] else dp[i]), collapse=" ")
  
  # define parameters
  paraJ <- c("real slope[J];", "real lasy[J];", "real uasy[J];", "real<lower=0> ed[J];", "real assym[J];")[!isfix]
  para <- c("real slope;", "real lasy;", "real uasy;", "real<lower=0> ed;", "real assym;")[!isfix]
  parac <- paste(sapply(1:length(pnl), function(i) if (pnl[i]) paraJ[i] else para[i]), collapse=" ")
  
  # random effects
  rpara <- c("real rslope[K];", "real rlasy[K];", "real ruasy[K];", "real red[K];", "real rassym[K];")
  rparasig <- c("real<lower=0> sigmasq_slope;", "real<lower=0> sigmasq_lasy;", "real<lower=0> sigmasq_uasy;", "real<lower=0> sigmasq_ed;", "real<lower=0> sigmasq_assym;")
  rparac <- paste(c(rpara[pnlr], rparasig[pnlr]), collapse=" ")
  
  trap <- c("slope", "lasy", "uasy", "ed", "assym")[!isfix]
  traj <- c("slope[idc[i]]", "lasy[idc[i]]", "uasy[idc[i]]", "ed[idc[i]]", "assym[idc[i]]")[!isfix]
  tra <- character(length(fix))
  tra[!isfix] <- sapply(1:length(pnl), function(i) if (pnl[i]) traj[i] else trap[i])
  tra[isfix] <- fix[isfix]
  if (is.null(random)){
    if (fct$name %in% c("LL.5", "LL.4", "LL.3")) trans <- paste("sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") / (1 + exp(-exp(", tra[1], ") * (log(x[i]/", tra[4], "))))^exp(", tra[5], ");}", collapse="")
    if (fct$name %in% c("L.5", "L.4", "L.3")) trans <- paste("sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") / (1 + exp(-exp(", tra[1], ") * (x[i] - ", tra[4], ")))^exp(", tra[5], ");}", collapse="")
    if (fct$name %in% c("W1.4", "W1.3")) trans <- paste("sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") * exp(-exp(-exp(", tra[1], ") * (log(x[i]) - log(", tra[4], "))));}", collapse="")
    if (fct$name %in% c("W2.4", "W2.3")) trans <- paste("sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") * (1 - exp(-exp(-exp(", tra[1], ") * (log(x[i]) - log(", tra[4], ")))));}", collapse="")
    if (fct$name %in% c("LN.4", "LN.3")) trans <- paste("sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") * normal_cdf(exp(", tra[1], ") * (log(x[i]) - log(", tra[4], ")), 0, 1);}", collapse="")
    if (attr(fct, "class") == "fp-logistic") trans <- paste("sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") / (1 + exp(-exp(", tra[1], ") * log(x[i] + 1)^", p1," + ", tra[4], "* log(x[i] + 1)^", p2,"));}", collapse="")
  } else {
    stra <- paste(c("real<lower=0> sigma_slope;", "real<lower=0> sigma_lasy;", "real<lower=0> sigma_uasy;", "real<lower=0> sigma_ed;", "real<lower=0> sigma_assym;")[pnlr], collapse=" ")
    strasq <- paste(c("sigma_slope <- sqrt(sigmasq_slope);",
                      "sigma_lasy <- sqrt(sigmasq_lasy);",
                      "sigma_uasy <- sqrt(sigmasq_uasy);",
                      "sigma_ed <- sqrt(sigmasq_ed);",
                      "sigma_assym <- sqrt(sigmasq_assym);")[pnlr], collapse=" ")
    trar <- c(paste("(", tra[1], " + rslope[idr[i]] )"), 
              paste("(", tra[2], " + rlasy[idr[i]] )"), 
              paste("(", tra[3], " + ruasy[idr[i]] )"), 
              paste("(", tra[4], " + red[idr[i]] )"), 
              paste("(", tra[5], " + rassym[idr[i]] )"))
    trc <- sapply(1:length(tra), function(i) if (pnlr[i]) trar[i] else tra[i])    
    if (fct$name %in% c("LL.5", "LL.4", "LL.3")) trans <- paste(stra, strasq, "sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") / (1 + exp(-exp(", trc[1], ") * (log(x[i]/", trc[4], "))))^exp(", trc[5], ");}", collapse="")
    if (fct$name %in% c("L.5", "L.4", "L.3")) trans <- paste(stra, strasq, "sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") / (1 + exp(-exp(", trc[1], ") * (x[i] -", trc[4], ")))^exp(", trc[5], ");}", collapse="")
    if (fct$name %in% c("W1.4", "W1.3")) trans <- paste(stra, strasq, "sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") * exp(-exp(-exp(", trc[1], ") * (log(x[i]) - log(", trc[4], "))));}", collapse="")
    if (fct$name %in% c("W2.4", "W2.3")) trans <- paste(stra, strasq, "sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") * (1 - exp(-exp(-exp(", trc[1], ") * (log(x[i]) - log(", trc[4], ")))));}", collapse="")
    if (fct$name %in% c("LN.4", "LN.3")) trans <- paste(stra, strasq, "sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") * normal_cdf(exp(", trc[1], ") * (log(x[i]) - log(", trc[4], ")), 0, 1);}", collapse="")
    if (attr(fct, "class") == "fp-logistic") trans <- paste("sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") / (1 + exp(-exp(", trc[1], ") * log(x[i] + 1)^", p1," + ", trc[4], "* log(x[i] + 1)^", p2,"));}", collapse="")
  }
  
  ### population parameters
  popJ <- c("real pslope[J];", "real plasy[J];", "real puasy[J];", "real ped[J];", "real passym[J];")
  pop <- c("real pslope;", "real plasy;", "real puasy;", "real ped;", "real passym;")
  pop[!isfix] <- sapply(1:length(pnl), function(i) if (pnl[i]) popJ[!isfix][i] else pop[!isfix][i]) 
  popc <- paste(pop, collapse=" ")
  rnpara <- c("real rnslope;", "real rnlasy;", "real rnuasy;", "real rned;", "real rnassym;")
  rnparac <- paste(rnpara[pnlr], collapse=" ")
  ppp <- c("slope", "lasy", "uasy", "ed", "assym")[!isfix]
  ppj <- c("slope[j]", "lasy[j]", "uasy[j]", "ed[j]", "assym[j]")[!isfix]
  pp <- character(length(fix))
  pp[!isfix] <- sapply(1:length(pnl), function(i) if (pnl[i]) ppj[i] else ppp[i])
  pp[isfix] <- fix[isfix]
  po <- c("pslope", "plasy", "puasy", "ped", "passym")
  poj <- paste("for (j in 1:J)", c("pslope[j]", "plasy[j]", "puasy[j]", "ped[j]", "passym[j]"))
  pof <- character(length(fix))
  pof[!isfix] <- sapply(1:length(pnl), function(i) if (pnl[i]) poj[!isfix][i] else po[!isfix][i])
  pof[isfix] <- po[isfix]
  ppc <- paste(pof[1], " <- ", pp[1], if (pnlr[1]) " + rnslope", ";", pof[2], "<-", pp[2], if (pnlr[2]) " + rnlasy", ";", pof[3], " <- ", pp[3], if (pnlr[3]) " + rnuasy", ";", pof[4], " <- ", pp[4], if (pnlr[4]) " + rned", ";", pof[5], " <- ", pp[5], if (pnlr[5]) " + rnassym", ";")
  
  modp <- paste(c("rnslope ~ normal(0, sigma_slope);", "rnlasy ~ normal(0, sigma_lasy);", "rnuasy ~ normal(0, sigma_uasy);", "rned ~ normal(0, sigma_ed);", "rnassym ~ normal(0, sigma_assym);")[pnlr], collapse=" ") 
  
  
  #### model statements
  mod <- paste(c(paste("slope ~", priors$b, ";"), 
                 paste("lasy ~", priors$c, ";"), 
                 paste("uasy ~", priors$d, ";"), 
                 paste("ed ~", priors$e, ";"), 
                 paste("assym ~", priors$f, ";"))[!isfix], collapse=" ") 
  
  mody <- paste("sigmasq_y ~", priors$sy, "; y ~ normal(mu, sigma_y);")
  
  modr <- paste(c("rslope ~ normal(0, sigma_slope);", "rlasy ~ normal(0, sigma_lasy);", "ruasy ~ normal(0, sigma_uasy);", "red ~ normal(0, sigma_ed);", "rassym ~ normal(0, sigma_assym);")[pnlr], collapse=" ") 
  modrsig <- paste(c(paste("sigmasq_slope ~", priors$sb, ";"), 
                     paste("sigmasq_lasy ~", priors$sc, ";"), 
                     paste("sigmasq_uasy ~", priors$sd, ";"), 
                     paste("sigmasq_ed ~", priors$se, ";"), 
                     paste("sigmasq_assym ~",priors$sf, ";"))[pnlr], collapse=" ") 
  
  # new random effects for derived ed
  moded <- paste(c("rnslope ~ normal(0, sigma_slope);", "rnlasy ~ normal(0, sigma_lasy);", "rnuasy ~ normal(0, sigma_uasy);", "rned ~ normal(0, sigma_ed);", "rnassym ~ normal(0, sigma_assym);")[pnlr], collapse=" ") 
  
  stancode <- paste("data {", 
                    "int<lower=0> N; real y[N];",
                    "real<lower=0> x[N];", 
                    if (!is.null(idc)) dJ,
                    if (!is.null(random)) dK,
                    dpc,
                    "} ",
                    "parameters { real<lower=0> sigmasq_y;",
                    parac,
                    if (!is.null(random)) rparac, 
                    if (!is.null(random)) rnparac,  
                    "} ",
                    "transformed parameters {",
                    "real<lower=0> sigma_y;",
                    "real mu[N];",                   
                    trans,                      
                    "} ",
                    "model {",
                    mod,
                    mody,                     
                    if (!is.null(random)) modr,
                    if (!is.null(random)) modrsig,    
                    if (!is.null(random)) modp,                    
                    "}",
                    "generated quantities {",
                    "real residuals[N];", 
                    "real log_lik[N];",
                    popc,
                    "for (i in 1:N){",
                    "residuals[i] <- y[i] - mu[i];",
                    "log_lik[i] <- normal_log(y[i], mu[i], sigma_y);",
                    "}",
                    ppc,                                     
                    "}",
                    sep="")
  # save example data
  # fileConn<-file("~/Desktop/projects/jax_arsenic_toxicology_analysis/temp_files/stan_example_2.stan")
  # writeLines(stancode, fileConn, sep="\n")
  # close(fileConn)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  fit <- stan(model_code = stancode, 
              model_name = "test_arsenic",
              data = stan_dat
              , ...)
  out <- list()
  out$call <- callDetail
  out$data <- stan_dat
  out$model <- stancode
  out$stan <- fit
  out$fixedpars <- trap
  out$pars <- c("pslope", "plasy", "puasy", "ped", "passym")
  if (!is.null(random)) out$random <- c("slope", "lasy", "uasy", "ed", "assym")[pnlr] 
  out$curves <- list(pars=pnl, J=J, names=curvenames)
  out$fixed <- fix
  out$fct <- fct
  class(out) <- "standrc"
  return(out)
}


robert_ED <- function(object, ..., respLev=NULL){ 
  if (is.null(respLev)) stop("Please provide response levels!")
  plist <- extract(object$stan, pars=object$pars)
  fix <- object$fixed
  fid <- logical(length=length(fix))
  fid[is.na(fix)] <- object$curves$pars
  samp <- lapply(1:object$curves$J, function(i){
    sapply(1:length(fix), function(j){
      if (fid[j]) plist[[j]][,i] else plist[[j]]
    })
  })
  EDlist <- lapply(respLev, function(p){
    smat <- cbind(sapply(1:length(samp), function(i){
      apply(samp[[i]], 1, function(x){
        xt <- x
        xt[1] <- -exp(xt[1])
        xt[5] <- exp(xt[5])
        p <- 100-p
        if (object$fct$name %in% c("LL.5", "LL.4", "LL.3")){
          tempVal <- log((100 - p)/100)
          value <- xt[4] * (exp(-tempVal/xt[5]) - 1)^(1/xt[1])
        }
        if (object$fct$name %in% c("L.5", "L.4", "L.3")){
          tempVal <- 100/p
          value <- xt[4] + (log(tempVal^(1/xt[5]) - 1))/xt[1]
        }
        if (object$fct$name %in% c("W1.4", "W1.3")){
          tempVal <- log(-log((100 - p)/100))
          value <- exp(tempVal/xt[1] + log(xt[4]))
        }
        if (object$fct$name %in% c("W2.4", "W2.3")){
          p <- 100 - p
          tempVal <- log(-log((100 - p)/100))
          value <- exp(tempVal/xt[1] + log(xt[4]))
        }
        if (object$fct$name %in% c("LN.4", "LN.3")){
          tempVal <- 1 - (100-p)/100
          value <- exp(qnorm(tempVal)/xt[1]) * xt[4]
        }
        return(value)
      })
    }))
    colnames(smat) <- object$curves$names
    return(smat)
  })
  names(EDlist) <- respLev
  class(EDlist) <- "EDsamp"
  return(EDlist)  
}

print.EDsamp <- function(x, ...){
  out <- lapply(x, function(xl){
    dat <- as.data.frame(t(apply(cbind(xl), 2, function(xc) quantile(xc, c(0.5, 0.025, 0.975)))))
    colnames(dat) <- c("median", "2.5%", "97.5%")
    return(dat)
  })
  print(out)
}

plot_new <- function(x, ..., ndose=25, logx=FALSE, lim=NULL){

    if (is.null(x$data$total)) xcc <- as.character(x$call$formula[[2]]) else xcc <- as.character("p")
  if (x$curves$J > 1){
    dframe <- data.frame(x$data$y, x$data$x, x$curves$names[x$data$idc])
    names(dframe) <- c(xcc,
                       as.character(x$call$formula[[3]]), 
                       as.character(x$call$curveid[[3]]))
  } else {
    dframe <- data.frame(x$data$y, x$data$x)
    names(dframe) <- c(xcc,
                       as.character(x$call$formula[[3]]))
  }
  if (!is.null(x$data$total)) dframe[,1] <- dframe[,1]/x$data$total
  
  if (x$curves$J > 1){
    curvn <- paste(",", as.character(x$call$curveid[[3]]), "=x$curves$names")
    dframe[,3] <- factor(dframe[,3], levels=x$curves$names)
  } else {
    curvn=NULL
  }
  if (is.null(lim)){
    minx <- min(x$data$x)
    maxx <- max(x$data$x)
  } else {
    if (length(lim) != 2) stop("Please provide limits as vector with 2 elements.")
    minx <- lim[1]
    maxx <- lim[2]
  }
  if (logx){    
    newd <- eval(parse(text=paste("expand.grid(", as.character(x$call$formula[[3]]), "=exp(seq(log(", minx, "), log(", maxx, "), length=ndose))", curvn, ")")))
  } else {
    newd <- eval(parse(text=paste("expand.grid(", as.character(x$call$formula[[3]]), "=seq(", minx, ", ", maxx, ", length=ndose)", curvn, ")")))    
  }
  pm <- predict(x, newdata=newd)
  newd$p <- apply(pm, 2, function(x) mean(x, na.rm=TRUE))
  # newd$pmin <- apply(pm, 2, function(x) quantile(x, 0.025, na.rm=TRUE))
  # newd$pmax <- apply(pm, 2, function(x) quantile(x, 0.975, na.rm=TRUE))
  
  if (logx) lxt <- "+ coord_trans(x='log')" else lxt <- NULL
#browser()
  chard <-  paste(ggplot(dframe, aes(x= dose , y= response , colour=  )) +  geom_point() + geom_line(data=newd, aes(y=p)))
  #eval(parse(text=paste("ggplot(dframe, aes(x=",  as.character(x$call$formula[[3]]),", y=", xcc, ", colour=", as.character(x$call$curveid[[3]]),")) +  geom_point() + geom_line(data=newd, aes(y=p))", lxt)) ) 
return (list("command" =chard ,
             "dframe" = dframe,
             "newd" = newd))
}

predict_harmony <- function(object, ..., newdata=NULL){ 
  # browser()
  x <- newdata[,as.character(object$call$formula[3])]
  idc <- newdata[,as.character(object$call$curveid)[3]]
  
  plist <- extract(object$stan, pars=object$pars)
  # this should have fixed upper and lower?
  fix <- object$fixed
  fid <- logical(length=length(fix))
  # will have (for b,c,d,e,f), TRUE FALSE FALSE TRUE FALSE
  fid[is.na(fix)] <- object$curves$pars
  samp <- lapply(1:object$curves$J, function(i){
    sapply(1:length(fix), function(j){
      if (fid[j]) plist[[j]][,i] else plist[[j]]
    })
  })

  pred <- sapply(1:length(x), function(i){
    apply(samp[[as.numeric(idc)[i]]], 1, function(xp){
      # this is custom function for LL.5, LL.4, LL.3
      # check predict standrc for additional functions
      return(xp[2] + (xp[3] - xp[2]) / (1 + exp(-exp(xp[1]) * (log(x[i]/ xp[4]))))^exp(xp[5]))
    })
  })
  return(pred)
}

predict_harmony_stan <- function(object, ..., newdata=NULL){ 
  
  # extract pred columns
  x <- newdata[,as.character(object$call$formula[3])]
  idc <- newdata[,as.character(object$call$curveid)[3]]

  plist <- rstan::extract(object$stan, pars=object$pars)
  fix <- object$fixed
  fid <- logical(length=length(fix))
  fid[is.na(fix)] <- object$curves$pars
  samp <- lapply(1:object$curves$J, function(i){
    sapply(1:length(fix), function(j){
      if (fid[j]) plist[[j]][,i] else plist[[j]]
    })
  })

  pred <- sapply(1:length(x), function(i){
    apply(samp[[as.numeric(idc)[i]]], 1, function(xp){
      # this is LL.4() specific extraction
      return(xp[2] + (xp[3] - xp[2]) / (1 + exp(-exp(xp[1]) * (log(x[i]/ xp[4]))))^exp(xp[5]))
    })
  })
  return(pred)
}

### Custom extraction functions


prepare_ed50s_drc <- function(model_, original_data, type='none', pool = 'FALSE'){
  # extract ed50 with whatever default calculation
  ed50_drc <- as.data.frame(drc::ED(model_, 50, 
                                    interval = type,
                                    display = FALSE, 
                                    lref = 0,
                                    uref = 5,
                                    pool=pool, 
                                    type="absolute")
                            )
  # regex the ugly rownames
  ed50_drc$individual <- str_match(rownames(ed50_drc), ":\\s*(.*?)\\s*:")[,2]
  # for output
  df = data.frame(
    "individual_id" = as.character(ed50_drc$individual), 
    "dose" = ed50_drc$Estimate
    )
  # retrieve the response per ed50 estimation with the model
  df$pred <- predict(model_, data.frame(
    dose = ed50_drc$Estimate,
    individual_id = str_match(rownames(ed50_drc), ":\\s*(.*?)\\s*:")[,2]
  ))
  df$group <- "drc"
  df$ed50 <- df$dose
  #df$dose <- NULL
  return (df[,c('individual_id','ed50','pred','group')])
}
# returns individual id to original value
map_do_to_original <- function(df_new, df_org){
  lookup_df <- df_org[,c('individual','individual_id')]
  return(plyr::join(df_new, distinct(lookup_df), by = "individual_id"))
}

prepare_ed50s_stan <- function(model_, original_data){
  ed50_stan <- ED(model_, respLev=c(50))
  df = data.frame("individual_id" = sort(unique(original_data$individual_id)), "dose" = colMedians(ed50_stan$`50`))
  pred_stan <- predict_harmony_stan(model_, newdata = df)
    df <- df %>% 
  rename(
    ed50 = dose
    )
  df$pred <- apply(pred_stan, 2, function(x) mean(x, na.rm=TRUE))
  df$group <- "stan"
  return (df)

}

plot_posterior_growth_ci <- function(fit, dat, cred = TRUE){
  # browser()
  Y_mean <- rstan::extract(fit, "Y_mean")
  Y_mean_mean <- apply(Y_mean$Y_mean, 2, mean)
  Y_pred <- rstan::extract(fit, "Y_pred")
  Y_pred_mean <- apply(Y_pred$Y_pred, 2, mean)
  
  if (cred){
    Y_mean_cred <- apply(Y_mean$Y_mean, 2, quantile, c(0.05, 0.95))
    Y_pred_cred <- apply(Y_pred$Y_pred, 2, quantile, c(0.05, 0.95))
    plot(dat$y ~ dat$x, xlab="x", ylab="Y", 
       ylim=c(-1, 101), main="Non-linear Growth Curve")
    lines(dat$x, Y_mean_mean)
    points(dat$x, Y_pred_mean, pch=19)
    lines(dat$x, Y_mean_cred[1,], col=4)
    lines(dat$x, Y_mean_cred[2,], col=4)
    lines(dat$x, Y_pred_cred[1,], col=2)
    lines(dat$x, Y_pred_cred[2,], col=2)
    legend(x="bottomright", bty="n", lwd=2, lty=c(NA, NA, 1, 1,1),
           legend=c("observation", "prediction", "mean prediction",
                    "90% mean cred. interval", "90% pred. cred. interval"),
           col=c(1,1,1,4,2),  pch=c(1, 19, NA, NA, NA))
  }else{
    plot(dat$y ~ dat$x, xlab="x", ylab="Y", 
       ylim=c(-1, 101), main="Non-linear Growth Curve")
    lines(dat$x, Y_mean_mean)
    points(dat$x, Y_pred_mean, pch=19)
    legend(x="bottomright", bty="n", lwd=2, lty=c(NA, NA, 1, 1,1),
           legend=c("observation", "prediction", "mean prediction"),
           col=c(1,1,1,4,2),  pch=c(1, 19, NA, NA, NA))
  }
}

build_standrc_object_single <- function(stan_dat, stancode, fit){
  out <- list()
  out$data <- stan_dat
  out$model <- stancode
  out$stan <- fit
  out$fixedpars <- c("slope","lasy","uasy","ed")
  out$pars <- c("pslope", "plasy", "puasy", "ped", "passym")
  out$curves <- list(pars=c(FALSE,FALSE,FALSE,FALSE),
                     J=1,
                     names=NULL
                     )
  out$fixed <- c(NA,NA,NA,NA,0)
  out$fct <- LL.4()
  class(out) <- "standrc"
  return(out)
}

prepare_data_for_stan <- function(ft, dat, log_ = FALSE, scale_ = TRUE){
  dat$response <- as.numeric(unlist(dat[,c(ft)]))
  dat <- dat[!is.na(dat$response), ] 
  if (log_) {
      dat$response <- log(dat$response)
  }
  if (scale_){
    dat$response <- scales::rescale(dat$response, 
                                    to = c(0.01,100), 
                                    from = range(dat$response), # add outlier reduction here?
                                    na.rm = TRUE, 
                                    finite = TRUE)
  }
  return(dat)
}

# extracts upper/lower estimates from drc model, single-curve
extract_upper_lower <- function(m_4param){
  est_a <- coef(m_4param)[2]
  est_b <- coef(m_4param)[3]
  if(est_b > est_a){
    return (list("c" = est_a, "d" = est_b))
  }else{return (list("c" = est_b, "d" = est_a))}
}

# batch correction lme4 auto script
batch_correct <- function(formula_,feature, data, r.eff, target_col, debug = FALSE){
  if (debug) browser()
  pheno_raw <- as.numeric(unlist(data[,c(feature)]))
  mn.pheno <- mean(pheno_raw, na.rm=TRUE)
  sd.pheno <- sd(pheno_raw, na.rm=TRUE)
  data$pheno <- (pheno_raw - mn.pheno) / sd.pheno
  fit <- lmer(formula_, data, REML=FALSE)
  random_effects <- ranef(fit)
  effects <- random_effects[[r.eff]]
  batch.ef <- effects[as.character(data[,c(target_col)][[1]]),]
  pheno_sub <- data$pheno - batch.ef
  pheno_adj <- (sd.pheno * pheno_sub) + mn.pheno
  return (pheno_adj)
}