#' Make a network object for contrast-level data containing data, priors, and a JAGS model file
#' 
#' This is similar to the function \code{\link{network.data}}, except it uses contrast-level data instead of arms-level data. Contrast-level format uses treatment differences relative to the control arm.
#' Note that in two arm trials there is only one contrast value per trial, but in three arm trials there are two contrast values relative to the control arm.
#'
#' @param Outcomes A vector of Contrast-level outcomes. Outcome is assumed to be normally distributed. If there are three arms in a trial, need to include two contrast values for that trial. See parkinsons_contrast data for an example.
#' @param Treat A vector of treatments for each arm. Treatments should have positive integer values starting from 1 to total number of treatments.
#' @param SE A vector of standard error for each contrasts. 
#' @param na A vector of number of arms in each study.
#' @param V Needed if you have multi-arm trials. Length of this vector should be number of studies. If the study is multi-arm trial, need to specify variance of the baseline treatment in that trial. Denote it with NA if the study only has two-arm trials.
#' @param type Type of model fitted: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @param rank.preference Set it equal to "higher" if higher values are preferred (i.e. assumes events are good). Set it equal to "lower" if lower values are preferred (i.e. assumes events are bad). Default is "higher".
#' @param mean.d Prior mean for the relative effect
#' @param prec.d Prior precision for the relative effect
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal. It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm, dwish) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter.
#' @references A.J. Franchini, S. Dias, A.E. Ades, J.P. Jansen, N.J. Welton (2012), \emph{Accounting for correlation in network meta-analysis with multi-arm trials}, Research Synthesis Methods 3(2):142-160. [\url{https://doi.org/10.1002/jrsm.1049}] 
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. [\url{https://doi.org/10.1177/0272989X12458724}]
#' @examples
#' network <- with(parkinsons_contrast, {
#'  contrast.network.data(Outcomes, Treat, SE, na, V)
#' })
#' @export

contrast.network.data <- function(Outcomes, Treat, SE, na, V = NULL, type = "random", rank.preference = "higher",
                                  mean.d = 0, prec.d = 0.0001, hy.prior = list("dunif", 0, 100)){
  
  if(missing(Outcomes) || missing(Treat) || missing(SE) || missing(na)){
    stop("Outcomes, Treat, SE, and na have to be all specified")
  }
  
  if(any(na == 1)) stop("study cannot have only 1 arm")
  if(is.unsorted(na)) stop("please sort the studies so that studies with higher number of arms are at the end (in a increasing order)")
  
  if(max(na) >= 3 & is.null(V)){
    stop("Need to specify variance of the baseline treatment if you have multi-arm trials")
  } else if(max(na) < 3 & !is.null(V)){
    stop("Since you do not have multi-arm trial, V is not needed")
  }
  
  if(!type %in% c("fixed", "random")){
    stop("type has to be either fixed or random")
  }
  
  if(!rank.preference %in% c("higher", "lower")){
    stop("rank preference has to be either higher or lower")
  }
  
  Outcomes <- contrast.make.matrix(Outcomes, na -1)
  SE <- contrast.make.matrix(SE, na -1)
  Treat <- contrast.make.matrix(Treat, na)
  
  # Attach NA column for the first column
  Outcomes <- cbind(NA, Outcomes)
  SE <- cbind(NA, SE)
  
  na_count <- as.vector(table(na))
  ntreat <- unique(as.vector(Treat))
  ntreat <- length(ntreat[!is.na(ntreat)])
  nstudy <- sum(na_count)
  
  network <- list(Outcomes = Outcomes, Treat = Treat, SE = SE, na = na, na_count = na_count, ntreat = ntreat, nstudy = nstudy, type = type, mean.d = mean.d, prec.d = prec.d, hy.prior = hy.prior, response = "normal", Treat.order = 1:ntreat, rank.preference = rank.preference, V = V)
  
  code <- contrast.network.rjags(network)
  network$code <- code
  
  class(network) <- "contrast.network.data"
  return(network)
  
}

contrast.network.rjags <- function(network){
  
  with(network, {
    
    code <- paste0("model\n{",
                   "\n\tfor(i in 1:", na_count[1], ") {",
                   "\n\t\ty[i,2] ~ dnorm(delta[i,2], prec[i,2])",
                   "\n\t\tresdev[i] <- (y[i,2] - delta[i,2]) * (y[i,2] - delta[i,2]) * prec[i,2]",
                   "\n\t}")
    
    if(length(na_count) > 1){
      
      for(i in 2:length(na_count)){
        
        code <- paste0(code, "\n\tfor(i in ", cumsum(na_count)[i-1] + 1, ":", cumsum(na_count)[i], ") {", 
                       "\n\t\tfor(k in 1:(na[i]-1)) {",
                       "\n\t\t\tfor(j in 1:(na[i]-1)) {",
                       "\n\t\t\t\tSigma", i-1, "[i,j,k] <- V[i]*(1-equals(j,k)) + Var[i,k+1] * equals(j,k)",
                       "\n\t\t\t}",
                       "\n\t\t}",
                       "\n\t\tOmega", i-1, "[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma", i-1, "[i,,])",
                       "\n\t\ty[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]], Omega", i-1, "[i, 1:(na[i]-1), 1:(na[i]-1)])",
                       "\n\t\tfor(k in 1:(na[i]-1)){",
                       "\n\t\t\tydiff[i,k] <- y[i,(k+1)] - delta[i,(k+1)]",
                       "\n\t\t\tz[i,k] <- inprod(Omega", i-1, "[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])",
                       "\n\t\t}",
                       "\n\t\tresdev[i] <- inprod(ydiff[i,1:(na[i]-1)], z[i, 1:(na[i]-1)])",
                       "\n\t}")
      }  
    }
    
    if(type == "random"){
      code <- paste0(code, "\n\tfor(i in 1:", nstudy, ") {",
                     "\n\t\tw[i,1] <- 0",
                     "\n\t\tdelta[i,1] <- 0",
                     "\n\t\tfor(k in 2:na[i]) {",
                     "\n\t\t\tVar[i,k] <- pow(se[i,k], 2)",
                     "\n\t\t\tprec[i,k] <- 1/Var[i,k]",
                     "\n\t\t}", 
                     "\n\t\tfor(k in 2:na[i]) {",
                     "\n\t\t\tdelta[i,k] ~ dnorm(md[i,k], taud[i,k])",
                     "\n\t\t\tmd[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k]",
                     "\n\t\t\ttaud[i,k] <- tau * 2 * (k-1)/k",
                     "\n\t\t\tw[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]])",
                     "\n\t\t\tsw[i,k] <- sum(w[i,1:(k-1)])/ (k-1)",
                     "\n\t\t}",
                     "\n\t}")
    } else if(type == "fixed"){
      code <- paste0(code, "\n\tfor(i in 1:", nstudy, ") {",
                     "\n\t\tfor(k in 2:na[i]) {",
                     "\n\t\t\tVar[i,k] <- pow(se[i,k],2)",
                     "\n\t\t\tprec[i,k] <- 1/Var[i,k]",
                     "\n\t\t\tdelta[i,k] <- d[t[i,k]] - d[t[i,1]]",
                     "\n\t\t}",
                     "\n\t}")
    }
    
    code <- paste0(code, "\n\ttotresdev <- sum(resdev[])",
                   "\n\td[1] <- 0",
                   "\n\tfor(k in 2:", ntreat, ") {",
                   "\n\t\td[k] ~ dnorm(", mean.d, ",", prec.d, ")",
                   "\n\t}")
    
    if(type == "random"){
      code <- paste0(code, contrast.hy.prior.rjags(hy.prior))
    }
    code <- paste0(code, rank.rjags(rank.preference, ntreat))
    
    code <- paste0(code, "\n}")
    return(code)
  })
}

#' Run the model using the network object
#' 
#' This is similar to the function \code{\link{network.run}}, except it uses contrast-level data instead of arms-level data.
#'
#' @param network contrast level network object created from \code{\link{contrast.network.data}} function
#' @param inits Initial values for the parameters being sampled. If left unspecified, program will generate reasonable initial values.
#' @param n.chains Number of chains to run
#' @param max.run Maximum number of iterations that user is willing to run. If the algorithm is not converging, it will run up to \code{max.run} iterations before printing a message that it did not converge
#' @param setsize Number of iterations that are run between convergence checks. If the algorithm converges fast, user wouldn't need a big setsize. The number that is printed between each convergence checks is the gelman-rubin diagnostics and we would want that to be below the conv.limit the user specifies.
#' @param n.run Final number of iterations that the user wants to store. If after the algorithm converges, user wants less number of iterations, we thin the sequence. If the user wants more iterations, we run extra iterations to reach the specified number of runs
#' @param conv.limit Convergence limit for Gelman and Rubin's convergence diagnostic. Point estimate is used to test convergence of parameters for study effect (eta), relative effect (d), and heterogeneity (log variance (logvar)).
#' @param extra.pars.save Parameters that user wants to save besides the default parameters saved. See code using \code{cat(network$code)} to see which parameters can be saved.
#' @return
#' \item{data_rjags}{Data that is put into rjags function jags.model}
#' \item{inits}{Initial values that are either specified by the user or generated as a default}
#' \item{pars.save}{Parameters that are saved. Add more parameters in extra.pars.save if other variables are desired}
#' \item{burnin}{Half of the converged sequence is thrown out as a burnin}
#' \item{n.thin}{If the number of iterations user wants (n.run) is less than the number of converged sequence after burnin, we thin the sequence and store the thinning interval}
#' \item{samples}{MCMC samples stored using jags. The returned samples have the form of mcmc.list and can be directly applied to coda functions}
#' \item{max.gelman}{Maximum Gelman and Rubin's convergence diagnostic calculated for the final sample}
#' \item{deviance}{Contains deviance statistics such as pD (effective number of parameters) and DIC (Deviance Information Criterion)}
#' \item{rank.tx}{Rank probability calculated for each treatments. \code{rank.preference} parameter in \code{\link{network.data}} is used to define whether higher or lower value is preferred. The numbers are probabilities that a given treatment has been in certain rank in the sequence.}
#' @examples
#' network <- with(parkinsons_contrast, {
#'  contrast.network.data(Outcomes, Treat, SE, na, V)
#' })
#' result <- contrast.network.run(network)
#' @export

contrast.network.run <- function(network, inits = NULL, n.chains = 3, max.run = 100000, setsize = 10000, n.run = 50000,
                                 conv.limit = 1.05, extra.pars.save = NULL){
  
  if (!inherits(network, "contrast.network.data")) {
    stop('Given network is not contrast.network.data. Run contrast.network.data function first')
  }
  
  if(max.run < setsize){
    stop("setsize should be smaller than max.run")
  }
  
  with(network, {
    
    data <- list(y = Outcomes, t = Treat, se = SE, na = na)
    
    if(!is.null(V)){
      data$V <- V
    }
    
    if(type == "random"){
      data$hy.prior.1 <- hy.prior[[2]]
      data$hy.prior.2 <- hy.prior[[3]]
    }
    
    pars.save <- c("d", "totresdev", "resdev", "prob")
    
    if(type == "random"){
      pars.save <- c(pars.save, "sd", "delta")  
    }
    
    if(!is.null(extra.pars.save)) {
      extra.pars.save.check(extra.pars.save, pars.save)
      pars.save <- c(pars.save, extra.pars.save)
    }
    
    if(is.null(inits)){
      inits <- contrast.inits(network, n.chains)
    }
    samples <- jags.fit(network, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)
    result <- list(network = network, data.rjags = data, inits = inits, pars.save = pars.save)
    result <- c(result, samples)
    
    result$deviance <- calculate.contrast.deviance(result)
    
    result$rank.tx <- rank.tx(result)
    class(result) <- "contrast.network.result"
    return(result)
  })
}

pick.summary.variables.contrast <- function(result, extra.pars = NULL, only.pars = NULL){
  samples <- result[["samples"]]
  varnames <- dimnames(samples[[1]])[[2]]
  varnames.split <- sapply(strsplit(varnames, "\\["), '[[', 1)
  varnames.split <- gsub("[[:digit:]]","",varnames.split)
  
  if(!is.null(only.pars)){
    if(!all(only.pars %in% varnames.split)){
      stop(paste0(only.pars, "was not sampled"))
    }
  }
  if(is.null(only.pars)){
    pars <- c("d", "sd")
  } else{
    pars <- only.pars
  }
  if(!is.null(extra.pars)){
    if(!extra.pars %in% varnames.split){
      stop(paste0(extra.pars, " is not saved in result"))
    }
    pars <- c(pars, extra.pars)
  }
  summary.samples <- lapply(samples, function(x){x[,varnames.split %in% pars, drop = F]})
  summary.samples <- coda::mcmc.list(summary.samples)
  summary.samples
}



#' Summarize result run by \code{\link{contrast.network.run}}
#'
#' This function uses summary function in coda package to summarize mcmc.list object. Monte carlo error (Time-series SE) is also obtained using the coda package and is printed in the summary as a default.
#'
#' @param object Result object created by \code{\link{contrast.network.run}} function
#' @param ... Additional arguments affecting the summary produced
#' @examples
#' network <- with(parkinsons_contrast, {
#'  contrast.network.data(Outcomes, Treat, SE, na, V)
#' })
#' result <- contrast.network.run(network) 
#' summary(result)
#' @export

summary.contrast.network.result <- function(object, ...){
  
  if(!inherits(object, "contrast.network.result")) {
    stop('This is not the output from contrast.network.run. Need to run contrast.network.run function first')
  }
  summary.samples <- pick.summary.variables.contrast(object, ...)
  
  rval <- list("summary.samples"= summary(summary.samples),
               "deviance" = unlist(object$deviance[1:3]),
               "total_n" = sum(object$network$na))
  class(rval) <- 'summary.contrast.network.result'
  rval
}


#' Plot traceplot and posterior density of the result using contrast data
#'
#' This function uses plotting function in coda package to plot mcmc.list object
#'
#' @param x Result object created by \code{\link{contrast.network.run}} function
#' @examples
#' network <- with(parkinsons_contrast, {
#'  contrast.network.data(Outcomes, Treat, SE, na, V)
#' })
#' result <- contrast.network.run(network)
#' plot(result)
#' @export

plot.contrast.network.result <- function(x) {
  
  if(!inherits(x, "contrast.network.result")) {
    stop('This is not the output from contrast.network.run. Need to run contrast.network.run function first')
  }
  plot(x$samples)
}


#' Find deviance statistics such as DIC and pD.
#'
#' Calculates deviance statistics. This function automatically called in \code{\link{contrast.network.run}} and the deviance statistics are stored after sampling is finished.
#'
#' @param result Object created by \code{\link{contrast.network.run}} function
#' @return
#' \item{Dbar}{Overall residual deviance}
#' \item{pD}{Sum of leverage_arm (i.e. total leverage)}
#' \item{DIC}{Deviance information criteria (sum of Dbar and pD)}
#' \item{resdev_study}{Posterior mean of the residual deviance in each study}
#' \item{devtilda_study}{Deviance at the posterior mean of the fitted values}
#' \item{leverage_study}{Difference between resdev_study and devtilda_study for each trial}
#' @examples
#' network <- with(parkinsons_contrast, {
#'  contrast.network.data(Outcomes, Treat, SE, na, V)
#' })
#' result <- contrast.network.run(network)
#' calculate.contrast.deviance(result)
#' @references A.J. Franchini, S. Dias, A.E. Ades, J.P. Jansen, N.J. Welton (2012), \emph{Accounting for correlation in network meta-analysis with multi-arm trials}, Research Synthesis Methods 3(2):142-160. [\url{https://doi.org/10.1002/jrsm.1049}] 
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. [\url{https://doi.org/10.1177/0272989X12458724}]
#' @export

calculate.contrast.deviance <- function(result){
  
  network <- result$network
  samples <- result$samples
  
  totresdev <- lapply(samples, function(x){ x[,"totresdev"]})
  Dbar <- mean(unlist(totresdev))
  
  #posterior mean of residual deviance 
  resdev <- lapply(samples, function(x) { x[,grep("resdev\\[", dimnames(samples[[1]])[[2]])]})
  resdev <- do.call(rbind, resdev)
  resdev_study <- apply(resdev, 2, mean)
  
  #devtilda - deviance at the posterior mean of the fitted values
  ybar <- lapply(samples, function(x){ x[,grep("delta\\[", dimnames(samples[[1]])[[2]])] })
  ybar <- do.call(rbind, ybar)
  ybar <- apply(ybar, 2, mean)
  ybar_study <- devtilda_study <- rep(NA, network$nstudy)
  
  with(network, {
    
    # 2 arm
    for(i in 1:na_count[1]){
      r_value <- Outcomes[i,2]
      se_value <- SE[i,2]
      ybar_study[i] <- ybar[which(paste("delta[", i, ",", 2, "]", sep = "") == names(ybar))]
      devtilda_study[i] <- ifelse(se_value != 0, (r_value - ybar_study[i])^2 / se_value^2, 0)
    }
    
    # 3 arm or more
    if(length(na_count) > 1){
      for(ii in 2:length(na_count)){
        for(i in (cumsum(na_count)[ii-1]+1): cumsum(na_count)[ii]){
          Sigma <- matrix(V[i], na[i] - 1, na[i] - 1)
          diag(Sigma) <- SE[i, 2:na[i]]
          omega_value <- solve(Sigma)
          r_value <- Outcomes[i,2:na[i]]
          ybar_arm <- ybar[grepl(paste0(i, ","), names(ybar), fixed=TRUE)]
          ybar_arm <- ybar_arm[!grepl(paste0(",", 1), names(ybar_arm), fixed = TRUE)] #get rid of the delta[,1] column if it exists
          devtilda_study[i] <- (r_value - ybar_arm) %*% omega_value %*% (r_value - ybar_arm)
        }   
      }  
    }
    leverage_study <- resdev_study - devtilda_study
    pD <- sum(leverage_study, na.rm = TRUE)
    DIC <- Dbar + pD
    
    return(list(Dbar = Dbar, pD = pD, DIC = DIC, resdev_study = resdev_study, devtilda_study = devtilda_study, leverage_study = leverage_study))
  })
}


#' Make a contrast network deviance plot
#'
#' This makes a contrasrt network deviance plot which plots residual deviance (resdev_study) vs. all study.
#' @param result Object created by \code{\link{contrast.network.run}} function
#' @examples
#' network <- with(parkinsons_contrast, {
#' contrast.network.data(Outcomes, Treat, SE, na, V)
#' })
#' result <- contrast.network.run(network)
#' contrast.network.deviance.plot(result)
#' @export

contrast.network.deviance.plot <- function(result){
  deviance <- result$deviance
  dev_vector <- deviance$resdev_study
  dev_vector <- dev_vector[!is.na(dev_vector)]
  plot(1:result$network$nstudy, dev_vector, xlab = "Study", ylab = "Residual Deviance", main = "Per-study residual deviance")
}

#' Make a leverage plot
#'
#' This function makes a leverage vs. square root of residual deviance plot
#'
#' @param result Object created by \code{\link{contrast.network.run}} function
#' @export

contrast.network.leverage.plot <- function(result){
  deviance <- result$deviance
  resdev <- sqrt(deviance$resdev_study)
  leverage <- deviance$leverage_study
  plot(resdev, leverage, xlim = c(0, max(c(resdev, 2.5))), ylim = c(0, max(c(leverage,4))),
       xlab = "Square root of residual deviance", ylab = "Leverage", main = "Leverage versus residual deviance")
}

contrast.inits <- function(network, n.chains){
  
  with(network, {
    delta <- as.vector(t(Outcomes))
    delta <- delta[!is.na(delta)]
    
    Treat <- as.vector(t(Treat))
    Treat <- Treat[!is.na(Treat)]
    
    ends <- cumsum(na) # End row of trials
    starts <- c(1, ends[-length(ends)] + 1) # Start row of trials
    b.Treat <- rep(NA, length(na))
    b.id <- rep(F, sum(na))
    for (i in 1:length(na)){
      limits <- starts[i]:ends[i]
      b.Treat[i] <- min(Treat[limits])
      b.id[limits[b.Treat[i] == Treat[limits]]] <- T
    }
    
    # design matrix
    base.tx <- Treat[b.id]    # base treatment for N studies
    end.Study <- c(0, cumsum(na))  # end row number of each trial
    rows <- end.Study - seq(0, nstudy)   # end number of each trial not including base treatment arms
    design.mat <- matrix(0, sum(na) - nstudy, ntreat) # no. non-base arms x #txs
    for (i in seq(nstudy)){
      studytx <- Treat[(end.Study[i]+1):end.Study[i+1]]  #treatments in ith Study
      nonbase.tx <- studytx[studytx!=base.tx[i]]    #non-baseline treatments for ith Study
      design.mat[(1+rows[i]):rows[i+1],base.tx[i]] <- -1
      for (j in seq(length(nonbase.tx)))
        design.mat[j+rows[i],nonbase.tx[j]] <- 1
    }
    design.mat <- design.mat[,-1,drop=F]
    
    fit <- summary(lm(delta ~ design.mat - 1))
    d <- se.d <- rep(NA, ntreat)
    d[-1] <- coef(fit)[,1]
    se.d[-1] <- coef(fit)[,2]
    resid.var <- fit$sigma^2
    
    initial.values = list()
    for(i in 1:n.chains){
      initial.values[[i]] = list()
    }
    
    if(!is.nan(fit$fstat[1])){
      for(i in 1:n.chains){
        random.d = rnorm(length(d))
        initial.values[[i]][["d"]] <- d + se.d * random.d
        
        if(type == "random"){
          
          df <- fit$df[2]
          random.ISigma <- rchisq(1, df)
          sigma2 <- resid.var * df/random.ISigma
          
          if(hy.prior[[1]] == "dunif"){
            if(sqrt(sigma2) > network$hy.prior[[3]]){
              stop("data has more variability than your prior does")
            }
          }
          
          if(hy.prior[[1]] == "dgamma"){
            initial.values[[i]][["tau"]] <- 1/sigma2
          } else if(hy.prior[[1]] == "dunif" || hy.prior[[1]] == "dhnorm"){
            initial.values[[i]][["sd"]] <- sqrt(sigma2)
          }
        }
      }
    }
    return(initial.values)
  })
}


contrast.hy.prior.rjags <- function(hy.prior){
  
  code <- ""
  distr <- hy.prior[[1]]
  if (distr == "dunif") {
    code <- paste0(code,
                   "\n\tsd ~ dunif(hy.prior.1, hy.prior.2)",
                   "\n\ttau <- pow(sd,-2)")
  } else if(distr == "dgamma"){
    code <- paste0(code,
                   "\n\tsd <- pow(tau, -0.5)",
                   "\n\ttau ~ dgamma(hy.prior.1, hy.prior.2)")
  } else if(distr == "dhnorm"){
    code <- paste0(code,
                   "\n\tsd ~ dnorm(hy.prior.1, hy.prior.2)T(0,)",
                   "\n\ttau <- pow(sd, -2)")
  }
  return(code)
}


contrast.make.matrix <- function(vec, na){
  
  nstudy <- length(na)
  mat <- matrix(NA, nstudy, max(na))
  Study <- rep(1:nstudy, na)
  
  arms_index <- NULL
  for(i in 1:length(na)){
    arms_index <- c(arms_index, seq(na[i]))
  }
  for(i in 1:length(vec)){
    mat[Study[i], arms_index[i]] <- vec[i]  
  }
  mat
}
