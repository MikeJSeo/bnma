
#' Make a network object containing data, priors, and a JAGS model file
#'
#' This function makes a network object that can be used to run network meta-analysis using \code{\link{nodesplit.network.run}}.
#' User needs to specify Outcomes, Study, Treat, N or SE, and response.
#' Prior parameters are filled in automatically based on the data type if not specified.
#' The input data should be arm-level so that we have observations for each treatment in each study.
#' The input data is preprocessed to fit the format necessary to run model in JAGS.
#'
#' @param Outcomes Arm-level outcomes. If it is a multinomial response, the matrix would have dimensions treatment arms (row) by multinomial categories (column). If it is binomial or normal, it would be a vector.
#' @param Study A vector of study indicator for each arm
#' @param Treat A vector of treatment indicator for each arm
#' @param N A vector of total number of observations in each arm. Used for binomial and multinomial responses
#' @param SE A vector of standard error for each arm. Used only for normal response
#' @param response Specification of the outcomes type. Must specify one of the following: "normal", "binomial", or "multinomial"
#' @param Treat.order Treatment order which determines how treatments are compared. The first treatment that is specified is considered to be the baseline treatment. Default order is alphabetical. If the treatments are coded 1, 2, etc, then the treatment with a value of 1 would be assigned as a baseline treatment.
#' @param pair Define a pair to split. It has to be a vector of length 2 with treatment names
#' @param type Type of model fitted: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @param dic This is an indicator for whether user wants to calculate DIC. Model stores less information if you set it to FALSE.
#' @return Creates list of variables that are used to run the model using \code{\link{nodesplit.network.run}}
#' \item{data}{Data combining all the input data. User can check this to insure the data is correctly specified. For modelling purposes, character valued studies or treatment variables are changed to numeric values based on alphabetical order.}
#' \item{nrow}{Total number of arms in the meta-analysis}
#' \item{ncat}{Number of columns in the Outcomes. Will equal 1 for binary and normal and number of categories for multinomial}
#' \item{nstudy}{Number of study}
#' \item{na}{Number of arms for each study}
#' \item{ntreat}{Number of treatment}
#' \item{b.id}{Indicator in sequence of all treatments for which treatment is base treatment in Study}
#' \item{t}{\code{Treat} transformed into a matrix which has dimensions number of study by max number of arms in studies}
#' \item{r}{\code{Outcomes} made into an array that is suitable for use in rjags code. For multinomial, it has 3 dimensions: number of study by max number of arms in studies by number of categories.}
#' \item{code}{Rjags model file code that is generated using information provided by the user. To view model file inside R in a nice format, use \code{cat(network$code).}}
#' @examples
#' ###Thrombolytic data example
#' network <- with(thrombolytic,{
#'  nodesplit.network.data(Outcomes, Study, Treat, N, response = "binomial", pair = c(3,9)))
#' })
#' network
#' @references S. Dias, N.J. Welton, D.M. Caldwellb, A.E. Ades (2010), \emph{Checking consistency in mixed treatment}, Statistics in Medicine 29(7-8, Sp. Iss. SI): 932-944. [\url{https://doi.org/10.1002/sim.3767}]
#' @export

nodesplit.network.data <- function(Outcomes = NULL, Study = NULL, Treat = NULL, N = NULL, SE = NULL, response = NULL,  Treat.order = NULL,
                                  pair = NULL, type = "random", dic = TRUE){
  
  if(response == "multinomial"){
    stop("Not yet implemented")
  }
  
  if(is.null(pair)){
    stop("Pair has to be specified")
  }
  
  if(!all(pair %in% Treat)){
    stop("Pair name has to be exactly same as one of treatment names specified")
  }
  
  # rename the variables and order them based on specified treatment order
  network <- preprocess.data(Outcomes = Outcomes, Study = Study, Treat = Treat, N = N, SE = SE, response = response, Treat.order = Treat.order)
  
  # find characteristic values associated with the network, i.e. number of studies, number of treatments, etc
  characteristics <- find.characteristics(network)
  network <- append(network, characteristics)
  
  # change input dimensions for Outcomes, Study, Treat, and N/SE to fit JAGS coding format
  new.inputs <- change.dimensions(network)
  network <- append(network, new.inputs)
  
  # Nodesplit parameters
  pair_numeric <- c(which(network$Treat.order == pair[1]), which(network$Treat.order == pair[2]))
  pair_numeric <- pair_numeric[order(pair_numeric)]
  
  network <- append(network, list(pair = pair_numeric))
  nodesplit.inputs <- find.nodesplit.parameters(network)
  network <- append(network, nodesplit.inputs)
  
  # generate JAGS code
  network$type <- type
  code <- nodesplit.network.rjags(network)
  network$code <- code
  network$dic <- dic
  
  class(network) <- "nodesplit.network.data"
  return(network)
}


#' Run the model using the nodesplit network object
#' 
#' This is similar to the function \code{\link{network.run}}, except this is used for the nodesplitting model.
#'
#' @param network network object created from \code{\link{nodesplit.network.data}} function
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
#' @examples
#' ###Thrombolytic data example
#' network <- with(thrombolytic,{
#'  nodesplit.network.data(Outcomes, Study, Treat, N, response = "binomial", pair = c(3,9)))
#' })
#' result <- nodesplit.network.run(network)
#' @export

nodesplit.network.run <- function(network, inits = NULL, n.chains = 3, max.run = 100000, setsize = 10000, n.run = 50000,
                            conv.limit = 1.05, extra.pars.save = NULL){
  
  if (!inherits(network, "nodesplit.network.data")) {
    stop('Given network is not nodesplit.network.data. Run nodesplit.network.data function first')
  }
  
  if(max.run < setsize){
    stop("setsize should be smaller than max.run")
  }
  
  with(network, {
    
    data <- list(r = r, t = t, na = na)
    
    if(response == "binomial" || response == "multinomial"){
      data$n <- n
    } else if(response == "normal"){
      data$se <- se
    }
    data <- append(data, list(pair = pair, split = split, m =  m, bi = bi, si = si))
    
    pars.save <- c("direct","d", "sd", "diff", "prob", "oneminusprob")
    
    if(dic == TRUE){
    pars.save <- c(pars.save, "totresdev")
    if(response == "binomial"){
      pars.save <- c(pars.save, "rhat", "dev")
    } else if(response == "normal"){
      pars.save <- c(pars.save, "theta", "dev")
    }
    # else if(response == "multinomial"){
    #   pars.save <- c(pars.save, "rhat", "dev")
    # }
    }
    
    if(!is.null(extra.pars.save)) {
      extra.pars.save.check(extra.pars.save, pars.save)
      pars.save <- c(pars.save, extra.pars.save)
    }
    
    if(is.null(inits)){
      inits <- list()
      for(i in 1:n.chains){
        inits[[i]] <- list(direct=0,  d= c(NA, rep(0, ntreat - 1)), sd=1, mu=rep(0,nstudy))  
      }
    }

    samples <- jags.fit(network, data, pars.save, inits = NULL, n.chains, max.run, setsize, n.run, conv.limit)
    
    result <- list(network = network, data.rjags = data, inits = inits, pars.save = pars.save)
    result <- c(result, samples)
    
    result$deviance <- calculate.deviance(result)

    class(result) <- "nodesplit.network.result"
    return(result)
  })
}

nodesplit.network.rjags <- function(network){
  
  response <- network$response
  
  code <- paste0("model{\n")
  code2 <- if(response == "binomial"){
    nodesplit.model.binomial(network)
  } else if(response == "normal"){
    nodesplit.model.normal(network)
  }
  #} else if(response == "multinomial"){
  #  nodesplit.model.multinomial(network)
  #}
  code <- paste0(code, code2, "\n}")
  return(code)
}


nodesplit.model.normal <- function(network){
  
  with(network, {
    
    code <- paste0("for(i in 1:", nstudy, "){",
                   "\n\tw[i,1] <- 0",
                   "\n\tj[i,1] <- 0",
                   "\n\tdelta[i, bi[i]] <- 0",
                   "\n\tmu[i] ~ dnorm(0, .0001)",
                   "\n\tfor(k in 1:na[i]) {",
                   "\n\t\tprec[i,k] <- 1/pow(se[i,k],2)",
                   "\n\t\tr[i,k] ~ dnorm(theta[i,t[i,k]], prec[i,k])", 
                   "\n\t\ttheta[i,t[i,k]] <- mu[i] + delta[i,t[i,k]]",
                   "\n\t\tindex[i,k] <- split[i] * (equals(t[i,k], pair[1]) + equals(t[i,k], pair[2]))",
                   "\n\t}")
    
    if(type == "random"){
      code <- paste0(code, 
                     "\n\tfor(k in 2:na[i]) {",
                     "\n\t\tdelta[i, si[i,k]] ~ dnorm(md[i, si[i,k]], taud[i, si[i,k]])",
                     "\n\t\tmd[i, si[i,k]] <- (d[si[i,k]] - d[bi[i]] + sw[i,k]) * (1 - index[i,m[i,k]]) + direct * index[i,m[i,k]]",
                     "\n\t\tj[i,k] <- k - (equals(1, split[i]) * step(k-3))",
                     "\n\t\ttaud[i, si[i,k]] <- tau * 2 * (j[i,k] -1) / j[i,k]",
                     "\n\t\tw[i,k] <- (delta[i, si[i,k]] - d[si[i,k]] + d[bi[i]]) * (1 - index[i,k])",
                     "\n\t\tsw[i,k] <- sum(w[i,1:(k-1)])/(j[i,k]-1)",
                     "\n\t}",
                     "\n}")
    } else if(type == "fixed"){
      
      code <- paste0(code, 
                     "\n\tfor(k in 2:na[i]) {",
                     "\n\t\tdelta[i,si[i,k]] <-(d[si[i,k]] - d[bi[i]] )*(1-index[i,m[i,k]]) + direct*index[i,m[i,k]]",
                     "\n\t}",
                     "\n}")
    }
    
    if(type == "random"){
      code <- paste0(code, "\nsd ~ dunif(0, 100)",
                     "\nvarr <- pow(sd,2)",
                     "\ntau <- 1/varr")
    }      
    
    code <- paste0(code,
                   "\nd[1] <- 0",
                   "\ndirect ~ dnorm(0, .0001)",
                   "\nfor(k in 2:", ntreat, ") { d[k] ~ dnorm(0, 0.0001) }",
                   "\nfor(c in 1:", ntreat -1, "){",
                   "\n\tfor(k in (c+1):", ntreat, ") {",
                   "\n\t\tlor[c,k] <- d[k] - d[c]",
                   "\n\t}",
                   "\n}",
                   "\ndiff <- direct - lor[pair[1], pair[2]]",
                   "\nprob <- step(diff)",
                   "\noneminusprob <- 1 - prob")
                   
    return(code)
  })
}

nodesplit.model.binomial <- function(network){
  
  with(network, {
    
  code <- paste0("for(i in 1:", nstudy, "){",
                 "\n\tw[i,1] <- 0",
                 "\n\tj[i,1] <- 0",
                 "\n\tdelta[i, bi[i]] <- 0",
                 "\n\tmu[i] ~ dnorm(0, .0001)",
                 "\n\tfor(k in 1:na[i]) {",
                 "\n\t\tr[i,k] ~ dbin(p[i,t[i,k]], n[i,k])",
                 "\n\t\tlogit(p[i,t[i,k]]) <- mu[i] + delta[i, t[i,k]]",
                 "\n\t\tindex[i,k] <- split[i] * (equals(t[i,k], pair[1]) + equals(t[i,k], pair[2]))",
                 "\n\t\trhat[i,k] <- p[i,t[i,k]] * n[i,k]",
                 "\n\t\tdev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k]))
                                  + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))",
                 "\n\t}",
                 "\n\tresdev[i]<-sum(dev[i,1:na[i]])"
                 )
  
  if(type == "random"){
    code <- paste0(code, 
                   "\n\tfor(k in 2:na[i]) {",
                   "\n\t\tdelta[i, si[i,k]] ~ dnorm(md[i, si[i,k]], taud[i, si[i,k]])",
                   "\n\t\tmd[i, si[i,k]] <- (d[si[i,k]] - d[bi[i]] + sw[i,k]) * (1 - index[i,m[i,k]]) + direct * index[i,m[i,k]]",
                   "\n\t\tj[i,k] <- k - (equals(1, split[i]) * step(k-3))",
                   "\n\t\ttaud[i, si[i,k]] <- tau * 2 * (j[i,k] -1) / j[i,k]",
                   "\n\t\tw[i,k] <- (delta[i, si[i,k]] - d[si[i,k]] + d[bi[i]]) * (1 - index[i,k])",
                   "\n\t\tsw[i,k] <- sum(w[i,1:(k-1)])/(j[i,k]-1)",
                   "\n\t}",
                   "\n}"
                   )
  } else if(type == "fixed"){
    code <- paste0(code, 
                   "\n\tfor(k in 2:na[i]) {",
                   "\n\t\tdelta[i,si[i,k]] <-(d[si[i,k]] - d[bi[i]] )*(1-index[i,m[i,k]]) + direct*index[i,m[i,k]]",
                   "\n\t}",
                   "\n}")
  }
  
  if(type == "random"){
    code <- paste0(code, "\nsd ~ dunif(0, 100)",
                   "\nvarr <- pow(sd,2)",
                   "\ntau <- 1/varr")
  }      
  
  code <- paste0(code,          
                 "\nd[1] <- 0",
                 "\ndirect ~ dnorm(0, .0001)",
                 "\nfor(k in 2:", ntreat, ") { d[k] ~ dnorm(0, 0.0001) }",
                 "\ntotresdev<-sum(resdev[])",
                 "\nfor(c in 1:", ntreat -1, "){",
                 "\n\tfor(k in (c+1):", ntreat, ") {",
                 "\n\t\tlor[c,k] <- d[k] - d[c]",
                 "\n\t}",
                 "\n}",
                 "\ndiff <- direct - lor[pair[1], pair[2]]",
                 "\nprob <- step(diff)",
                 "\noneminusprob <- 1 - prob"
                 )
  return(code)
  })
}


find.nodesplit.parameters <- function(network){
  
  with(network,{
    # Calculate split (1 if node to split is present) and b (baseline position)
    checkPair <- PairXY(t, pair)
    
    # Build vector bi[i] with baseline treatment: t[i, b[i]]
    bi <- Basetreat(t, checkPair[,"b"])
    
    # Indexes to sweep non-baseline arms only
    m <- NonbaseSweep(checkPair, na)
    
    # Build matrix si[i,k] with non-baseline treatments: t[i, m[i,k]]
    si <- Sweeptreat(t,m)
    
    return(list(split = checkPair[,"split"], m = m, bi = bi, si = si))
  })
}                             



PairXY <- function(treat, pair)
  # From appendix of Dias et al. 2010
  # Check if pair(X,Y) in row i of data 
  # and give baseline for data row i
{
  N <- nrow(treat)
  out <- cbind(split=rep(0,N), b=rep(0,N))
  for (i in 1:N) {
    # returns positions of matches to elements of pair in t[i,]
    # or zero if not present
    pos <- match(pair, treat[i,], nomatch=0)   # lenght = length(pair) = 2
    out[i,1] <- ifelse(prod(pos)>0, 1, 0)      # 1 if pair in line i, 0 o.w.
    out[i,2] <- ifelse(prod(pos)==0, 1, pos[1])
  }
  out
}

NonbaseSweep <- function(index, na)
  # From appendix of Dias et al. 2010
  # gives na-1 indexes to sweep non-baseline arms only
{
  N <- NROW(na)
  C <- max(na)
  out <- matrix(nrow=N, ncol=C)
  for (i in 1:N) {
    for (k in 2:na[i]) {
      out[i,k] <- k - (index[i,"b"] >= k)
    }
  }
  out
}
#
Sweeptreat <- function(treat, m)
  # From appendix of Dias et al. 2010
  # Builds matrix with non-baseline treatments
{
  N <- NROW(treat)
  C <- NCOL(m)
  out <- matrix(nrow=N, ncol=C)
  for (i in 1:N) {
    for (k in 2:C) {
      out[i,k] <- treat[i,m[i,k]]
    }
  }
  out
}

#
Basetreat <- function(treat, b)
  # From appendix of Dias et al. 2010
  # Builds vector with baseline treatments
{
  N <- nrow(treat)
  out <- rep(0,N)
  for (i in 1:N) {
    out[i] <- treat[i,b[i]]
  }
  out
}


pick.summary.variables.nodesplit <- function(result, extra.pars = NULL, only.pars = NULL){
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
    pars <- c("direct","d", "sd", "diff", "prob", "oneminusprob")
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



#' Summarize result run by \code{\link{nodesplit.network.run}}
#'
#' This function uses summary function in coda package to summarize mcmc.list object. Monte carlo error (Time-series SE) is also obtained using the coda package and is printed in the summary as a default.
#'
#' @param object Result object created by \code{\link{nodesplit.network.run}} function
#' @param ... Additional arguments affecting the summary produced
#' @examples
#' ###Thrombolytic data example
#' network <- with(thrombolytic,{
#'  nodesplit.network.data(Outcomes, Study, Treat, N, response = "binomial", pair = c(3,9)))
#' })
#' result <- nodesplit.network.run(network)
#' summary(result)
#' @export

summary.nodesplit.network.result <- function(object, ...){
  
  
  if(!inherits(object, "nodesplit.network.result")) {
    stop('This is not the output from nodesplit.network.run. Need to run nodesplit.network.run function first')
  }
  
  summary.samples <- pick.summary.variables.nodesplit(object, ...)
  
  rval <- list("summary.samples" = summary(summary.samples),
               "deviance" = unlist(object$deviance[1:3]),
               "Inconsistency estimate" = summary(object$samples)$statistics["diff","Mean"], 
               "p_value" = min(summary(object$samples)$statistics["prob", "Mean"], summary(object$samples)$statistics["oneminusprob", "Mean"])
               )
  class(rval) <- 'summary.nodesplit.network.result'
  rval
}
