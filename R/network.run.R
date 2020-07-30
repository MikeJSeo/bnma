#' Run the model using the network object
#' 
#' This is the core function that runs the model in our program. Before running this function, we need to specify data, prior, JAGS code, etc. using \code{\link{network.data}}.
#'
#' @param network Network object created from \code{\link{network.data}} function
#' @param inits Initial values for the parameters being sampled. If left unspecified, program will generate reasonable initial values.
#' @param RNG.inits List of .RNG.name and .RNG.seed that control the JAGS RNGs. Please refer to jags.model function in rjags for more information.
#' @param n.chains Number of chains to run
#' @param max.run Maximum number of iterations that user is willing to run. If the algorithm is not converging, it will run up to \code{max.run} iterations before printing a message that it did not converge
#' @param setsize Number of iterations that are run between convergence checks. If the algorithm converges fast, user wouldn't need a big setsize. The number that is printed between each convergence checks is the gelman-rubin diagnostics and we would want that to be below the conv.limit the user specifies.
#' @param n.run Final number of iterations that the user wants to store. If after the algorithm converges, user wants less number of iterations, we thin the sequence. If the user wants more iterations, we run extra iterations to reach the specified number of runs
#' @param conv.limit Convergence limit for Gelman and Rubin's convergence diagnostic. Point estimate is used (instead of 95 percent C.I.) to test convergence of parameters for study effect (eta), relative effect (d), and heterogeneity (log variance (logvar)).
#' @param extra.pars.save Parameters that user wants to save besides the default parameters saved. See code using \code{cat(network$code)} to see which parameters can be saved.
#' @return
#' \item{data_rjags}{Data that is put into rjags function \code{\link{jags.model}}}
#' \item{inits}{Initial values that are either specified by the user or generated as a default}
#' \item{RNG.inits}{List of .RNG.name and .RNG.seed used for reproducibility}
#' \item{pars.save}{Parameters that are saved. Add more parameters in extra.pars.save if other variables are desired}
#' \item{burnin}{Half of the converged sequence is thrown out as a burnin}
#' \item{n.thin}{If the number of iterations user wants (n.run) is less than the number of converged sequence after burnin, we thin the sequence and store the thinning interval}
#' \item{samples}{MCMC samples stored using jags. The returned samples have the form of mcmc.list and can be directly applied to coda functions}
#' \item{max.gelman}{Maximum Gelman and Rubin's convergence diagnostic calculated for the final sample}
#' \item{deviance}{Contains deviance statistics such as pD (effective number of parameters) and DIC (Deviance Information Criterion)}
#' \item{rank.tx}{Rank probability calculated for each treatments. \code{rank.preference} parameter in \code{\link{network.data}} is used to define whether higher or lower value is preferred. The numbers are probabilities that a given treatment has been in certain rank in the sequence.}
#' @examples
#' #parkinson's example (normal)
#' network <- with(parkinsons,{
#'   network.data(Outcomes, Study, Treat, SE = SE, response = "normal")
#' })
#' \donttest{
#' result <- network.run(network)
#' }
#' @export

network.run <- function(network, inits = NULL, RNG.inits = NULL, n.chains = 3, max.run = 100000, setsize = 10000, n.run = 50000,
                        conv.limit = 1.05, extra.pars.save = NULL){
  
  if (!inherits(network, "network.data")) {
    stop('Given network is not network.data. Run network.data function first')
  }
  
  if(max.run < setsize){
    stop("setsize should be smaller than max.run")
  }
  
  with(network, {
    
    data <- if(response == "binomial"){
      list(na = na, t = t, r = r, n = n)
    } else if(response == "normal"){
      list(na = na, t = t, r = r, se = se)
    } else if(response == "multinomial"){
      list(na = na, t = t)
    }    
    
    if(response == "multinomial"){
      data$r <- r
      data$n <- n
    }
    
    data <- append(data, prior.data)
  
    # add covariate info
    if(!is.null(covariate)){
      for(i in seq(dim(covariate)[2])){
        data[[paste("mx",i, sep = "")]] = network[[paste("mx",i, sep = "")]]
        data[[paste("x",i, sep = "")]] = network[[paste("x",i, sep = "")]]
      }
    }
    
    # add baseline info
    if(baseline != "none"){
      data$mx_bl = mx_bl
    }
    
    # add placebo event rate when calculating RR, RD, NNT
    if(!is.null(mean.A) & !is.null(prec.A)){
      data$mean.A <- mean.A
      data$prec.A <- prec.A
    }
    
    ########## parameters to save in the model
    pars.save <-
      if(response == "binomial" || response == "normal"){
        c("Eta", "d", "sd", "logvar","prob","delta")
      } else if(response == "multinomial"){
        c("Eta", "d", "sigma", "sigma_transformed","prob","delta")
      }
    
    if(response == "normal" & !is.null(mean.A) & !is.null(prec.A)){
      pars.save <- c(pars.save, "T")
    }
    
    if(response == "binomial" & !is.null(mean.A) & !is.null(prec.A)){
      pars.save <- c(pars.save, "T", "RD", "RR", "NNT")
    }
    
    if(type == "fixed"){
      pars.save <- pars.save[!pars.save %in% c("sd", "sigma", "logvar", "sigma_transformed", "delta")]
    }
    
    if(!is.null(extra.pars.save)) {
      extra.pars.save.check(extra.pars.save, pars.save)
      pars.save <- c(pars.save, extra.pars.save)
    }
    
    pars.save <- c(pars.save, "totresdev")
    if(response == "binomial"){
      pars.save <- c(pars.save, "rhat", "dev")
    } else if(response == "normal"){
      pars.save <- c(pars.save, "theta", "dev")
    } else if(response == "multinomial"){
      pars.save <- c(pars.save, "rhat", "dev")
    }
    
    if(baseline != "none"){
      pars.save <- c(pars.save, "b_bl")
      if(baseline %in% c("common", "exchangeable")){
        pars.save <- c(pars.save, "B")
      }
      if(baseline == "exchangeable"){
        if(response == "multinomial"){
          pars.save <- c(pars.save, "sigmaB")
        } else{
          pars.save <- c(pars.save, "sdB")  
        }
      }
    }
    if(!is.null(covariate)){
      for(i in seq(dim(covariate)[2])){
        pars.save = c(pars.save, paste("beta",i,sep = ""))
      }
    }
    pars.save <- unique(pars.save)
    
    if(is.null(inits)){
      if(!any(is.na(network$data))){
        inits <- network.inits(network, n.chains)
      }
    }
    
    if(is.null(RNG.inits)){
      RNG.inits <- list()
      for(i in 1:n.chains){
        RNG.inits[[i]] <- list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = i)
      }
    }
    
    inits <- mapply(c, inits, RNG.inits, SIMPLIFY = FALSE)
    
    samples <- jags.fit(network, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)
    
    result <- list(network = network, data.rjags = data, inits = inits, RNG.inits= RNG.inits, pars.save = pars.save)
    result <- c(result, samples)
    
    result$deviance <- calculate.deviance(result)
    result$rank.tx <- rank.tx(result)
    
    class(result) <- "network.result"
    return(result)
  })
}


jags.fit <- function(network, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit) {
  
  mod = rjags::jags.model(textConnection(network$code), data = data, inits = inits, n.chains = n.chains, n.adapt = 0)
  
  adapted <- FALSE
  count <- 0
  while(!adapted){
    adapted <- rjags::adapt(mod, setsize, end.adaptation = FALSE)
    count <- count + 1
    if(count == 100){
      stop("algorithm has not adapted")
    }
  }
  
  if(class(network) == "network.data"){
    
    conv.save <- if(network$response == "multinomial"){
      c("d", "Eta", "sigma_transformed")
    } else if(network$response == "binomial" || network$response == "normal"){
      c("d", "Eta", "logvar")
    }
    if(network$type == "fixed"){
      conv.save <- conv.save[!conv.save %in% c("logvar", "sigma_transformed")]
    }  
  } else if(class(network) == "contrast.network.data" || class(network) == "ume.network.data"){
    conv.save <- pars.save
  } else if(class(network) == "nodesplit.network.data"){
    conv.save <- c("d", "sd", "diff")
  }
  
  
  samples <- rjags::coda.samples(model = mod, variable.names = pars.save, n.iter = setsize)
  varnames <- dimnames(samples[[1]])[[2]]
  varnames.split <- sapply(strsplit(varnames, "\\["), '[[', 1)
  conv.save.variables <- varnames.split %in% conv.save
  
  max.gelman <- find.max.gelman(samples, conv.save.variables)
  print(max.gelman)
  check <- max.gelman > conv.limit
  
  if(check) {
    count <- 1
    while (check & count < max.run/setsize) {
      samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = setsize)
      samples <- add.mcmc(samples, samples2)
      
      count <- count + 1
      
      max.gelman <- find.max.gelman(samples, conv.save.variables)
      check <- max.gelman > conv.limit
      print(max.gelman)
    }
  }
  
  start <- mcpar(samples[[1]])[1]
  end <- mcpar(samples[[1]])[2]
  mid <- (end + start-1)/2
  burnin <- ceiling(end - mid)
  samples <- window(samples, mid+1, end, 1) #keep the last half of the converged sequence
  samples <- new.mcmc(samples)
  
  n.thin <- 1
  if(check == TRUE){
    print("code didn't converge according to gelman-rubin diagnostics")
  } else if(n.run < burnin){
    n.thin <- ceiling(burnin/n.run)
    extra.run <- n.run * n.thin - burnin
    if(extra.run != 0){
      samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = extra.run)
      samples <- add.mcmc(samples, samples2)
    }
    samples <- window(samples, 1, dim(samples[[1]])[1], n.thin)
  } else if(n.run > burnin){
    extra.run <- n.run - burnin
    samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = extra.run)
    samples <- add.mcmc(samples, samples2)
  }
  max.gelman <- find.max.gelman(samples, conv.save.variables)
  print(max.gelman)
  
  out <-list(burnin = burnin, n.thin = n.thin, samples = samples, max.gelman = max.gelman)
  return(out)
}

extra.pars.save.check <- function(extra.pars.save, pars.save){
  if(!is.atomic(extra.pars.save) || !is.vector(extra.pars.save)) stop("extra pars should be a vector of strings")
  for(i in 1:length(extra.pars.save)){
    if(!is.character(extra.pars.save[i])) stop("extra pars should be a vector of strings")
    if(extra.pars.save[i] %in% pars.save) stop(paste0(extra.pars.save[i], " is already one of default parameters to save") )
  }
}

new.mcmc <- function(x){
  n.chains <- length(x)
  n.var <- nvar(x)
  newobjects <- vector("list", length = n.chains)
  
  for(i in 1:n.chains){
    newobjects[[i]] <- matrix(NA, nrow = 0, ncol = n.var, dimnames = list(NULL, dimnames(x[[1]])[[2]]))
    newobjects[[i]] <- x[[i]]
    newobjects[[i]] <- mcmc(newobjects[[i]])
  }
  mcmc.list(newobjects)
}

add.mcmc <- function(x, y){
  
  n.chains <- length(x)
  n.var <- nvar(x)
  newobjects <- vector("list", length = n.chains)
  
  for(i in 1:n.chains){
    newobjects[[i]] <- matrix(NA, nrow = 0, ncol = n.var, dimnames = list(NULL, dimnames(x[[1]])[[2]]))
    newobjects[[i]] <- rbind(x[[i]], y[[i]])
    newobjects[[i]] <- mcmc(newobjects[[i]])
  }
  mcmc.list(newobjects)
}

find.max.gelman <- function(samples, index){
  
  samples2 <- lapply(samples, function(x){ x[,index]})
  samples2 <- lapply(samples2, function(x) { x[,colSums(abs(x)) != 0] })
  
  max(gelman.diag(samples2, multivariate = FALSE)$psrf[,1]) #look at point estimate instead of 95% C.I.
}

find.max.gelman.variable <- function(samples, index){
  samples2 <- lapply(samples, function(x){ x[,index]})
  samples2 <- lapply(samples2, function(x) { x[,colSums(abs(x)) != 0] })
  
  names(which.max(gelman.diag(samples2, multivariate = FALSE)$psrf[,1])) # find the biggest
}