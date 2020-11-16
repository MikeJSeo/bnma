network.inits <- function(network, n.chains){
  
  response <- network$response
  
  inits <- if(response == "multinomial"){
    multinomial.inits(network, n.chains)
  } else if(response == "binomial"){
    binomial.inits(network, n.chains)
  } else if(response == "normal"){
    normal.inits(network, n.chains)
  }
  return(inits)
}

normal.inits <- function(network, n.chains){
  
  with(network,{
    Eta <- Outcomes[b.id]
    se.Eta <- SE[b.id]
    delta <- Outcomes - rep(Eta, times = na)
    delta <- delta[!b.id,] #eliminate base-arm
    
    inits <- make.inits(network, n.chains, delta, Eta, se.Eta)
    return(inits)
  })
}

binomial.inits <- function(network, n.chains){
  
  with(network,{
    
    Outcomes <- Outcomes + 0.5 # ensure ratios are always defined
    N <- N + 1
    p <- Outcomes/N
    logits <- log(p/(1-p))
    se.logits <- sqrt(1/Outcomes + 1/(N - Outcomes))
    
    Eta <- logits[b.id]
    se.Eta <- se.logits[b.id]
    delta <- logits - rep(Eta, times = na)
    delta <- delta[!b.id,]
    
    inits = make.inits(network, n.chains, delta, Eta, se.Eta)
    return(inits)  
  })
}

make.inits <- function(network, n.chains, delta, Eta, se.Eta){
  
  with(network,{
    
    # dependent variable for regression
    y <- delta
    
    # design matrix
    base.tx <- Treat[b.id]    # base treatment for N studies
    end.Study <- c(0, cumsum(na))  # end row number of each trial
    rows <- end.Study - seq(0, nstudy)   # end number of each trial not including base treatment arms
    design.mat <- matrix(0, sum(na) - nstudy, ntreat) # no. non-base arms x #txs
    for (i in seq(nstudy)){
      studytx <- Treat[(end.Study[i]+1):end.Study[i+1]]  #treatments in ith Study
      nonbase.tx <- studytx[studytx!=base.tx[i]]    #non-baseline treatments for ith Study
      design.mat[(1+rows[i]):rows[i+1],base.tx[i]] <- -1
      for (j in seq(length(nonbase.tx))){
        design.mat[j+rows[i],nonbase.tx[j]] <- 1
      }
    }
    design.mat <- design.mat[,-1,drop=F]
    
    fit <- summary(lm(y ~ design.mat - 1))
    d <- se.d <- rep(NA, ntreat)
    
    # in case regression fails due to lack of data
    if(length(coef(fit)[,1]) == (ntreat -1)){
      d[-1] <- coef(fit)[,1]
      se.d[-1] <- coef(fit)[,2]  
    } else{
      d[-1] <- rep(0, ntreat -1)
      se.d[-1] <- rep(0.1, ntreat -1) 
    }
        
    resid.var <- fit$sigma^2
    
    # covariate
    if(!is.null(covariate)) {
      x.cen = matrix(0, nrow = sum(na), ncol = dim(covariate)[2])
      for(i in 1:dim(covariate)[2]){
        x.cen[,i] <- rep(covariate[,i], times = na)
      }
      x.cen <- x.cen[-seq(dim(x.cen)[1])[b.id],,drop=F]
      x.cen <- scale(x.cen, scale = FALSE)
      
      slope <- se.slope <- array(NA, c(ntreat, dim(covariate)[2]))
      for(i in 1:dim(covariate)[2]){
        fit2 <- if(covariate.model == "common" || covariate.model == "exchangeable"){
          summary(lm(y ~ x.cen[,i] -1))
        } else if(covariate.model == "independent"){
          summary(lm(y ~ design.mat:x.cen[,i] - 1))
        }
        
        if(length(coef(fit2)[,1]) == (ntreat-1) & covariate.model=="independent"){
          slope[-1,i] <- coef(fit2)[,1]
          se.slope[-1,i] <- coef(fit2)[,2]  
        } else if(length(coef(fit2)[,1]) == 1 & covariate.model %in% c("common", "exchangeable")){
          slope[-1,i] <- rep(coef(fit2)[,1], ntreat-1)
          se.slope[-1,i] <- rep(coef(fit2)[,2], ntreat-1)  
        } else{
          slope[-1,i] <- rep(0, ntreat-1)
          se.slope[-1,i] <- rep(0.1, ntreat-1)  
        }
      }
    }
    
    # baseline
    if(baseline != "none"){
      baseline.cen <- rep(Eta, na)
      baseline.cen <- baseline.cen[-seq(length(baseline.cen))[b.id]]
      baseline.cen <- scale(baseline.cen, scale = FALSE)
      
      baseline.slope <- baseline.se.slope <- rep(NA, ntreat)
      fit3 <- if(baseline == "common" || baseline == "exchangeable"){
        summary(lm(y ~ baseline.cen -1))
      } else if(baseline == "independent"){
        summary(lm(y ~ design.mat:baseline.cen - 1))
      }
      
      if(length(coef(fit3)[,1]) == (ntreat-1) & baseline=="independent"){
        baseline.slope[-1] <- coef(fit3)[,1]
        baseline.se.slope[-1] <- coef(fit3)[,2] 
      } else if(length(coef(fit3)[,1]) == 1 & baseline %in% c("common", "exchangeable")){
        baseline.slope[-1] <- rep(coef(fit3)[,1], ntreat-1)
        baseline.se.slope[-1] <- rep(coef(fit3)[,2], ntreat-1) 
      } else{
        baseline.slope[-1] <- rep(0, ntreat - 1)
        baseline.se.slope[-1] <- rep(0.1, ntreat - 1)
      }
    }
    
    ############## Generate initial values
    initial.values = list()
    for(i in 1:n.chains){
      initial.values[[i]] = list()
    }
    for(i in 1:n.chains){
      random.Eta <- rnorm(length(Eta))
      initial.values[[i]][["Eta"]] <- Eta + se.Eta * random.Eta
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
            if(sqrt(sigma2) > network$prior.data$hy.prior.2){
              stop("data has more variability than your prior does")
            }
          }
          
          if(hy.prior[[1]] == "dgamma"){
            initial.values[[i]][["prec"]] <- 1/sigma2
          } else if(hy.prior[[1]] == "dunif" || hy.prior[[1]] == "dhnorm"){
            initial.values[[i]][["sd"]] <- sqrt(sigma2)
          }
          
          # generate values for delta
          delta = matrix(NA, nrow = nrow(t), ncol = ncol(t))
          for(j in 2:ncol(delta)){
            diff_d <- ifelse(is.na(d[t[,1]]), d[t[,j]], d[t[,j]] - d[t[,1]])
            for(ii in 1:nrow(delta)){
              if(!is.na(diff_d[ii])) delta[ii,j] = rnorm(1, mean = diff_d[ii], sd = sqrt(sigma2))
            }
          }
          initial.values[[i]][["delta"]] <- delta
        }
      }
    }
    
    if (!is.null(covariate)) {
      if(!is.nan(fit2$fstat[1])){
        for(i in 1:n.chains){
          random.slope <- matrix(rnorm(dim(slope)[1]*dim(slope)[2]),dim(slope))
          for(j in 1:dim(covariate)[2]){
            initial.values[[i]][[paste("beta", j, sep = "")]] = slope[,j] + se.slope[,j] * random.slope[,j]
          }
        }
      }
    }
    
    if(baseline != "none"){
      if(!is.nan(fit3$fstat[1])){
        for(i in 1:n.chains){
          random.baseline = rnorm(length(baseline.slope))
          initial.values[[i]][["b_bl"]] = baseline.slope + baseline.se.slope * random.baseline
        }
      }
    }
    return(initial.values)
  })
  
  
}

############################################ multinomial inits functions

multinomial.inits <- function(network, n.chains)
{
  with(network,{
    
    Dimputed = Outcomes + 0.5
    
    logits <- as.matrix(log(Dimputed[, -1]) - log(Dimputed[, 1]))
    se.logits <- as.matrix(sqrt(1/Dimputed[, -1] + 1/Dimputed[, 1]))
    
    Eta <- se.Eta <- matrix(NA, nstudy, ncat)
    Eta[,2:ncat] <- logits[b.id,]
    se.Eta[,2:ncat] <- se.logits[b.id,]
    
    delta <- logits - apply(as.matrix(Eta[, -1]), 2, rep, times = na)
    rows.of.basetreat <- seq(dim(as.matrix(delta))[1])*as.numeric(b.id)
    delta <- delta[-rows.of.basetreat,,drop=F]   # Eliminate base treatment arms
    
    ###################### Using delta, Eta, and se.Eta make initial values
    
    y <- delta            # dependent variable for regression (part of Delta)
    d <- se.d <- matrix(NA, length(unique(Treat)), ncat - 1)
    resid.var <- rep(NA, ncat -1)
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
    
    for(k in 1:(ncat - 1)){
      fit <- summary(lm(y[,k] ~ design.mat - 1))
      
      if(length( coef(fit)[1:(ntreat-1), 1]) == (ntreat-1)){
        d[-1,k] <- coef(fit)[1:(ntreat-1), 1]
        se.d[-1,k] <- coef(fit)[1:(ntreat-1), 2]
      } else{
        d[-1,k] <- rep(0, ntreat - 1)
        se.d[-1,k] <- rep(0.1, ntreat -1)
      }
      resid.var[k] <- fit$sigma^2
    }
    
    # covariate
    if(!is.null(covariate)){
      x.cen <- matrix(0, nrow = sum(na), ncol = dim(covariate)[2])
      for(i in 1:dim(covariate)[2]){
        x.cen[,i] <- rep(covariate[,i], na)
      }
      x.cen <- x.cen[-seq(dim(x.cen)[1])[b.id],,drop=F]
      x.cen <- scale(x.cen, scale = FALSE)
      slope <- se.slope <- array(NA, c(ntreat, dim(covariate)[2], ncat - 1))
      
      for(i in 1:dim(covariate)[2]){
        for(k in 1:(ncat-1)){
          
          fit2 <- if(covariate.model == "independent" || covariate.model == "exchangeable"){
            summary(lm(y[,k] ~ design.mat:x.cen[,i] - 1))
          } else if(covariate.model == "common"){
            summary(lm(y[,k] ~ x.cen[,i] - 1))
          }
          
          if(length(coef(fit2)[,1]) == (ntreat-1) & covariate.model=="independent"){
            slope[-1,i,k] <- coef(fit2)[,1]
            se.slope[-1,i,k] <- coef(fit2)[,2]  
          } else if(length(coef(fit2)[,1]) == 1 & covariate.model %in% c("common", "exchangeable")){
            slope[-1,i,k] <- rep(coef(fit2)[,1], ntreat -1)
            se.slope[-1,i,k] <- rep(coef(fit2)[,2], ntreat -1)
          } else{
            slope[-1,i,k] <- rep(0, ntreat -1)
            se.slope[-1,i,k] <- rep(0.1, ntreat -1)  
          }
        }
      }
    }
    
    # baseline
    if(baseline != "none"){
      baseline.cen <- apply(as.matrix(Eta[, -1]), 2, rep, times = na)
      baseline.cen <- baseline.cen[-seq(dim(baseline.cen)[1])[b.id],]
      baseline.cen <- scale(baseline.cen, scale = FALSE)
      
      baseline.slope <- baseline.se.slope <- matrix(nrow = ntreat, ncol = ncat -1)
      
      for(k in 1:(ncat -1)){
        fit3 <- if(baseline == "common" || baseline == "exchangeable"){
          summary(lm(y[,k] ~ baseline.cen[,k] -1))
        } else if(baseline == "independent"){
          summary(lm(y[,k] ~ design.mat:baseline.cen[,k] - 1))
        }
        
        if(length(coef(fit3)[,1]) == (ntreat-1) & baseline == "independent"){
          baseline.slope[-1, k] <- coef(fit3)[,1]
          baseline.se.slope[-1, k] <- coef(fit3)[,2]
        } else if(length(coef(fit3)[,1]) ==1 & baseline %in% c("common", "exchangeable")){
          baseline.slope[-1, k] <- rep(coef(fit3)[,1], ntreat-1)
          baseline.se.slope[-1, k] <- rep(coef(fit3)[,2], ntreat-1)
        } else{
          baseline.slope[-1, k] <- rep(0, ntreat - 1)
          baseline.se.slope[-1, k] <- rep(0.1, ntreat - 1)
        }
      
      }
    }
    
    ################################################
    initial.values = list()
    for(i in 1:n.chains){
      initial.values[[i]] = list()
    }
    
    for(i in 1:n.chains){
      random.Eta <- matrix(rnorm(dim(Eta)[1]*dim(Eta)[2]),dim(Eta)[1],dim(Eta)[2])
      initial.values[[i]][["Eta"]] <- Eta + se.Eta * random.Eta
    }
    
    if(!is.nan(fit$fstat[1])){
      for(i in 1:n.chains){
        random.d = matrix(rnorm(dim(d)[1]*dim(d)[2]),dim(d)[1],dim(d)[2])
        initial.values[[i]][["d"]] = d + se.d * random.d
        
        if(type == "random"){
          df <- fit$df[2]
          random.ISigma <- rchisq(1, df)
          sigma2 <- resid.var * df/random.ISigma
          
          initial.values[[i]][["prec"]] <- 1/sigma2 * diag(ncat - 1)
          
          if(max(na) == 2){
            delta <- array(NA, dim = c(nstudy, max(na), ncat))
            for(j in 2:max(na)){
              for(m in 1:(ncat-1)){
                diff_d <- ifelse(is.na(d[t[,1],m]), d[t[,j],m], d[t[,j],m] - d[t[,1],m])
                for(ii in 1:nstudy){
                  if(!is.na(diff_d[ii])) delta[ii,j,m+1] <- rnorm(1, mean = diff_d[ii], sd = sqrt(sigma2))
                }
              }
            }
            initial.values[[i]][["delta"]] <- delta
          }
        }
      }
    }
    
    if (!is.null(covariate)) {
      if(!is.nan(fit2$fstat[1])){
        for(i in 1:n.chains){
          random.slope <- array(rnorm(dim(slope)[1]*dim(slope)[2]*dim(slope)[3]),dim(slope))
          for(j in 1:dim(covariate)[2]){
            initial.values[[i]][[paste("beta", j, sep = "")]] = slope[,j,] + se.slope[,j,] * random.slope[,j,]
          }
        }
      }
    }
    
    if(baseline != "none"){
      if(!is.nan(fit3$fstat[1])){
        for(i in 1:n.chains){
          random.baseline = matrix(rnorm(dim(baseline.slope)[1]*dim(baseline.slope)[2]),dim(baseline.slope))
          initial.values[[i]][["b_bl"]] = baseline.slope + baseline.se.slope * random.baseline
        }
      }
    }
    
    return(initial.values)
  })
}

