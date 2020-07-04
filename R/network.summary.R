pick.summary.variables <- function(result, extra.pars = NULL, only.pars = NULL){
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
    pars <- c("d", "sd", "sigma", "b_bl", "beta", "B", "sdB")
  } else{
    pars <- only.pars
  }
  if(!is.null(extra.pars)){
    pars <- c(pars, extra.pars)
  }
  summary.samples <- lapply(samples, function(x){x[,varnames.split %in% pars, drop = F]})
  summary.samples <- coda::mcmc.list(summary.samples)
  summary.samples
}

#' Summarize result run by \code{\link{network.run}}
#'
#' This function uses summary function in coda package to summarize mcmc.list object. Monte carlo error (Time-series SE) is also obtained using the coda package and is printed in the summary as a default.
#'
#' @param object Result object created by \code{\link{network.run}} function
#' @param ... Additional arguments affecting the summary produced
#' @return Returns summary of the network model result
#' @examples
#' network <- with(statins, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial",
#'  Treat.order = c("Placebo", "Statin"), covariate = covariate, covariate.type = "discrete")
#' })
#' \donttest{
#' result <- network.run(network)
#' summary(result)
#' }
#' @export

summary.network.result <- function(object, ...){

  if(!inherits(object, "network.result")) {
    stop('This is not the output from network.run. Need to run network.run function first')
  }
  summary.samples <- pick.summary.variables(object, ...)

  rval <- list("summary.samples"= summary(summary.samples),
               "Treat.order" =  object$network$Treat.order,
               "deviance" = unlist(object$deviance[1:3]),
               "total_n" = sum(object$network$na))
  class(rval) <- 'summary.network.result'
  rval
}

#' Plot traceplot and posterior density of the result
#'
#' This function uses plotting function in coda package to plot mcmc.list object
#'
#' @param x Result object created by \code{\link{network.run}} function
#' @param ... Additional arguments affecting the plot produced
#' @return None
#' @examples
#' network <- with(statins, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial",
#'  Treat.order = c("Placebo", "Statin"), covariate = covariate, covariate.type = "discrete")
#' })
#' \donttest{
#' result <- network.run(network)
#' plot(result, only.pars = "sd")
#' }
#' @export

plot.network.result <- function(x, ...) {
  summary.samples <- pick.summary.variables(x, ...)
  plot(summary.samples)
}

#' Use coda package to plot Gelman-Rubin diagnostic plot
#'
#' This function plots Gelman-Rubin diagnostic using coda package.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param extra.pars Extra parameters that the user wants to plot other than the default parameters.
#' @param only.pars Parameters that user wants to display. This gets rids of other default parameters user doesn't want to show.
#' @return None
#' @examples
#' network <- with(statins, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial",
#'  Treat.order = c("Placebo", "Statin"), covariate = covariate, covariate.type = "discrete")
#' })
#' \donttest{
#' result <- network.run(network)
#' network.gelman.plot(result)
#' }
#' @export

network.gelman.plot <- function(result, extra.pars = NULL, only.pars = NULL){
  summary.samples <- pick.summary.variables(result, extra.pars, only.pars)
  summary.samples <- mcmc.list(lapply(summary.samples, function(x) { x[,colSums(abs(x)) != 0] }))

  for(v in 1:nvar(summary.samples)){
    gelman.plot(summary.samples[,v,drop=FALSE])
  }
}

#' Use coda package to find Gelman-Rubin diagnostics
#'
#' This function uses coda package to find Gelman-Rubin diagnostics.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param extra.pars Extra parameters that the user wants to display other than the default parameters.
#' @param only.pars Parameters that user wants to display. This gets rids of other default parameters user doesn't want to show.
#' @return Returns gelman-rubin diagnostics
#' @examples
#' network <- with(statins, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial",
#'  Treat.order = c("Placebo", "Statin"), covariate = covariate, covariate.type = "discrete")
#' })
#' \donttest{
#' result <- network.run(network)
#' network.gelman.diag(result, extra.pars = "Eta")
#' }
#' @export

network.gelman.diag <- function(result, extra.pars = NULL, only.pars = NULL){
  summary.samples <- pick.summary.variables(result, extra.pars, only.pars)
  summary.samples <- mcmc.list(lapply(summary.samples, function(x) { x[,colSums(abs(x)) != 0] }))

  gelman.diag(summary.samples, multivariate = FALSE)$psrf
}

#' Generate autocorrelation diagnostics using coda package
#'
#' This function generates autocorrelation diagnostics using coda package. User can specify lags and parameters to display.
#' Note that to display extra parameters that are not saved, user needs to first specify parameters in \code{extra.pars.save} parameter in \code{\link{network.run}} function.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param lags A vector of lags at which to calculate the autocorrelation
#' @param extra.pars Extra parameters that the user wants to display other than the default parameters.
#' @param only.pars Parameters that user wants to display. This gets rids of other default parameters user doesn't want to show.
#' @return Returns autocorrelation diagnostics
#' @examples
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' \donttest{
#' result <- network.run(network)
#' network.autocorr.diag(result, only.pars = "d")
#' }
#' @export

network.autocorr.diag <- function(result, lags = c(0,1,5,10,50), extra.pars = NULL, only.pars = NULL){
  summary.samples <- pick.summary.variables(result, extra.pars, only.pars)
  summary.samples <- mcmc.list(lapply(summary.samples, function(x) { x[,colSums(abs(x)) != 0] }))

  autocorr.diag(summary.samples, lags = lags)
}

#' Generate autocorrelation plot using coda package
#'
#' This function plots autocorrelation using coda package.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param extra.pars Extra parameters that the user wants to plot other than the default parameters.
#' @param only.pars Parameters that user wants to display. This gets rids of other default parameters user doesn't want to show
#' @return None
#' @examples
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' \donttest{
#' result <- network.run(network)
#' network.autocorr.plot(result)
#' }
#' @export

network.autocorr.plot <- function(result, extra.pars = NULL, only.pars = NULL){
  summary.samples <- pick.summary.variables(result, extra.pars, only.pars)
  summary.samples <- mcmc.list(lapply(summary.samples, function(x) { x[,colSums(abs(x)) != 0] }))
  autocorr.plot(summary.samples)
}

#' Find relative effects for base treatment and comparison treatments
#'
#' This function calculates relative effects for base treatment and comparison treatments.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param base.treatment Base treatment user wants for the relative effects. Base treatment is initially set by \code{Treat.order} parameter in \code{\link{network.data}} (first one in the list). If set to null, default is to use base treatment.
#' @param comparison.treatments Treatments that user wants to compare against base treatment. If set to null, all the treatments besides base treatment is considered as comparison treatments.
#' @param base.category Base category user wants for the relative effects. Only used for multinomial data.
#' @param comparison.categories Category that user wants to compare against base.category. Only used for multinomial data.
#' @param covariate Covariate value at which to compute relative effects. Only used if covariate value is specified in the model.
#' @return
#' This returns a mcmc.list sample of relative effects for the base treatment specified. This allows user to obtain relative effects of different base.treatment after the sampling has been done.
#' For a simple summary, use \code{\link{relative.effects.table}}.
#' @examples
#' network <- with(parkinsons, {
#'  network.data(Outcomes, Study, Treat, SE = SE, response = "normal")
#' })
#' \donttest{
#' result <- network.run(network)
#' summary(relative.effects(result, base.treatment = "Placebo"))
#' }
#' @seealso \code{\link{relative.effects.table}}
#' @export

relative.effects <- function(result, base.treatment = NULL, comparison.treatments = NULL, base.category = NULL, comparison.categories = NULL, covariate = NULL){

  network <- result$network

  if(!is.null(covariate)){
    stopifnot(length(covariate) == dim(network$covariate)[2])
  }

  Treat.order <- network$Treat.order
  if(!is.null(base.treatment)){
    stopifnot(base.treatment %in% Treat.order)
  } else{
    base.treatment <- Treat.order[1]
  }
  if(!is.null(comparison.treatments)){
    stopifnot(comparison.treatments %in% Treat.order)
    stopifnot(!comparison.treatments %in% base.treatment)
  } else{
    comparison.treatments <- Treat.order[-which(Treat.order == base.treatment)]
  }
  if(!is.null(covariate)){
    summary.samples <- pick.summary.variables(result, only.pars = c("d", "beta"))
  } else{
    summary.samples <- pick.summary.variables(result, only.pars = c("d"))
  }
  vars <- dimnames(summary.samples[[1]])[[2]]

  if(network$response != "multinomial"){
    effects <- matrix(0, nrow = network$ntreat, ncol = length(comparison.treatments))
    effects[which(Treat.order == base.treatment),] = -1

    col_name = NULL
    for(i in 1:ncol(effects)){
      effects[which(comparison.treatments[i] == Treat.order),i] = 1
      col_name <- c(col_name, paste0("d_treatment", base.treatment, comparison.treatments[i]))
    }

    if(!is.null(covariate)){
      cov_matrix <-  covariate_centerered  <- NULL
      for(i in 1:length(covariate)){
        cov <- effects
        covariate_centered <- covariate[i] - network[[paste0("mx",i)]]
        cov <- cov * covariate_centered
        cov_matrix <- rbind(cov_matrix, cov)
      }
      effects <- rbind(cov_matrix, effects)
    }
    colnames(effects) <- col_name
    rownames(effects) <- vars

    samples <- as.mcmc.list(lapply(summary.samples, function(chain){
      samples <- chain %*% effects
      colnames(samples) <- colnames(effects)
      mcmc(samples, start = start(chain), end = end(chain), thin = thin(chain))
    }))
  } else{
    vars_d <- vars[grep("d\\[", vars)]
    categories_row <- as.numeric(substr(vars_d, nchar(vars_d[1])-1, nchar(vars_d[1])-1))
    categories_row <- categories_row+1
    ncat <- network$ncat

    if(!is.null(base.category)){
      stopifnot(base.category %in% 1:ncat)
    } else{
      base.category <- 1
    }
    if(!is.null(comparison.categories)){
      stopifnot(comparison.categories %in% 1:ncat)
      stopifnot(!comparison.categories %in% base.category)
    } else{
      comparison.categories <- (1:ncat)[-base.category]
    }

    effects <- matrix(0, nrow = network$ntreat*(network$ncat-1), length(vars), ncol = length(comparison.treatments) * length(comparison.categories))
    categories_column <- rep(comparison.categories, each = length(comparison.treatments))

    effects[which(rep(Treat.order, ncat-1) == base.treatment),] <- -1
    col_name <- NULL
    for(i in 1:ncol(effects)){
      effects[which(rep(Treat.order, ncat-1) == rep(comparison.treatments, length(comparison.categories))[i]),i] <- 1
      col_name <- c(col_name, paste0("d_treatment", base.treatment, rep(comparison.treatments, length(comparison.categories))[i]))
    }
    colnames(effects) <- col_name

    for(i in 1:ncol(effects)){
      effects[which(categories_row == base.category),i] <- -effects[which(categories_row == base.category),i]
      effects[which(categories_row != base.category & categories_row != rep(comparison.categories, each = length(comparison.treatments))[i]),i] <- 0
      colnames(effects)[i] <- paste0(colnames(effects)[i], "_category", base.category, rep(comparison.categories, each = length(comparison.treatments))[i])
    }

    if(!is.null(covariate)){
      cov_matrix <-  covariate_centerered  <- NULL
      for(i in 1:length(covariate)){
        cov <- effects
        covariate_centered <- covariate[i] - network[[paste0("mx",i)]]
        cov <- cov * covariate_centered
        cov_matrix <- rbind(cov_matrix, cov)
      }
      effects <- rbind(cov_matrix, effects)
    }
    rownames(effects) <- vars

    samples <- as.mcmc.list(lapply(summary.samples, function(chain){
      samples <- chain %*% effects
      colnames(samples) <- colnames(effects)
      mcmc(samples, start = start(chain), end = end(chain), thin = thin(chain))
    }))
  }
  samples
}

#' Make a summary table for relative effects
#'
#' This function creates a summary table of relative effects. Relative effects are in units of log odds ratio for binomial and multinomial data and real number scale for normal data.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param summary_stat Specifies what type of statistics user wants. Options are: "mean", "ci", "quantile", "sd", "p-value".
#' "ci" gives 95% confidence interval (0.025, 0.5, 0.975) and "quantile" gives specific quantile specified in probs parameter.
#' "p-value" is the probability relative effect (in binomial, log odds ratio) is less than 0.
#' @param probs Used only for the quantile summary. Specifies which quantile user wants the summary of (should be one numeric value between 0 to 1)
#' @param base.category Specifies for which base category user wants for the summary. Used only for multinoimal.
#' @return Returns relative effects table
#' @examples
#' #cardiovascular
#' network <- with(cardiovascular,{
#'  network.data(Outcomes, Study, Treat, N, response = "multinomial")
#' })
#' \donttest{
#' result <- network.run(network)
#' exp(relative.effects.table(result)) #look at odds ratio instead of log odds ratio
#' }
#' @seealso \code{\link{relative.effects}}
#' @export

relative.effects.table <- function(result, summary_stat = "mean", probs = NULL, base.category = NULL){

  stopifnot(summary_stat %in% c("mean", "quantile", "sd", "p-value", "ci"))

  if(!is.null(probs)){
    if(length(probs) != 1){
      stop("length of probs should be 1")
    }
  }

  Treat.order <- result$network$Treat.order

  ts <- 1:length(Treat.order)
  comps <- combn(ts, 2)

  if(result$network$response != "multinomial"){
    tbl <- matrix(NA, nrow = length(ts), ncol = length(ts), dimnames = list(Treat.order, Treat.order))

    for (i in 1:ncol(comps)) {
      comp <- comps[, i]
      samples <- as.matrix(relative.effects(result, base.treatment = Treat.order[comp[1]], comparison.treatments = Treat.order[comp[2]]))

      if(summary_stat == "mean"){
        tbl[comp[1], comp[2]] <- mean(samples)
        tbl[comp[2], comp[1]] <- -tbl[comp[1], comp[2]]
      } else if(summary_stat == "ci"){
        q <- round(quantile(samples, probs = c(0.025, 0.5, 0.975)), 6)
        tbl[comp[1], comp[2]] <- paste0("[", q[1], ",", q[2], ",", q[3], "]")
        tbl[comp[2], comp[1]] <- paste0("[", -q[3], ",", -q[2], ",", -q[1], "]")
      } else if(summary_stat == "quantile"){
        tbl[comp[1], comp[2]] <- round(quantile(samples, probs = probs), 6)
        tbl[comp[2], comp[1]] <- -tbl[comp[1], comp[2]]
      } else if(summary_stat == "sd"){
        tbl[comp[1], comp[2]] <- tbl[comp[2], comp[1]] <- sd(samples)
      } else if(summary_stat == "p-value"){
        tbl[comp[1], comp[2]] <- sum(samples < 0)/ dim(samples)[1]
        tbl[comp[2], comp[1]] <- 1 - tbl[comp[1], comp[2]]
      }
    }
  } else if(result$network$response == "multinomial"){
    ncat <- result$network$ncat
    tbl <- array(NA, dim = c(length(ts), length(ts), ncat -1), dimnames = list(Treat.order, Treat.order, NULL))

    for (i in 1:ncol(comps)) {
      comp <- comps[, i]
      samples <- as.matrix(relative.effects(result, base.treatment = Treat.order[comp[1]], comparison.treatments = Treat.order[comp[2]], base.category = base.category))

      if(summary_stat == "mean"){
        tbl[comp[1], comp[2],] <- apply(samples, 2, mean)
        tbl[comp[2], comp[1],]  <- -tbl[comp[1], comp[2],]
      } else if(summary_stat == "ci"){
        q <- round(apply(samples, 2, quantile, probs = c(0.025, 0.5, 0.975)), 6)
        q1 <- apply(q, 2, function(x){ paste0("[", x[1], ",", x[2], ",", x[3], "]")})
        q2 <- apply(q, 2, function(x){ paste0("[", -x[3], ",", -x[2], ",", -x[1], "]")})
        tbl[comp[1], comp[2],] <- q1
        tbl[comp[2], comp[1],] <- q2
      } else if(summary_stat == "quantile"){
        tbl[comp[1], comp[2],] <- apply(samples, 2, quantile, probs = probs)
        tbl[comp[2], comp[1],] <- -tbl[comp[1], comp[2],]
      } else if(summary_stat == "sd"){
        tbl[comp[1], comp[2],] <- tbl[comp[2], comp[1],] <- apply(samples, 2, sd)
      } else if(summary_stat == "p-value"){
        tbl[comp[1], comp[2],] <- apply(samples, 2, function(x){ sum(x <0) / length(x)})
        tbl[comp[2], comp[1],] <- 1 - tbl[comp[1], comp[2],]
      }
    }
  }
  tbl
}

#' Create a treatment rank table
#'
#' This function makes a table of ranking for each treament. Each number in the cell represents a probability certain treatment was in such rank.
#' This table is also stored as an output from \code{\link{network.run}}.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @return Returns a table of ranking
#' @examples
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' \donttest{
#' result <- network.run(network)
#' rank.tx(result)
#' }
#' @seealso \code{\link{network.rank.tx.plot}}
#' @export

rank.tx <- function(result){
  samples <- result[["samples"]]
  varnames <- dimnames(samples[[1]])[[2]]
  varnames.split <- sapply(strsplit(varnames, "\\["), '[[', 1)
  varnames.split <- gsub("[[:digit:]]","",varnames.split)

  rank.samples <- lapply(samples, function(x){x[,varnames.split %in% "prob"]})

  Treat.order <- result$network$Treat.order
  response <- result$network$response

  if(response != "multinomial"){
    prob.matrix <- matrix(NA, nrow = length(Treat.order), ncol = length(Treat.order), dimnames = list(paste0("rank ", 1:length(Treat.order)), paste0("treatment ", Treat.order)))
    for(i in 1:nrow(prob.matrix)){
      for(j in 1:ncol(prob.matrix)){
        prob.matrix[i,j] <- mean(unlist(lapply(rank.samples, function(x){ x[,paste0("prob[", i, ",", j, "]")]})))
      }
    }
  } else if(response == "multinomial"){
    ncat <- result$network$ncat
    prob.matrix <- array(NA, dim = c(length(Treat.order), length(Treat.order), ncat-1), dimnames = list(paste0("rank ", 1:length(Treat.order)), paste0("treatment ", Treat.order),  paste0("Category ", 1:(ncat-1))))
    for(i in 1:nrow(prob.matrix)){
      for(j in 1:ncol(prob.matrix)){
        for(k in 1:(ncat-1)){
          prob.matrix[i,j,k] <- mean(unlist(lapply(rank.samples, function(x){ x[,paste0("prob[", i, ",", j, ",", k, "]")]})))
        }
      }
    }
  }
  return(prob.matrix)
}



#' Create a treatment rank plot
#'
#' This plot displays how each treatment is ranked. For each rank, we show how likely each treatment will be at that rank.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param txnames Treatment names used in creating legend
#' @param catnames Category names. Only used in multinomial.
#' @param legend.position x,y position of the legend
#' @return None
#' @examples
#' network <-with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' \donttest{
#' result <- network.run(network)
#' network.rank.tx.plot(result, txnames = c("a", "b"))
#' }
#' @seealso \code{\link{rank.tx}}
#' @export

network.rank.tx.plot <- function(result, txnames = NULL, catnames = NULL, legend.position = c(1,1)){

  rank.table <- rank.tx(result)
  ntreat = dim(rank.table)[1]
  if (is.null(txnames)) txnames <- paste("Treatment", result$network$Treat.order)

  if(result$network$response != "multinomial"){
    plot(seq(ntreat),seq(ntreat),type="n",xaxt="n",ylim=c(0,1),pty="s",yaxt="n",ylab="Probability",xlab="Rank")
    axis(side=1,at=seq(ntreat))
    axis(side=2,at=seq(0,1,by=0.2))
    for (i in seq(ntreat)) {
      points(seq(ntreat), rank.table[,i],type="b",lty=i,col=i,pch=i)
    }
    legend(legend.position[1], legend.position[2],txnames,lty=1:ntreat,bty="n",cex=.75,col=1:ntreat)
  } else if(result$network$response == "multinomial"){
    ncat <- dim(rank.table)[3]
    if (is.null(catnames)) catnames <- paste("Outcome Category with base 1 and comparison", 1+seq(ncat))
    for (j in seq(ncat)) {
      plot(seq(ntreat),seq(ntreat),type="n",xaxt="n",ylim=c(0,1),pty="s",yaxt="n",ylab="Probability",xlab="Rank")
      axis(side=1,at=seq(ntreat))
      axis(side=2,at=seq(0,1,by=0.2))
      title(catnames[j])
      for (i in seq(ntreat)) {
        points(seq(ntreat), rank.table[,i,j],type="b",lty=i,col=i,pch=i)
      }
      legend(legend.position[1], legend.position[2],txnames,lty=1:ntreat,bty="n",cex=.75,col=1:ntreat)
    }
  }
}



#' Create a treatment cumulative rank plot
#'
#' This function creates a treatment cumulative rank plot. Rank preference can be specified by the \code{rank.preference} parameter in \code{\link{network.data}}
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param txnames Treatment names used in creating legend
#' @param catnames Category names. Only used in multinomial.
#' @param legend.position x, y position of the legend
#' @return None
#' @examples
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' \donttest{
#' result <- network.run(network)
#' network.cumrank.tx.plot(result, txnames = c("control", "beta blocker"))
#' }
#' @seealso \code{\link{rank.tx}}
#' @export

network.cumrank.tx.plot <- function(result, txnames = NULL, catnames = NULL, legend.position = c(1,1)){
  rank.table <- rank.tx(result)
  ntreat = dim(rank.table)[1]
  if (is.null(txnames)) txnames <- paste("Treatment", result$network$Treat.order)

  if(result$network$response != "multinomial"){
    x <- apply(rank.table,2,cumsum)
    plot(seq(ntreat),seq(ntreat),type="n",xaxt="n",ylim=c(0,1),yaxt="n",ylab="Cumulative Probability",xlab="Rank")
    axis(side=1,at=seq(ntreat))
    axis(side=2,at=seq(0,1,by=0.2))
    for (j in seq(ntreat))
      points(seq(ntreat), x[,j],type="l",lty=j,col=j)
    legend(legend.position[1], legend.position[2], txnames,lty=1:(ntreat),bty="n",cex=.75,col=1:(ntreat))
  } else if(result$network$response == "multinomial"){
    ncat <- dim(rank.table)[3]
    if (is.null(catnames)) catnames <- paste("Outcome Category with base 1 and comparison", 1+seq(ncat))

    for (i in seq(ncat))  {
      x = apply(rank.table[,,i],2,cumsum)
      plot(seq(ntreat),seq(ntreat),type="n",xaxt="n",ylim=c(0,1),yaxt="n",ylab="Cumulative Probability",xlab="Rank")
      axis(side=1,at=seq(ntreat))
      axis(side=2,at=seq(0,1,by=0.2))
      title(catnames[i])
      for (j in seq(ntreat))
        points(seq(ntreat), x[,j],type="l",lty=j,col=j)
      legend(legend.position[1], legend.position[2],txnames,lty=1:ntreat,bty="n",cex=.75,col=1:ntreat)
    }
  }
}

#' Calculate SUCRA
#'
#' SUCRA is the surface under the cumulative ranking distribution defined in Salanti et al. (2011)
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param txnames Treatment names used in creating legend
#' @param catnames Category names. Only used in multinomial.
#' @return Returns SUCRA for each treatment
#' @examples
#' ########### certolizumab (with baseline risk)
#' network <- with(certolizumab, {
#'  network.data(Outcomes, Study, Treat, N=N, response = "binomial", Treat.order,
#'  baseline = "common", hy.prior = list("dhnorm", 0, 9.77))
#' })
#' \donttest{
#' result <- network.run(network)
#' sucra(result)
#' }
#' @seealso \code{\link{rank.tx}}
#' @references G. Salanti, A.E. Ades, J.P.A. Ioannidisa (2011), \emph{Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial}, Journal of Clinical Epidemiology 64(2):163-71. [\url{https://doi.org/10.1016/j.jclinepi.2010.03.016}]
#' @export

sucra = function(result, txnames = NULL, catnames = NULL)
{
  rank.table <- rank.tx(result)
  ntreat = dim(rank.table)[1]
  if (is.null(txnames)) txnames <- paste("Treatment", result$network$Treat.order)

  if(result$network$response != "multinomial"){
    if(ntreat ==2){
      x <- rank.table[-ntreat,]
    } else{
      x <- apply(apply(rank.table[-ntreat,],2,cumsum),2,sum)/(ntreat-1)
    }
    names(x) <- txnames
  } else if(result$network$response == "multinomial"){
    ncat <- dim(rank.table)[3]
    if (is.null(catnames)) catnames <- paste("Outcome Category with base 1 and comparison", 1+seq(ncat))
    x <- array(NA,dim(rank.table)[2:3])
    for (i in seq(ncat)){
      if(ntreat ==2){
        x[,i] <-   rank.table[-ntreat,,i]
      } else{
        x[,i] <- apply(apply(rank.table[-ntreat,,i],2,cumsum),2,sum)/(ntreat-1)
      }
      dimnames(x) <- list(txnames,catnames)
    }
  }
  return(x)
}

#################### Deviance calculation and plots

#' Find deviance statistics such as DIC and pD.
#'
#' Calculates deviance statistics. This function automatically called in \code{\link{network.run}} and the deviance statistics are stored after sampling is finished.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @return
#' \item{Dbar}{Overall residual deviance}
#' \item{pD}{Sum of leverage_arm (i.e. total leverage)}
#' \item{DIC}{Deviance information criteria (sum of Dbar and pD)}
#' \item{data.points}{Total number of arms in the meta analysis}
#' \item{dev_arm}{Posterior mean of the residual deviance in each trial arm}
#' \item{devtilda_arm}{Deviance at the posterior mean of the fitted values}
#' \item{leverage_arm}{Difference between dev_arm and devtilda_arm for each trial}
#' \item{rtilda_arm}{Posterior mean of the fitted value for binomial and multinomial}
#' \item{ybar_arm}{Posterior mean of the fitted value for normal}
#' @examples
#' #parkinsons
#' network <- with(parkinsons, {
#'  network.data(Outcomes, Study, Treat, SE = SE, response = "normal")
#' })
#' \donttest{
#' result <- network.run(network)
#' calculate.deviance(result)
#' }
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. [\url{https://doi.org/10.1177/0272989X12458724}]
#' @export

calculate.deviance <- function(result){

  network <- result$network
  samples <- result$samples

  totresdev <- lapply(samples, function(x){ x[,"totresdev"]})
  Dbar <- mean(unlist(totresdev))

  ###### find residual deviance by arm
  if(network$response == "multinomial" & !is.null(network$miss.matrix)){
    dev <- list()
    for(ii in seq(network$npattern)){
      dev_each <- lapply(samples, function(x) { x[,grep(paste0("dev", ii, "\\["), dimnames(samples[[1]])[[2]])]})
      dev_each <- do.call(rbind, dev_each)
      dev_each <- apply(dev_each, 2, mean)

      n_value <- network[[paste0("n", ii)]]
      dev_matrix <- matrix(NA, nrow = dim(n_value)[1], ncol = dim(n_value)[2])

      for(i in 1:dim(dev_matrix)[1]){
        for(j in 1:dim(dev_matrix)[2]){
          ind <- which(paste0("dev", ii, "[", i, ",", j, "]") == names(dev_each))
          if(length(ind) != 0){
            dev_matrix[i,j] <- dev_each[ind]
          }
        }
      }
      dev[[paste0("dev", ii)]] <- dev_matrix
    }
    dev_arm <- do.call(rbind, dev)
  } else{
    dev <- lapply(samples, function(x) { x[,grep("dev\\[", dimnames(samples[[1]])[[2]])]})
    dev <- do.call(rbind, dev)
    dev <- apply(dev, 2, mean)

    dev_matrix <- matrix(NA, nrow =  network$nstudy, ncol = max(network$na))
    for(i in 1:dim(dev_matrix)[1]){
      for(j in 1:dim(dev_matrix)[2]){
        ind <- which(paste("dev[", i, ",", j, "]", sep = "") == names(dev))
        if(length(ind) != 0){
          dev_matrix[i,j] <- dev[ind]
        }
      }
    }
    dev_arm <- dev_matrix
  }

  ############find leverage
  if(network$response == "binomial"){

    rtilda <- lapply(samples, function(x){ x[,grep("rhat\\[", dimnames(samples[[1]])[[2]])] })
    rtilda <- do.call(rbind, rtilda)
    rtilda <- apply(rtilda, 2, mean)

    rtilda_arm <- devtilda_arm <- matrix(NA, nrow = network$nstudy, ncol = max(network$na))
    for(i in 1:network$nstudy){
      for(j in 1:network$na[i]){
        r_value <- network$r[i,j]
        n_value <- network$n[i,j]
        rtilda_arm[i,j] <- rtilda[which(paste("rhat[", i, ",", j, "]", sep = "") == names(rtilda))]

        devtilda_arm[i,j] <- ifelse(r_value != 0, 2 * r_value * (log(r_value)-log(rtilda_arm[i,j])), 0)
        devtilda_arm[i,j] <- devtilda_arm[i,j] + ifelse((n_value - r_value) != 0, 2 * (n_value-r_value) *(log(n_value-r_value) - log(n_value- rtilda_arm[i,j])), 0)
      }
    }
  } else if(network$response == "normal"){
    ybar <- lapply(samples, function(x){ x[,grep("theta\\[", dimnames(samples[[1]])[[2]])] })
    ybar <- do.call(rbind, ybar)
    ybar <- apply(ybar, 2, mean)

    ybar_arm <- devtilda_arm <- matrix(NA, nrow = network$nstudy, ncol = max(network$na))

    for(i in 1:network$nstudy){
      for(j in 1:network$na[i]){
        r_value <- network$r[i,j]
        se_value <- network$se[i,j]

        if(class(network) == "nodesplit.network.data"){
          ybar_arm[i,j] <- ybar[which(paste("theta[", i, ",", network$t[i,j], "]", sep = "") == names(ybar))]
        } else{
          ybar_arm[i,j] <- ybar[which(paste("theta[", i, ",", j, "]", sep = "") == names(ybar))]
        }
        devtilda_arm[i,j] <- ifelse(se_value != 0, (r_value - ybar_arm[i,j])^2 / se_value^2, 0)
      }
    }
  } else if(network$response == "multinomial"){
    if(is.null(network$miss.matrix)){ #complete dataset
      rtilda <- lapply(samples, function(x){ x[,grep("rhat\\[", dimnames(samples[[1]])[[2]])]})
      rtilda <- do.call(rbind, rtilda)
      rtilda <- apply(rtilda, 2, mean)

      rtilda_arm <- devtilda_category <- array(NA, dim = c(network$nstudy, max(network$na), network$ncat))
      for(i in 1:network$nstudy){
        for(j in 1:network$na[i]){
          for(k in 1:network$ncat){
            r_value <- network$r[i,j,k]
            rtilda_arm[i,j,k] <- rtilda[which(paste("rhat[", i, ",", j, ",", k, "]", sep = "") == names(rtilda))]
            devtilda_category[i,j,k] <- ifelse(r_value != 0,  2 * r_value * log(r_value/rtilda_arm[i,j,k]), 0)
          }
        }
      }
      devtilda_arm <- apply(devtilda_category, 1:2, sum)
    } else{ #incomplete datacase
      devtilda_value <- rtilda_arm <- list()
      for(ii in seq(network$npattern)){
        r_values <- network[[paste0("r",ii)]]
        devtilda_category <- rtilda_matrix <- array(NA, dim = dim(r_values))

        rtilda <- lapply(samples, function(x){ x[,grep(paste0("rhat", ii, "\\["), dimnames(samples[[1]])[[2]])] })
        rtilda <- do.call(rbind, rtilda)
        rtilda <- apply(rtilda, 2, mean)

        for(i in 1:dim(r_values)[1]){
          for(j in 1:dim(r_values)[2]){
            for(k in 1:dim(r_values)[3]){
              found <- which(paste("rhat", ii, "[", i, ",", j, ",", k, "]", sep = "") == names(rtilda))
              r_value <- r_values[i,j,k]
              if(!is.na(r_value) & length(found) != 0){
                rtilda_matrix [i,j,k] <- rtilda[found]
                devtilda_category[i,j,k] <- ifelse(r_value != 0,  2 * r_value * log(r_value/rtilda_matrix[i,j,k]), 0)
              }
            }
          }
        }
        devtilda_matrix <- apply(devtilda_category, 1:2, sum)
        rtilda_arm[[ii]] <- rtilda_matrix
        devtilda_value[[ii]] <- devtilda_matrix
      }
      devtilda_arm <- do.call(rbind, devtilda_value)
    }
  }
  leverage_arm <- dev_arm - devtilda_arm
  pD <- sum(leverage_arm, na.rm = TRUE)
  DIC <- Dbar + pD

  out <- list(Dbar = Dbar, pD = pD, DIC = DIC, data.points = sum(network$na), dev_arm = dev_arm, devtilda_arm = devtilda_arm, leverage_arm = leverage_arm)
  if(network$response == "binomial" || network$response == "multinomial"){
    out$rtilda_arm = rtilda_arm
  } else if(network$response == "normal"){
    out$ybar_arm = ybar_arm
  }
  return(out)
}


#' Make a deviance plot
#'
#' This makes a deviance plot which plots residual deviance (dev_arm) vs. all the arms for each study.
#' @param result Object created by \code{\link{network.run}} function
#' @return None
#' @examples
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' \donttest{
#' result <- network.run(network)
#' network.deviance.plot(result)
#' }
#' @export

network.deviance.plot <- function(result){
  deviance <- result$deviance
  dev_vector <- c(t(deviance$dev_arm))
  dev_vector <- dev_vector[!is.na(dev_vector)]
  plot(seq(sum(result$network$na)), dev_vector, xlab = "Arm", ylab = "Residual Deviance", main = "Per-arm residual deviance")
}

#' Make a leverage plot
#'
#' This function makes a leverage vs. square root of residual deviance plot
#'
#' @param result Object created by \code{\link{network.run}} function
#' @return None
#' @examples
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' \donttest{
#' result <- network.run(network)
#' network.leverage.plot(result)
#' }
#' @export

network.leverage.plot <- function(result){
  deviance <- result$deviance
  dev <- sqrt(apply(deviance$dev_arm, 1, mean, na.rm = TRUE))
  leverage <- apply(deviance$leverage, 1, mean, na.rm = TRUE)
  plot(dev, leverage, xlim = c(0, max(c(dev, 2.5))), ylim = c(0, max(c(leverage,4))),
       xlab = "Square root of residual deviance", ylab = "Leverage", main = "Leverage versus residual deviance")
  mtext("Per-study mean per-datapoint contribution")
}

#' Make a covariate plot
#'
#' This function makes a covariate plot of how the relative effect changes as the covariate value changes.
#' User needs to specify one base treatment and one comparison treatment to make this plot (base category and comparison category is also needed for multinomial).
#' The function uses the \code{\link{relative.effects}} to calculate the correct relative effect. 2.5\%, median, and 97.5\% C.I. are drawn.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param base.treatment Base treatment for relative effect
#' @param comparison.treatment Treatment comparing against base treatment
#' @param base.category Base category for multinomial data. Note that category in multinomial denotes which column it is in the Outcomes matrix. Thus, this should be a numeric value.
#' @param comparison.category Comparison category for multinomial data
#' @param covariate.name A vector of covariate names of the covariate that goes into x-axis label
#' @return None
#' @examples
#' ########### certolizumab (with covariate)
#' network <- with(certolizumab, {
#'  network.data(Outcomes, Study, Treat, N=N, response="binomial", Treat.order,
#'  covariate = covariate, hy.prior = list("dhnorm", 0, 9.77))
#' })
#' \donttest{
#' result <- network.run(network)
#' network.covariate.plot(result, base.treatment = "Placebo", comparison.treatment = "CZP",
#' covariate.name = "Disease Duration")
#' }
#' @export

network.covariate.plot <- function(result, base.treatment = NULL, comparison.treatment= NULL, base.category = NULL, comparison.category = NULL, covariate.name = NULL){

  if(is.null(network$covariate)){
    stop("need to provide covariate information to make this plot")
  }
  if(result$network$response != "multinomial"){
    if(is.null(base.treatment) || is.null(comparison.treatment)){
      stop("need to specify both base.treatment and comparison.treatment")
    }
  } else{
    if(is.null(base.treatment) || is.null(comparison.treatment) || is.null(base.category) || is.null(comparison.category)){
      stop("need to specify all base.treatment, comparison.treatment, base.category, and comparison.category")
    }
  }

  network <- result$network
  observed <- network$covariate
  xvals <- matrix(NA, nrow = dim(network$covariate)[2], ncol = 7)
  xlim <- matrix(NA, nrow = dim(network$covariate)[2], ncol = 2)
  covariate_mx <- NULL
  for(i in 1:dim(network$covariate)[2]){
    xlim[i,] <- c(min(observed[,i], na.rm = TRUE), max(observed[,i], na.rm = TRUE))
    xvals[i,] <- seq(xlim[i,1], xlim[i,2], length.out = 7)
    covariate_mx <- c(covariate_mx, network[[paste0("mx",i)]])
  }

  for(i in 1:dim(network$covariate)[2]){
    res <- lapply(xvals[i,], function(xval) {
      covariate <- covariate_mx
      covariate[i] <- xval
      if(network$response != "multinomial"){
        samples <- relative.effects(result, base.treatment, comparison.treatment, covariate = covariate)
      } else{
        samples <- relative.effects(result, base.treatment, comparison.treatment, base.category, comparison.category, covariate = covariate)
      }

      samples <- as.matrix(samples)
      stats <- t(apply(samples, 2, quantile, probs = c(0.025, 0.5, 0.975)))
      data.frame(median = stats[,"50%"], lower = stats[,"2.5%"], upper = stats[,"97.5%"])
    })
    res <- do.call(rbind,res)

    dim_names <- if(network$response != "multinomial"){
      dimnames(as.matrix(relative.effects(result, base.treatment, comparison.treatment)))[[2]]
    } else{
      dimnames(as.matrix(relative.effects(result, base.treatment, comparison.treatment, base.category, comparison.category)))[[2]]
    }

    ylim <- c(min(res), max(res))
    xlab_name <- ifelse(is.null(covariate.name), paste0("covariate ", i), covariate.name[i])

    plot(xvals[i,], res$median, type = "l", xlim = xlim[i,], ylim = ylim, main = "Treatment effect vs. covariate", xlab = xlab_name, ylab = dim_names)
    lines(xvals[i,], res$lower, lty = 2)
    lines(xvals[i,], res$upper, lty = 2)
  }
}

#' Calculate correlation matrix for multinomial heterogeneity parameter.
#'
#' This function calculates correlation matrix from the variance matrix for heterogeneity parameter. Only used for multinomial.
#' @param result Object created by \code{\link{network.run}} function
#' @return Returns correlation matrix
#' @examples
#' #cardiovascular
#' network <- with(cardiovascular, {
#'  network.data(Outcomes, Study, Treat, N, response = "multinomial")
#' })
#' \donttest{
#' result <- network.run(network)
#' variance.tx.effects(result)
#' }
#' @export

variance.tx.effects = function(result)
{
  if(result$network$response != "multinomial"){
    stop("this function is used only for multinomial response")
  }
  samples_sigma <- pick.summary.variables(result, only.pars = c("sigma"))
  samples_sigma <- do.call(rbind, samples_sigma)
  samples_sigma <- apply(samples_sigma, 2, mean)

  sigma_matrix <- matrix(samples_sigma, nrow = result$network$ncat-1)
  cor_matrix <- sigma_matrix/outer(sqrt(diag(sigma_matrix)),sqrt(diag(sigma_matrix)))

  return(list(sigma_matrix = sigma_matrix, cor_matrix = cor_matrix))
}

#' Draws forest plot
#'
#' Draws forest plot of pooled treatment effect. Reports odds ratio for binomial and multinomial outcomes and continuous scale for normal outcomes.
#'
#' @param result Object created by \code{\link{network.run}} function
#' @param level Confidence level. Default is 0.95 denoting 95 percent C.I.
#' @param ticks.position Position of the x-axis tick marks. If left unspecified, the function tries to set it at sensible values
#' @param label.multiplier This is a multiplying factor to move the position of the text associated with median[lower, upper] values. This number is multiplied by the range of x-axis and added to the x-axis limit. Default multiplier is set to 0.2.
#' @param label.margin This is how much margin space you specify to assign space for the median[lower, upper] values. Default margin is set to 10.
#' @param title Header name which you can modify
#' @param only.reference.treatment Indicator for plotting only the comparison to the reference treatment
#' @return None
#' @examples
#' network <- with(certolizumab, {
#'  network.data(Outcomes, Study, Treat, N=N, response="binomial", Treat.order,
#'  covariate = covariate, hy.prior = list("dhnorm", 0, 9.77))
#' })
#' \donttest{
#' result <- network.run(network)
#' network.forest.plot(result)
#' }
#' @references W. Viechtbauer (2010), \emph{Conducting meta-analyses in R with the metafor package}, Journal of Statistical Software, 36(3):1-48. [\url{https://doi.org/10.18637/jss.v036.i03}]
#' @export

network.forest.plot <- function(result, level = 0.95, ticks.position = NULL, label.multiplier = 0.2, label.margin = 10, title = "Network Meta-analysis Forest plot", only.reference.treatment = FALSE){

  ncat <- ifelse(result$network$response == "multinomial", result$network$ncat, 2)

  for(i in 1:(ncat-1)){

    if(i != 1) grid::grid.newpage()

    if(result$network$response == "multinomial"){
      lower <- relative.effects.table(result, summary_stat = "quantile", probs = (1- level)/2)[,,i]
      OR <- relative.effects.table(result, summary_stat = "quantile", probs = 0.5)[,,i]
      upper <- relative.effects.table(result, summary_stat = "quantile", probs = level + (1- level)/2)[,,i]
    } else{
      lower <- relative.effects.table(result, summary_stat = "quantile", probs = (1- level)/2)
      OR <- relative.effects.table(result, summary_stat = "quantile", probs = 0.5)
      upper <- relative.effects.table(result, summary_stat = "quantile", probs = level + (1- level)/2)
    }

    if(only.reference.treatment == TRUE){
      lower <- lower[1,-1]
      OR <- OR[1,-1]
      upper <- upper[1,-1]
    } else{
      lower <- -lower[lower.tri(lower)]
      OR <- -OR[lower.tri(OR)]
      upper <- -upper[lower.tri(upper)]
    }

    odds <- data.frame(lower = lower, OR = OR, upper = upper)

    if(result$network$response %in% c("binomial", "multinomial")){
      odds <- exp(odds) #report odds ratio instead of log odds ratio
    }

    Treat.order <- result$network$Treat.order
    ts <- 1:length(Treat.order)
    comps <- combn(ts, 2)

    name <- rep(NA, ncol(comps))
    for(j in 1:ncol(comps)){
      name[j] <- paste0(Treat.order[comps[2,j]]," vs ", Treat.order[comps[1,j]])
    }

    if(only.reference.treatment == TRUE){
      name <- name[1:(length(Treat.order)-1)]
    }
    odds$name <- name

    if(is.null(ticks.position)){
      if(result$network$response %in% c("binomial", "multinomial")){
        ticks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
      } else if(result$network$response == "normal"){
        ticks <- pretty(c(min(odds$lower, na.rm =TRUE), max(odds$upper, na.rm = TRUE)))
      }
    } else{
      ticks <- ticks.position
    }

    p <- ggplot(odds, aes(y = OR, x = name)) +
      geom_point() +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
      scale_x_discrete(limits = name) +
      geom_hline(yintercept = 1, linetype = 2) +
      coord_flip() +
      theme_bw() +
      theme(plot.margin = unit(c(1,label.margin,1,1), "lines"))

    if(result$network$response %in% c("binomial")){
      p <- p + labs(x = "Treatment comparison", y = "Odds Ratio", title = title) +
        scale_y_log10(breaks = ticks, labels = ticks)
    } else if(result$network$response %in% c("multinomial")){
      p <- p + labs(x = "Treatment comparison", y = "Odds Ratio", title = paste0(title, ": Multinomial Category ", (i+1), " vs 1")) +
        scale_y_log10(breaks = ticks, labels = ticks)
    } else if(result$network$response %in% c("normal")){
      p <- p + labs(x = "Treatment comparison", y = "Continuous Scale", title = title) +
        scale_y_continuous(breaks = ticks, labels = ticks)
    }

    #find actual xlim range; this part of code keeps changing with ggplot update..
    xlim.range <- ggplot_build(p)$layout$panel_params[[1]]$x.range

    p <- p + geom_text(aes(label = paste0(sprintf("%0.2f", round(OR, digits = 2)), " [", sprintf("%0.2f", round(lower, digits = 2)) , ", ", sprintf("%0.2f", round(upper, digits = 2)), "]")), y = xlim.range[2] + diff(xlim.range)*label.multiplier, x = 1:length(name))   # hjust = -1, vjust = 2)

    median_name_location <- ifelse(length(odds[,1]) <= 3, length(name) + 0.5, length(name) + 1)
    p <- p + geom_text(aes(label = "Median [95% Crl]"), y = xlim.range[2] + diff(xlim.range)*label.multiplier, x = median_name_location)

    gt <- ggplot_gtable(ggplot_build(p))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    grid::grid.draw(gt)
  }
}


#' Draws network graph using igraph package
#'
#' This function draws network graph using igraph package
#' @param network Object created by \code{\link{network.data}} function
#' @param label.dist distance of the label from the node. Default is 2.
#' @return None
#' @examples
#' #cardiovascular
#' network <- with(thrombolytic, {
#'  network.data(Outcomes, Study, Treat, N=N, response = "binomial")
#' })
#' draw.network.graph(network)
#' @export

draw.network.graph = function(network, label.dist = 2){

  if(class(network) == "contrast.network.data"){
    Treat <- c(t(network$Treat))[!is.na(c(t(network$Treat)))]
    Study <- rep(1:length(network$na), times = network$na)
  } else{
    Treat <- network$Treat.order[network$Treat]
    Study <- network$Study
  }
  pairs <- do.call(rbind, lapply(split(Treat, Study),
                                 function(x) t(combn(x,2))))
  pairs <- aggregate(rep(1, length(X1)) ~ X1 + X2, data = data.frame(pairs), sum)
  colnames(pairs)[3] <- "freq"
  g <- igraph::graph.edgelist(as.matrix(pairs[,1:2]), directed=FALSE)
  plot(g, edge.curved=FALSE, edge.width=pairs$freq, vertex.label.dist= label.dist)
}
