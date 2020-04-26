#' Make a network object containing data, priors, and a JAGS model file
#'
#' This function makes a network object that can be used to run network meta-analysis using \code{\link{network.run}}.
#' User needs to specify Outcomes, Study, Treat, N or SE, and response.
#' Prior parameters are filled in automatically based on the data type if not specified.
#' The input data should be arm-level so that we have observations for each treatment in each study.
#' The input data is preprocessed to fit the format necessary to run model in JAGS.
#'
#' @param Outcomes Arm-level outcomes. If it is a multinomial response, the matrix would have dimensions treatment arms (row) by multinomial categories (column). If it is binomial or normal, it would be a vector.
#' @param Study A vector of study indicator for each arm
#' @param Treat A vector of treatment indicator for each arm
#' @param N A vector of total number of observations in each arm. Used for binomial and multinomial responses.
#' @param SE A vector of standard error for each arm. Used only for normal response.
#' @param response Specification of the outcomes type. Must specify one of the following: "normal", "binomial", or "multinomial".
#' @param Treat.order Treatment order which determines how treatments are compared. The first treatment that is specified is considered to be the baseline treatment. Default order is alphabetical. If the treatments are coded 1, 2, etc, then the treatment with a value of 1 would be assigned as a baseline treatment.
#' @param type Type of model fitted: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @param rank.preference Set it equal to "higher" if higher values are preferred (i.e. assumes events are good). Set it equal to "lower" if lower values are preferred (i.e. assumes events are bad). Default is "higher".
#' @param baseline Three different assumptions for treatment x baseline risk interactions (slopes): "independent", "common", or "exchangeable". Default is "none" which doesn't incorporate baseline risk.
#' @param baseline.risk Two different assumptions for baseline risk: "independent" or "exchangeable". See Achana et al. (2012) for more information about baseline risk.
#' @param covariate A covariate matrix with each row representing each trial and column representing each covariate. This is a study-level data, meaning that the user doesn't need to repeatedly specify covariates for each arm.
#' @param covariate.type Should be a vector indicating the type of the covariate. Covariate can be either "continuous" or "discrete". If it continuous, covariates are centered. If the covariate is discrete it is not centered and it has to be in a dummy integer format (i.e. 0,1,2,...). The code doesn't factor the covariates for the user, so user needs to specify dummy variables if factor is needed.
#' @param covariate.model "independent" allows covariate effects for each treatment. "common" restricts same covariate effect for all treatment. Lastly, "exchangeable"  assumes that the covariate effects are different but related and strength is borrowed across them. We set "common" to be default. See Cooper et al. (2009) for more details on covariates.
#' @param mean.d Prior mean for the relative effect
#' @param prec.d Prior precision for the relative effect
#' @param mean.Eta Prior mean for the study effect (baseline risk)
#' @param prec.Eta Prior precision for the study effect (baseline risk)
#' @param hy.prior.Eta Between treatment heterogeneity in baseline risk (for exchangeable assumption only). Format of the parameter is same as hy.prior.
#' @param mean.bl Prior mean for the baseline slope
#' @param prec.bl Prior precision for the baseline slope
#' @param hy.prior.bl Between treatment heterogeneity in baseline slope (for exchangeable regression coefficient only). Format of the parameter is same as hy.prior.
#' @param mean.cov Prior mean for the covariate effect
#' @param prec.cov Prior precision for the covariate effect
#' @param hy.prior.cov Between treatment heterogeneity in covariate effect (for exchangeable regression coefficient only). Format of the parameter is same as hy.prior. Default is set to be dunif(0, 5) for binary, dunif(0, 100) for normal, and wishart with identity scale matrix and (# of categories - 1) degrees of freedom for multinomial.
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal and binomial response and wishart for multinomial response. It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm, dwish) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter. For wishart distribution, the last two parameter would be the scale matrix and the degrees of freedom.
#' @param mean.A Mean effect of 'standard' treatment (i.e. placebo) in a logit scale; this is used for binomial outcome when the risk difference, relative risk, or number needed to treat is needed; this should be informed from external evidence or can be found by meta-analyzing single proportions; For number needed to treat, we assume that events are "good". Reversal of sign is needed if the events are "bad".
#' @param prec.A Precision of 'standard' treatment in a logit scale
#' @return Creates list of variables that are used to run the model using \code{\link{network.run}}
#' \item{data}{Data combining all the input data. User can check this to insure the data is correctly specified. For modelling purposes, character valued studies or treatment variables are changed to numeric values based on alphabetical order.}
#' \item{nrow}{Total number of arms in the meta-analysis}
#' \item{ncat}{Number of columns in the Outcomes. Will equal 1 for binary and normal and number of categories for multinomial}
#' \item{nstudy}{Number of study}
#' \item{na}{Number of arms for each study}
#' \item{ntreat}{Number of treatment}
#' \item{b.id}{Indicator in sequence of all treatments for which treatment is base treatment in Study}
#' \item{t}{\code{Treat} transformed into a matrix which has dimensions number of study by max number of arms in studies}
#' \item{r}{\code{Outcomes} made into an array that is suitable for use in rjags code. For multinomial, it has 3 dimensions: number of study by max number of arms in studies by number of categories.}
#' \item{mx}{If the continuous covariate is included, it calculates the mean of the covariates which is used to center the covariates. The numeric indicator after mx refers to column number of the covariates if there are more than one covariates included. Discrete covariates are not centered.}
#' \item{mx_bl}{If the baseline effect is specified, it also calculates the mean baseline risk.}
#' \item{prior.data}{Prior data created using the user inputs or default values. If no user input is specifies for the prior, it uses default values.}
#' \item{code}{Rjags model file code that is generated using information provided by the user. To view model file inside R in a nice format, use \code{cat(network$code).}}
#' @examples
#' ###Blocker data example
#' blocker
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' network
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. [\url{https://doi.org/10.1177/0272989X12458724}]
#' @references F.A. Achana, N.J. Cooper, S. Dias, G. Lu, S.J.C. Rice, D. Kendrick, A.J. Sutton (2012), \emph{Extending methods for investigating the relationship between treatment effect and baseline risk from pairwise meta-analysis to network meta-analysis}, Statistics in Medicine 32(5):752-771. [\url{https://doi.org/10.1002/sim.5539}]
#' @references N.J. Cooper, A.J. Sutton, D. Morris, A.E. Ades, N.J. Welton (2009), \emph{Addressing between-study heterogeneity and inconsistency in mixed treatment comparisons: Application to stroke prevention treatments in individuals with non-rheumatic atrial fibrillation}, Statistics in Medicine 28:1861-1881. [\url{https://doi.org/10.1002/sim.3594}]
#' @export

network.data <- function(Outcomes = NULL, Study = NULL, Treat = NULL, N = NULL, SE = NULL, response = NULL, Treat.order = NULL, type = "random", rank.preference = "higher",
                         baseline = "none", baseline.risk = "independent", covariate = NULL, covariate.type = NULL, covariate.model = NULL,
                         mean.d = NULL, prec.d = NULL, mean.Eta = NULL, prec.Eta = NULL, hy.prior.Eta = NULL, mean.bl = NULL, prec.bl = NULL, hy.prior.bl = NULL,
                         mean.cov = NULL, prec.cov = NULL, hy.prior.cov = NULL, hy.prior = NULL, mean.A = NULL, prec.A = NULL) {

  # rename the variables and order them based on specified treatment order
  network <- preprocess.data(Outcomes = Outcomes, Study = Study, Treat = Treat, N = N, SE = SE, response = response, Treat.order = Treat.order, type = type, rank.preference = rank.preference,
                  baseline = baseline, baseline.risk = baseline.risk, covariate = covariate, covariate.type = covariate.type, covariate.model = covariate.model,
                  hy.prior.Eta = hy.prior.Eta, hy.prior.bl = hy.prior.bl, hy.prior.cov = hy.prior.cov, hy.prior = hy.prior, mean.A = mean.A, prec.A = prec.A)
  
  # find characteristic values associated with the network, i.e. number of studies, number of treatments, etc
  characteristics <- find.characteristics(network)
  network <- append(network, characteristics)
  
  # change input dimensions for Outcomes, Study, Treat, and N/SE to fit JAGS coding format
  new.inputs <- change.dimensions(network)
  network <- append(network, new.inputs)

  # generate default priors for heterogeneity if not specified
  hy.prior <- hy.prior.update(network, hy.prior.Eta, hy.prior.bl, hy.prior.cov, hy.prior)
  network <- append(network, hy.prior)
    
  # generate prior data that will be used for running JAGS model
  prior.data <- network.prior.default(network, mean.d, prec.d, mean.Eta, prec.Eta, hy.prior.Eta, mean.bl, prec.bl, hy.prior.bl, mean.cov, prec.cov, hy.prior.cov, hy.prior)
  network$prior.data <- prior.data
  
  # generate JAGS code
  code <- network.rjags(network)
  network$code <- code
  
  # calculate baseline log odds if baseline effect is specified
  mx_bl <- calculate.baseline.log.odds(network)
  network$mx_bl <- mx_bl
  
  # calculate covariate mean if covariate is specified
  cov.mean <- calculate.covariate.mean(network)
  network <- append(network, cov.mean)
  
  class(network) <- "network.data"
  return(network)
}


change.dimensions <- function(network){

  with(network, {
    n <- se <- NULL
    t <- make.byStudy.matrix(data[,"Treat"], data[,"Study"])
    if(response == "binomial" || response == "multinomial"){
      n <- make.byStudy.matrix(data[,"N"], data[,"Study"])
    } else if(response == "normal"){
      se <- make.byStudy.matrix(data[,"SE"], data[,"Study"])
    }
    r <- make.byStudy.Outcome(as.matrix(data[,1:ncat]), data[,"Study"], nstudy, na)
    if(response != "multinomial"){
      r <- r[,,1]
    }
    return(list(t = t, n = n, se = se, r = r))
  })
}


find.characteristics <- function(network){

  with(network,{
    ncat <- dim(data)[2] - 3
    nrow <- dim(data)[1]
    nstudy <- length(Study.order)
    ntreat <- length(Treat.order)
    na <- rle(data[,"Study"])$lengths
   
    ends <- cumsum(na) # End row of trials
    starts <- c(1, ends[-length(ends)] + 1) # Start row of trials
    b.Treat <- rep(NA, length(na))
    b.id <- rep(F, sum(na))
    for (i in 1:length(na)){
      limits <- starts[i]:ends[i]
      b.Treat[i] <- min(data[,"Treat"][limits])
      b.id[limits[b.Treat[i] == data[,"Treat"][limits]]] <- T
    }

    return(list(ncat = ncat, nrow = nrow, nstudy = nstudy, ntreat = ntreat, na = na, b.id = b.id))
  })
}


make.byStudy.Outcome = function(Outcomes, Study, nstudy, na){
  r = structure(.Data = rep(NA, nstudy*max(na)* dim(Outcomes)[2]), .Dim = c(nstudy, max(na), dim(Outcomes)[2]))

  arms_index = NULL
  for(i in 1:length(na)){
    arms_index = c(arms_index, seq(na[i]))
  }
  for(i in 1:dim(Outcomes)[1]){
    r[Study[i],arms_index[i],] = Outcomes[i,]
  }
  return(r)
}


make.byStudy.matrix = function(Treat, Study){

  #make the by-arms vector into a by-study matrix
  nstudy = length(unique(Study))
  na = rle(Study)$lengths
  Study = rep(1:nstudy, times = na)

  t = matrix(NA, nrow = nstudy, ncol = max(na))

  arms_index = NULL
  for(i in 1:length(na)){
    arms_index = c(arms_index, seq(na[i]))
  }

  for(i in 1:length(Treat)){
    t[Study[i], arms_index[i]] = Treat[i]
  }
  return(t)
}


calculate.baseline.log.odds <- function(network){
  
  mx_bl <-
  with(network,{
    if(baseline %in% c("independent", "common", "exchangeable")){
      if(response == "normal"){
        return(mean(r[,1][t[,1]==1], na.rm = TRUE))
      }
      if(response == "binomial"){
        rdummy = r[,1]
        rdummy[r[,1] == 0] = 0.5
        ndummy = n[,1]
        ndummy[r[,1] == 0 & !is.na(r[,1])] = ndummy[r[,1] == 0 & !is.na(r[,1])] + 1
        
        #take only the non-active control group (treatment A) to calculate the observed mean log odds
        p = (rdummy[!is.na(rdummy)]/ndummy[!is.na(rdummy)]) #[t[,1][!is.na(rdummy)]==1]
        lodds = log(p/(1-p))
        lodds <- lodds[is.finite(lodds)]
        return(mean(lodds, na.rm = TRUE))
      }
      if(response == "multinomial"){
        #the first response in the multinomial is the reference, we will call it J
        #take only the control group
        P_J = (r[,1,1]/n[,1])[t[,1]==1]
        P_j = matrix(0, ncol = ncat-1, nrow = sum(t[,1]==1))
        for(j in 2:ncat){
          P_j[,j-1] = (r[,1,j]/n[,1])[t[,1]==1]
        }
        lodds = log(P_j/P_J)
        return(apply(lodds, 2, mean, na.rm = TRUE))
      }
    }
  })
  return(mx_bl)
}

calculate.covariate.mean <- function(network){
  
  with(network,{
    # calculate mean of covariate
    store <- list()
    if(!is.null(covariate)){
      for(i in 1:dim(covariate)[2]){
        nam <- paste("mx",i, sep = "")
        nam <- assign(nam, mean(covariate[,i], na.rm = TRUE))
        
        store[[paste("mx",i, sep = "")]] <- ifelse(covariate.type[i] == "continuous", nam, 0)
        
        nam2 <- paste("x", i, sep = "")
        nam2 <- assign(nam2, covariate[,i])
        store[[paste("x", i, sep = "")]] <- nam2
      }
    }
    return(store)
  })
}
