

preprocess.data <- function(Outcomes = NULL, Study = NULL, Treat = NULL, N = NULL, SE = NULL, response = NULL, Treat.order = NULL, type = "random", rank.preference = "higher",
                            baseline = "none", baseline.risk = "independent", covariate = NULL, covariate.type = NULL, covariate.model = NULL,
                            hy.prior.Eta = NULL, hy.prior.bl = NULL, hy.prior.cov = NULL, hy.prior = NULL, A.probability = NULL){

  # This is first bit of our code that runs when the user enters the data using \code{\link{network.data}}.
  # The function preprocesses the given data to fit the data format necessary to run our NMA model in JAGS.
  # It is essential to know how the data is preprocessed to obtain correct analysis of interest.
  # Thus, we have exported the function so that user can have better understanding on what is going on underneath.
  
  
  network <- list(Outcomes = Outcomes, Study = Study, Treat = Treat, N = N, SE = SE, response = response, Treat.order = Treat.order, type = type, rank.preference = rank.preference,
                  baseline = baseline, baseline.risk = baseline.risk, covariate = covariate, covariate.type = covariate.type, covariate.model = covariate.model,
                  hy.prior.Eta = hy.prior.Eta, hy.prior.bl = hy.prior.bl, hy.prior.cov = hy.prior.cov, hy.prior = hy.prior, A.probability = A.probability)
          
  check.for.errors(network) # check for errors in the data specified

  # changes treatment and study names and order it based on treat.order or alphabetical order if left unspecified
  transform <- transform.data(network)
  data <- transform$data

  # If baseline.risk = "exchangeable", add a fictitious arm with overall reference treatment
  if(baseline.risk == "exchangeable"){
    data <- add.fictitious.row(network, data)
  }
  
  # add in outcomes, study, treat N, SE again
  Outcomes <- as.matrix(data[,1:(ncol(data) - 3)])
  Treat <- data[,"Treat"]
  Study <- data[,"Study"]
  store <- list(Outcomes = Outcomes, Treat = Treat, Study = Study)
  if(response == "binomial" || response == "multinomial"){
    N <- data[,"N"]
  } else if(response == "normal"){
    SE <- data[,"SE"]
  }
  
  # covariate modifications
  if(!is.null(covariate)){
    if(is.null(covariate.model)){
      covariate.model <- "common"
    }
    covariate <- as.matrix(covariate)
    if(is.null(covariate.type)){
      covariate.type = rep("continuous", dim(covariate)[2])
    }
  }
  
  list(Outcomes = Outcomes, Treat = Treat, Study = Study, N = N, SE = SE,
       data = data, Treat.order = transform$Treat.order, Study.order = transform$Study.order,
       response = response, type = type, rank.preference = rank.preference, baseline = baseline, baseline.risk = baseline.risk,
       covariate = covariate, covariate.type = covariate.type, covariate.model = covariate.model)

}



add.fictitious.row <- function(network, data){

  with(network,{
    no_reference <- vector(mode = "integer")
    for(i in seq(length(unique(Study)))){
      if(!1 %in% data[data[,"Study"] == i, "Treat"]){
        no_reference <- c(no_reference, i)
      }
    }

    store <- vector(mode = "integer")
    ncol <- dim(as.matrix(Outcomes))[2]
    if(length(no_reference) != 0){
      for(i in 1:length(no_reference)){
        store <- rbind(store, c(rep(NA, ncol), 1, no_reference[i], 1))
      }
    }

    if(length(store) != 0){
      colnames(store) <- colnames(data)
      data <- rbind(data, store)
    }

    #reorder the data
    ordering <- order(data[,"Study"], data[,"Treat"])
    data <- data[ordering,]
    return(data)
  })
}


transform.data <- function(network){

  with(network, {

    Orig_Study <- Study
    Orig_Treat <- Treat

    #relabel Study names into to a numeric sequence (1:nstudy)
    na <- rle(Study)$lengths
    Study.order <- unique(Study)
    names(Study.order) <- 1:length(Study.order)
    Study <- rep(1:length(unique(Study)), times = na)

    #relabel the treatment according to treatment order specified
    if(is.null(Treat.order)){
      Treat.order <- sort(unique(Treat))
    }
    Treat <- relabel.vec(Treat, Treat.order)
    names(Treat.order) <- 1:length(Treat.order)

    data <- if(response == "normal"){
      cbind(Outcomes, SE, Study, Treat)
    } else {
      cbind(Outcomes, N, Study, Treat)
    }

    ordering <- order(Study, Treat)
    data <- data[ordering,]

    return(list(data = data, Treat.order = Treat.order, Study.order = Study.order))
  })
}


relabel.vec <- function(x, order)
{
  old.x <- x
  x <- rep(NA, length(old.x))
  for (i in seq(length(order))) x[old.x == order[i]] <- i #relabel studies in numerical order starting with one
  return(x)
}


check.for.errors <- function(network){

  with(network, {

    if(is.null(response)){
      stop("Response has to be specified: binomial, multinomial, or normal")
    }
    
    if(is.null(Study) || is.null(Treat) || is.null(Outcomes)){
      stop("Study, Treat, and Outcomes have to be all specified")
    }

    if(response == "multinomial" || response == "binomial"){
      if(is.null(N)){
        stop("If the response is multinomial or binomial, N has to be specified")
      }
    } else if (response == "normal"){
      if(is.null(SE)){
        stop("If the response is normal, SE has to be specified")
      }
    }

    if(is.null(baseline) || !baseline %in% c("none", "independent", "exchangeable", "common")){
      stop("baseline has to be none, independent, exchangeable, or common")
    }

    if(is.null(baseline.risk) || !baseline.risk %in% c("independent", "exchangeable")){
      stop("baseline risk has to be independent or exchangeable")
    }

    if(!type %in% c("fixed", "random")){
      stop("type has to be either fixed or random")
    }

    if(!is.null(covariate)){
      if(!is.matrix(covariate) && !is.vector(covariate)){
        stop("covariate has to be vector if there is only one and it is a matrix if there is more than one")
      }

      if(!all(apply(as.matrix(covariate), 2, is.numeric))){
        stop("covariate has to be numeric type")
      }
    }

    if(!is.null(Treat.order)){
      if(length(unique(Treat.order)) != length(unique(Treat))){
        stop("Need to specify treatment order for all treatments")
      }
      if(!all(Treat.order %in% Treat)){
        stop("Treat.order names have to be specified correctly.")
      }
    }

    na <- rle(Study)$lengths
    if(any(na == 1)) stop("study cannot have only 1 arm or arms have to be next to each other in each study")
    
    if(length(unique(paste(Study, Treat))) != length(Treat)){
      stop("There are repetitive treatments in a single study. Each study should have only one treatment")
    }
    
    
    # check heterogeneity priors
    if(!is.null(hy.prior)){
      check.hy.prior(hy.prior, response)
    }
    if(!is.null(hy.prior.bl)){
      check.hy.prior(hy.prior.bl, response)
    }
    if(!is.null(hy.prior.cov)){
      check.hy.prior(hy.prior.cov, response)
    }
    if(!is.null(hy.prior.Eta)){
      check.hy.prior(hy.prior.Eta, response)
    }
    

  })
}


check.hy.prior <- function(hy.prior, response){
  
  if(length(hy.prior) != 3) stop("length of the hy.prior has to be 3 (Distribution, and two parameters for the distribution)")
  
  distr = hy.prior[[1]]
  stopifnot(distr %in% c("dunif", "dgamma", "dhnorm", "dwish"))
  
  if(response == "normal" || response == "binomial"){
    stopifnot(distr %in% c("dunif", "dgamma", "dhnorm"))
  } else if(response == "multinomial"){
    stopifnot(distr %in% c("dwish"))
  }
}
