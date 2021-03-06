#' bnma: A package for network meta analysis using Bayesian methods
#'
#' A package for running Bayesian network meta analysis
#'
#' Network meta-analysis or mixed treatment comparison (MTC) is a method that allows simultaneous comparison of more than two treatments.
#' We use a Bayesian approach to combine both direct and indirect evidence as in Dias et al. 2013a.
#' This package is a user friendly application that can run network meta analysis models without having to code a JAGS model.
#' The program takes the input data and transforms it to a suitable format of analysis, generates a JAGS model and reasonable
#' initial values and runs the model through the rjags package.
#' The focus of this package was inclusion of multinomial response and various options for adding covariates and/or baseline risks effects.
#' Also, while sampling, the package uses Gelman-Rubin convergence criteria to decide whether to continue sampling or not.
#' Furthermore, package includes different models such as contrast based models and unrelated mean effects (UME) model and nodesplitting model to test for inconsistency.
#'
#' @docType package
#' @name bnma-package
#' @references A.J. Franchini, S. Dias, A.E. Ades, J.P. Jansen, N.J. Welton (2012), \emph{Accounting for correlation in network meta-analysis with multi-arm trials}, Research Synthesis Methods 3(2):142-160. \doi{10.1002/jrsm.1049}
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. \doi{10.1177/0272989X12458724}
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013b), \emph{Heterogeneity-Subgroups, Meta-Regression, Bias, and Bias-Adjustment}, Medical Decision Making 33(5):618-640. \doi{10.1177/0272989X13485157}
#' @references S. Dias, N.J. Welton, D.M. Caldwellb, A.E. Ades (2010), \emph{Checking consistency in mixed treatment}, Statistics in Medicine 29(7-8, Sp. Iss. SI): 932-944. \doi{10.1002/sim.3767}
#' @references S. Dias, N.J. Welton, A.J. Sutton, D.M. Caldwell, G. Lu, and A.E. Ades (2013), \emph{Evidence synthesis for decision making 4: inconsistency in networks of evidence based on randomized controlled trials}, Medical Decision Making 33(5):641-656. \doi{10.1177/0272989X12455847}
#' @references C.H. Schmid, T.A. Trikalinos, I. Olkin (2014), \emph{Bayesian network meta-analysis for unordered categorical outcomes with incomplete data}, Research Synthesis Methods 5(2):162-185. \doi{10.1002/jrsm.1103}
#' @references A. Gelman, D.B. Rubin (1992), \emph{Inference from iterative simulation using multiple sequences}, Statistical Science 7(4):457-472. \doi{10.1214/ss/1177011136}
#' @references D.J. Spiegelhalter, N.G. Best, and B.P. Carlin (1998), \emph{Bayesian deviance, the effective nunmber of parameters, and the comparison of arbitrarily complex models}, Technical report, MRC Biostatistics Unit, Cambridge, UK.
#' @references F.A. Achana, N.J. Cooper, S. Dias, G. Lu, S.J.C. Rice, D. Kendrick, A.J. Sutton (2012), \emph{Extending methods for investigating the relationship between treatment effect and baseline risk from pairwise meta-analysis to network meta-analysis}, Statistics in Medicine 32(5):752-771. \doi{10.1002/sim.5539}
#' @references F.A. Achana, N.J. Cooper, S. Bujkiewicz, S.J. Hubbard, D. Kendrick, D.R. Jones, A.J. Sutton (2014), \emph{Network meta-analysis of multiple outcomes measures accounting for borrowing of information across outcomes}, BMC Medical Research Methodology 14:92. \doi{10.1186/1471-2288-14-92}
#' @references G. Salanti, A.E. Ades, J.P.A. Ioannidisa (2011), \emph{Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial}, Journal of Clinical Epidemiology 64(2):163-171. \doi{10.1016/j.jclinepi.2010.03.016}
#' @references G. van Valkenhoef, G. Lu, B. de Brock, H. Hillege, A.E. Ades, and N.J. Welton (2012), \emph{Automating network meta-analysis}, Research Synthesis Methods 3(4):285-299. \doi{10.1002/jrsm.1054}
#' @references N.J. Cooper, A.J. Sutton, D. Morris, A.E. Ades, N.J. Welton (2009), \emph{Addressing between-study heterogeneity and inconsistency in mixed treatment comparisons: Application to stroke prevention treatments in individuals with non-rheumatic atrial fibrillation}, Statistics in Medicine 28:1861-1881. \doi{10.1002/sim.3594}
#' @references W. Viechtbauer (2010), \emph{Conducting meta-analyses in R with the metafor package}, Journal of Statistical Software, 36(3):1-48. \doi{10.18637/jss.v036.i03}
#' @seealso \code{\link{network.data}}, \code{\link{network.run}}
NULL

#' Beta blockers to prevent mortality after myocardial infarction
#'
#' A dataset of 22 trials investigating beta blockers versus control to prevent mortality after myocardial
#' infarction. Control is coded as 1 and beta blocker treatment is coded as 2.
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. \doi{10.1177/0272989X12458724}
#' @format A list of Outcomes, Treat, Study, and N.
"blocker"

#' Trials of certolizumab pegol (CZP) for the treatment of rheumatoid arthritis in patients
#'
#' A dataset of 12 trials for investigating CZP for the treatment for those who had failed on disease-modifying antirheumatic drugs, including methotrexate (MTX).
#' Data provides the number of patients who have improved and there are 6 different treatments with placebo. Mean disease duration (years) is provided as a covariate.
#'
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013b), \emph{Heterogeneity-Subgroups, Meta-Regression, Bias, and Bias-Adjustment}, Medical Decision Making 33(5):618-640. \doi{10.1177/0272989X13485157}
#' @format A list of Outcomes, Treat, Study, N, covariate, and Treat.order
"certolizumab"

#' Trials of statins for cholesterol lowering vs. placebo or usual care
#'
#' A dataset of 19 trials of statins for cholesterol lowering vs. placebo.
#' Each trial has a subgroup indicator for primary prevention (patients included had no previous heart disease) or
#' secondary prevention (patients had previous heart disease). Dummy variable is coded such that covariate is equal to 1
#' if a study is a secondary prevention study and 0 if a study is a primary prevention study.
#'
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013b), \emph{Heterogeneity-Subgroups, Meta-Regression, Bias, and Bias-Adjustment}, Medical Decision Making 33(5):618-640. \doi{10.1177/0272989X13485157}
#' @format A list of Outcomes, Treat, Study, N, covariate, and Treat.order
"statins"

#' Dopamine agonists as adjunct therapy in Parkinson's disease
#'
#' A dataset of 7 studies investigating the mean lost work-time reduction in patients given
#' 4 dopamine agonists and placebo as adjunct therapy for Parkinson's disease.
#' There is placebo and four active drugs coded 2 to 5.
#'
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. \doi{10.1177/0272989X12458724}
#' @format A list of Outcomes, Treat, Study, N, covariate, and Treat.order
"parkinsons"


#' Dopamine agonists as adjunct therapy in Parkinson's disease
#'
#' A contrast level (i.e. treatment difference) dataset of 7 studies investigating the mean lost work-time reduction in patients given
#' 4 dopamine agonists and placebo as adjunct therapy for Parkinson's disease. Placebo is coded as 1, and four active drugs are coded 2 to 5.
#' There is placebo, coded as 1, and four active drugs coded 2 to 5.
#'
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. \doi{10.1177/0272989X12458724}
#' @format A list of Outcomes, Treat, SE, na, and V
"parkinsons_contrast"


#' Trials of low dose and high dose statins for cardiovascular disease vs. placebo
#'
#' A dataset of 17 studies investigating dosage of statin for cardiovascular disease.
#' There are two treatments and a placebo. High dose statin is coded as 3, low dose statin as 2, and placebo is coded as 1 and treated as a baseline treatment.
#' Outcomes are reported as three mutually exclusive unordered outcomes.
#' First column of the outcome is the patients who are still alive (ALIVE). Second column is fatal non-cardiovascular disease (FnCVD).
#' And, the last column is fatal cardiovascular disease (FCVD).
#'
#' @references C.H. Schmid, T.A. Trikalinos, I. Olkin (2014), \emph{Bayesian network meta-analysis for unordered categorical outcomes with incomplete data}, Research Synthesis Methods 5(2):162-185. \doi{10.1002/jrsm.1103}
#' @format A list of Outcomes, Treat, Study, and N
"cardiovascular"


#' Smoking cessation counseling programs
#'
#' Twenty-four studies, including 2 three-arm trials, compared 4 smoking cessation counseling programs and recorded the number of
#' individuals with successful smoking cessation at 6 to 12 month. Counseling programs include 1 = no intervention, 2 = self-help,
#' 3 = individual counseling, and 4 = group counseling. 
#'
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. \doi{10.1177/0272989X12458724}
#' @format A list of Outcomes, Treat, Study, and N
"smoking"


#' Thrombolytic drugs and percutaneous transluminal coronary angioplasty
#'
#' A dataset consisting of 50 trials comparing 8 thrombolytic drugs and percutaneous transluminal coronary angioplasty, following
#' acute myocardial infarction. Data consist of the number of deaths in 30 or 35 days and the number of patients in each treatment arm.
#' There are 9 treatments in total: streptokinase (1), alteplase (2), accelerated alteplase (3), streptokinase + alteplase (4),
#' reteplase (5), tenecteplase (6), percutaneous transluminal coronary angioplasty (7), urokinase (8), anistreptilase (9)
#' 
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. \doi{10.1177/0272989X12458724}
#' @format A list of Outcomes, Treat, Study, and N
"thrombolytic"


#' @importFrom graphics axis legend lines mtext plot points title abline curve text
NULL

#' @importFrom stats coef end lm quantile rchisq rnorm sd start window aggregate
NULL

#' @importFrom utils combn
NULL

#' @import coda
NULL

#' @import ggplot2
NULL

#' @import grid
NULL

