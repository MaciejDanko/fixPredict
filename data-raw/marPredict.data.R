library(marPredict)
library(lme4)

################################################################################
#  A function to simulate random intercept data
################################################################################

#' Simple simulations for Poisson or binomial random intercept models 
#' (experimental function!)
#' 
#' @description 
#' This is an experimental function, for more complicated models please use the 
#' generic \code{\link{simulate}} function or call directly the 
#' \code{lme4::simulate.formula} function.
#' 
#' @param formula a formula used to generate the data; 
#' only fixed effect factors should be specified here.
#' @param theta a theta parameter for random intercept.
#' @param coef a vector of factor coefficients for resulting model matrix.
#' @param coef.c a value for coefficient for a continouse variable. 
#' @param func.c a function that describes contribution of the continouse 
#' variable to the link predictor. 
#' It takes two arguments: \code{C1} - a vector of values of continouse 
#' variable and \code{coef.c} - the parameter related to this variable.
#' @param s,vcovmat a residual standard deviation and a variance-covariance 
#' matrix for the fixed model coefficients.
#' @param n.levels a vector specifying a number of levels for each variable 
#' in the formula.
#' @param n.ID a number of levels (IDs) for random intercept.
#' @param n.pop a resulting sample size.
#' @param family a distribution family; either 'Poisson' or 'binomial'.
#' @param link a specification for the model link function. 
#' See \code{\link{family}{stats}}.
#' @export
#' @seealso fixPredict
#' @return a data frame with simulated data set
#' @import stats lme4 MASS
#' @author Maciej J. Danko
sim_glmer_data <- function(formula = ~ X1 * X2, theta=0.75,
                           coef = c(5, 0.5, 0.7, -0.5, -0.7, 0.2, 0.3, 0.3, 0.2),
                           coef.c = NULL,
                           func.c = function(C1, coef.c) C1*coef.c,
                           s=0, vcovmat=diag(rep(s*s,length(coef))), 
                           n.levels=c(3,3), n.ID = 100, n.pop=1e4,
                           family=c('Poisson','binomial'), link){
  
  myletters<-function(x, n) 
    sapply(x, function(k) paste(rep(letters[k],n),sep='',collapse=''))
  
  formula <- stats::delete.response(formula)
  terms.nam <- attr(stats::terms(formula),'term.labels')
  if (any(terms.nam %in% c('ID','Y','linkY','IDv','C1'))) 
    stop(
      '"ID", "Y", "linkY", "IDv", and "C1" cannot be used as names for factors',
         call.=FALSE)
  rnd.term <- grepl('|',terms.nam, fixed=TRUE)
  family <- tolower(family[1])
  if (any(rnd.term)) 
    warning('Please remove random terms from the formula.',call.=FALSE)
  rnd.term.nam <- terms.nam[which(rnd.term)]
  terms.nam <- terms.nam[which(!rnd.term)]
  variables.nam <-unique(unlist(
    regmatches(terms.nam,gregexpr(':', terms.nam, 
                                  fixed = TRUE),invert = TRUE)))
  if (length(n.levels)!=length(variables.nam))
    stop(paste('n.levels should be a vector of length',
               length(variables.nam)), call.=FALSE)
  if (length(rnd.term.nam)) formula <- 
    update(formula, paste('~. - (',rnd.term.nam,')'))
  # generate data - fixed effects
  data <- lapply(seq_along(n.levels), function(k)
    if (n.levels[k] > 1){ 
      sample(myletters(seq_len(n.levels[k]),k), n.pop, replace=TRUE)
    } else {
      stats::runif(n.pop)
    })
  names(data)<-variables.nam
  data <- as.data.frame(data)
  
  # buildmodel matrix
  fixed.mm<-stats::model.matrix(formula, data)
  if(ncol(fixed.mm)!=length(coef)){
    stop(paste('coefi should be a vector of',ncol(fixed.mm),'elements:',
               paste(colnames(fixed.mm),collapse = ', ')))
  }
  names(coef)<-colnames(fixed.mm)
  
  # generate data - random effects
  IDnames <- paste('ID',seq_len(n.ID),sep='')
  ID <- sample(IDnames, n.pop, replace=TRUE)
  rnd_eff<-stats::rnorm(length(IDnames), 0, theta)
  names(rnd_eff) <- IDnames
  rnd_eff <- rnd_eff - mean(rnd_eff)
  IDv <- rnd_eff[ID]
  
  # linear link predictor
  if (length(vcovmat) && sum(diag(vcovmat))) {
    coefMat <- 
      mvrnorm(n.pop, Sigma=as.matrix(vcovmat), mu =coef)
    linkY <- rowSums(coefMat * fixed.mm) + IDv
  } else
    linkY <-  fixed.mm%*%coef + IDv 
  
  if (length(coef.c)){
    if (length(coef.c)>1) {
      coef.c <- coef.c[1]
      warning('currently only one continouse variablecan be simulated.', 
              call. = FALSE)
    }
    if(!length(func.c)) stop('c.func must bespeciffied',call.=FALSE)
    C1 <- runif(n.pop)-0.5
    linkY <- linkY + func.c(C1, coef.c)
    data$C1 <- C1
  } 
  
  if (family=='poisson') {
    if (missing(link)) {
      Y <- stats::rpois(n=rep(1,length(linkY)),
                        lambda=poisson()$linkinv(linkY))
    } else 
      Y <- stats::rpois(n=rep(1,length(linkY)),
                        lambda=poisson(link)$linkinv(linkY))
  } else if (family=='binomial') {
    if (missing(link)){
      Y <- stats::rbinom(n=rep(1,length(linkY)),size=1, 
                         prob=binomial()$linkinv(linkY))
    } else
      Y <- stats::rbinom(n=rep(1,length(linkY)),size=1, 
                         prob=binomial(link)$linkinv(linkY))
  } else stop('Not supported family of distributions.',call.=FALSE)
  data <- data.frame(Y=Y, linkY=linkY, data, ID, IDv)[order(ID), ]
  attr(data,'coef') <- coef
  attr(data,'theta') <- theta
  rownames(data)<-NULL
  return(data)
}

################################################################################
#  Artificial random intercept data set # 1
################################################################################

set.seed(5)
data1 <- sim_glmer_data(formula = ~ X1 * X2, theta = 0.75,
                        coef= c(5, 0.5, 0.7, 0.5, -0.7, 0.2, 0.3, 0.3, 0.2),
                        n.levels=c(3,3), n.ID = 100, n.pop=1e3)
head(data1)

# fit the model
fit1 <- lme4::glmer(Y~X1*X2+(1|ID),family=poisson(),data=data1)
blmeco::dispersion_glmer(fit1)

# comapre the data-generating coeficients and the fitted coeficients
rbind(cbind(data=attr(data1,'coef'),
            fit=round(fixef(fit1),2)),
      theta=round(
        c(data=attr(data1,'theta'),
          fit=getME(fit1,'theta')),2))

# The model can precisely recapture the coefficients used to 
# generate the data. 
marPredict_e1 <- data1 
usethis::use_data(marPredict_e1, overwrite = TRUE)

################################################################################
#  Artificial random intercept data set # 2
################################################################################
set.seed(5)
data2 <- sim_glmer_data(formula = ~ X1, theta = 0.75,
                        coef= c(3, 0.5, -0.7),
                        coef.c = 5,
                        func.c = function(C1, coef.c) 
                          abs(C1)^1.5 * coef.c - 1.5,
                        n.levels=3, n.ID = 100, n.pop=1e3)

marPredict_e2 <- data2 
usethis::use_data(marPredict_e2, overwrite = TRUE) 

################################################################################
#  Artificial random intercept data set # 3
################################################################################
set.seed(5)
data3 <- sim_glmer_data(formula = ~ X1, theta = 0.9,
                        coef= c(0.5, 1.5, -1.25),
                        coef.c = 2,
                        func.c = function(C1, coef.c) C1 * coef.c,
                        n.levels=3, n.ID = 100, n.pop=1e3,
                        family='binomial')

marPredict_e3 <- data3 
usethis::use_data(marPredict_e3, overwrite = TRUE)
