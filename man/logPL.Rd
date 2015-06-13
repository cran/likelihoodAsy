\name{logPL}
\alias{logPL}
\title{
Profile likelihood computation
}
\description{
  This function evaluates the profile likelihood for a subset of the model parameter. 
  The result is optionally returned with a minus sign, so the function can be used directly as input 
  to a general-purpose optimizer.  
}
\usage{
logPL(psival, data, thetainit, floglik, fscore=NULL, indpsi, minus=FALSE, 
      trace=FALSE)
}
\arguments{
  \item{psival}{
    A numerical vector containing the value of the parameter of interest.
    }  
	
  \item{data}{
    The data as a list. All the elements required to compute the likelihood function 
    at a given parameter value should be included in this list. The required format of such list 
    will be determined by the user-provided function \code{floglik}.
    }
  \item{thetainit}{
   A numerical vector with the size of the entire model parameter, that will be used as starting 
   point in the constrained optimization performed to obtain the maximum likelihood estimate under 
   the null. The specific meaning of \code{thetainit} is determined by the specification of
    \code{floglik}.  
  }     
  \item{floglik}{
    A function which returns the log likelihood function at a given parameter value. 
    In particular, for a certain parameter value contained in a numerical vector \code{theta},
    a call \code{floglik(theta, data)} should return a scalar numerical value, the log likelihood
    function at \code{theta}. Note that the parameter of interest should be a subset of the coordinates
    of \code{theta}.
    }
 \item{fscore}{
    An optional function which returns the score function at a given parameter value. It must return a numerical
    vector of the same length as \code{thetainit}. For a certain parameter value contained in a numerical 
    vector \code{theta}, a call \code{fscore(theta, data)} should return the gradient of the log likelihood function 
    at \code{theta}. Default is \code{NULL}, implying that numerical differentiation will be employed.
    }
\item{indpsi}{
    A vector of integers in the range \code{1:length(theta)} containing the indexes of the parameter of interest, so
    that the parameter of interest will be given by \code{theta[indpsi]}. 
   }
\item{minus}{
    Logical. Should the profile likelihood be multiplied by -1? This may be useful for usage with optimizers. 
    Default is \code{FALSE}. 
   }   
\item{trace}{
     Logical. When set to \code{TRUE} will cause the printing of the MPL value, which can be useful to monitor 
     optimization. Default is \code{FALSE}.
    }   
}
 \value{
  A scalar value, minus the profile likelihood at \code{psival}.
}
\details{
  This function is designed to be used with external functions, such as optimizers and evaluators over a grid of points. 
  } 
\references{
Severini, T.A. (2000). Likelihood Methods in Statistics. Oxford University Press.
}
\examples{
# A negative binomial example, taken from Venables and Ripley (2002, MASS4 book)
library(MASS)
# The quine data are analysed in Section 7.4
data(quine)
# We fit a model with just the main effects
quine.nb1 <- glm.nb(Days ~ Eth + Sex + Age + Lrn, data = quine) 
# The data list includes the design matrix and the response vector
quinedata<-list(X=model.matrix(quine.nb1), y=quine$Days)      
# Let us define the various functions
# Log likelihood, log link
logLikNbin <- function(theta,data) 
{
  y <- data$y
  X <- data$X
  eta <- X \%*\% theta[1:ncol(X)] 
  mu <- exp(eta)
  alpha <- theta[ncol(X)+1]
  l <- sum(lgamma(y + alpha) + y * log(mu) - (alpha + y) * log(alpha + mu) 
            - lgamma(alpha) + alpha * log(alpha))
  return(l)
}


# Score function
gradLikNbin <- function(theta,data) 
{
  y <- data$y
  X <- data$X
  eta <- X \%*\% theta[1:ncol(X)] 
  mu <- exp(eta)
  alpha <- theta[ncol(X)+1]
  g <-rep(0,ncol(X)+1)
  g[1:ncol(X)] <- t(y - (alpha+y)*mu / (alpha+mu)) \%*\% X
  g[ncol(X)+1] <- sum(digamma( y + alpha) - log(alpha + mu) - (alpha + y) / (alpha + mu) 
                  - digamma(alpha) + 1 + log(alpha))
  return(g)
}
# Data generator
genDataNbin<- function(theta,data)
{
  out <- data
  X <- data$X
  eta<- X \%*\% theta[1:ncol(X)] 
  mu <- exp(eta)
  out$y <- rnegbin(length(data$y), mu=mu, theta=theta[ncol(X)+1])
  return(out)
}		
# First we refine the maximum likelihood estimates
mleFull <- optim( c(coef(quine.nb1),quine.nb1$theta), logLikNbin, gr=gradLikNbin,
           method="BFGS", data=quinedata, control=list(fnscale=-1), hessian=TRUE) 
 # Then we can plot the profile likelihood
list.psi <- seq(0.90, 1.70, l=30)
list.prof <- sapply(list.psi, logPL, data=quinedata, thetainit=mleFull$par, floglik=logLikNbin, 
                    fscore=gradLikNbin, indpsi=8, trace=FALSE) 
plot(list.psi, list.prof-max(list.prof), type="l", xlab=expression(psi), ylab="Log likelihood")
}
\keyword{htest}
