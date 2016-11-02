\name{logMPL}
\alias{logMPL}
\title{
Modified profile likelihood computation
}
\description{
  This function evaluates the Modified Profile Likelihood (MPL) for a subset of the 
  model parameter. The result is optionally returned with a minus sign, so the function 
  can be used directly as input to a general-purpose optimizer.  
}
\usage{
logMPL(psival, data, mle, floglik, fscore=NULL, indpsi, datagen, R=500, seed=NULL, 
       minus=FALSE, onestep=FALSE, jhat=NULL, trace=FALSE)
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
  \item{mle}{
   A numerical vector, containing the maximum likelihood estimate of the entire model parameter. 	
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
    vector of the same length of \code{mle}. For a certain parameter value contained in a numerical vector \code{theta},
    a call \code{fscore(theta, data)} should return the gradient of the log likelihood function at \code{theta}. Default 
    is \code{NULL}, implying that numerical differentiation will be employed.
    }
\item{indpsi}{
    A vector of integers in the range \code{1:length(theta)} containing the indexes of the parameter of interest, so
    that the parameter of interest will be given by \code{theta[indpsi]}. 
   }
\item{datagen}{
    A function which simulates a data set. A call  \code{datagen(theta, data)} will generate a copy of the \code{data} list, 
     with the values of the response variable replaced by a set of values simulated from the parametric statistical model assumed 
     for the response variable. 
    }
\item{R}{
    The number of Monte Carlo replicates used for computing the modified profile likelihood. A positive integer, default 
    is \code{500}.
    }
\item{seed}{
    Optional positive integer, the random seed for the Monte Carlo computation. Default is \code{NULL}.
    }  
\item{minus}{
    Logical. Should the modified profile likelihood be multiplied by -1? This may be useful for usage with 
    optimizers.  Default is \code{FALSE}. 
   }   
\item{onestep}{
    Logical. If set to \code{TRUE} the constrained estimate of the nuisance parameter is replaced by a one-step approximation around the 
    maximum likelihood estimate. Default is \code{FALSE}. 
   }      
\item{jhat}{
    A squared matrix with dimension equal to \code{length(mle)} containing the observed information matrix evaluated at \code{mle}.
    It is employed only when \code{onestep=TRUE}. Default is \code{NULL}. 
   }    
\item{trace}{
     Logical. When set to \code{TRUE} will cause the printing of the MPL value, which can be useful to monitor 
     optimization. Default is \code{FALSE}.
    }   
}
 \value{
  A scalar value, minus the modified profile likelihood at \code{psival}.
}
\details{
  The function  implements the Modified Profile Likelihood employing the approximation to sample space
  derivatives proposed in Skovgaard (1996). The function 
  is designed to be used with external functions, such as optimizers and evaluators over a grid of points.  
}

\references{
Severini, T.A. (2000). Likelihood Methods in Statistics. Oxford University Press.

Skovgaard, I.M. (1996) An explicit large-deviation approximation to one-parameter tests. 
\emph{Bernoulli}, \bold{2}, 145--165.

}
\examples{
# Approximating the conditional likelihood for logistic regression
# Let us define the various functions	
# Log likelihood for logistic regression
loglik.logit<- function(theta, data) 
{
  y <- data$y
  den <- data$den
  X <- data$X
  eta <- X \%*\% theta
  p <- plogis(eta)
  l <- sum(y * log(p) + (den - y) * log(1-p))
  return(l)
}
# Score function
grad.logit<- function(theta, data) 
{
  y <- data$y
  den <- data$den
  X <- data$X
  eta <- X \%*\% theta
  p <- plogis(eta)
  out <- t(y - p * den) \%*\% X
  return(drop(out))
}
# Data generator
gendat.logit<- function(theta, data)
{
  X <- data$X
  eta <- X \%*\% theta
  p <- plogis(eta)
  out <- data
  out$y <- rbinom(length(data$y), size = data$den, prob = p)
  return(out) 
}		
# Famous crying babies data
data(babies)	
mod.glm <- glm(formula = cbind(r1, r2) ~ day + lull - 1, family = binomial, 
               data = babies)
data.obj <- list(y = babies$r1, den = babies$r1 + babies$r2, 
                 X = model.matrix(mod.glm))	
# Numerical optimization of profile and modified profile log likelihoods
max.prof <- nlminb(0, logPL, data=data.obj, thetainit=coef(mod.glm), 
                  floglik=loglik.logit, fscore=grad.logit, indpsi=19, minus=TRUE, trace=FALSE)
max.mpl <- nlminb(0, logMPL, data=data.obj, mle=coef(mod.glm), 
                  floglik=loglik.logit, fscore=grad.logit, datagen=gendat.logit,
                  indpsi=19, R=50, seed=2020, minus=TRUE, trace=FALSE)
c(max.prof$par, max.mpl$par)                 
# We can plot the profile likelihood and the modified profile likelihood
# R=50 suffices for the modified profile likelihood as the model is a full exp. family
psi.vals <- seq(-0.3, 3.7, l=20)
obj.prof <- sapply(psi.vals, logPL, data=data.obj, thetainit=coef(mod.glm), 
                floglik=loglik.logit, fscore=grad.logit, indpsi=19, trace=FALSE)
obj.mpl <- sapply(psi.vals, logMPL, data=data.obj, mle=coef(mod.glm), 
                floglik=loglik.logit, fscore=grad.logit, datagen=gendat.logit,
                indpsi=19, trace=FALSE, R=50, seed=2020)
par(pch="s")
plot(psi.vals, obj.prof - max(obj.prof), type="l", xlab=expression(psi), 
     ylab="log likelihood", lwd=2, las=1)
lines(psi.vals, obj.mpl - max(obj.mpl), col="red", lwd=2)
legend("topright", col=c(1, 2), lty=1, lwd=2, legend=c("Profile","MPL"), bty="n")
  
}
\keyword{htest}
