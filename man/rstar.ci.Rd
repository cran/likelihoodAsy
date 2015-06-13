\name{rstar.ci}
\alias{rstar.ci}
\title{
Confidence intervals on a scalar function of interest by the r* statistic
}
\description{
  This function obtains confidence intervals for a scalar function of interest, based on the r* statistic.
}
\usage{
	rstar.ci(data, thetainit, floglik, fscore=NULL, fpsi, datagen, R=1000, seed=NULL, 
	         ronly=FALSE, psidesc=NULL, constr.opt="solnp", lower=NULL, upper=NULL, 
	         control=list(...), ...)
}
\arguments{
  \item{data}{
    The data as a list. All the elements required to compute the likelihood function 
    at a given parameter value should be included in this list.
    }
  \item{thetainit}{
    A numerical vector containing the initial value for the parameter of the model. It will be used 
    as starting point in the numerical optimization of the likelihood function.
    }
  \item{floglik}{
    A function which returns the log likelihood function at a given parameter value. 
    In particular, for a certain parameter value contained in a numerical vector \code{theta},
    a call \code{floglik(theta, data)} should return a scalar numerical value, the log likelihood
    function at \code{theta}.
    }
 \item{fscore}{
    An optional function which returns the score function at a given parameter value. It must return a 
    numerical vector of the same length of \code{thetainit}. For a certain parameter value contained in a 
    numerical vector \code{theta}, a call \code{fscore(theta, data)} should return the gradient of the 
    log likelihood function at \code{theta}. Default is \code{NULL}, implying that numerical differentiation 
    will be employed.
    }
\item{fpsi}{
    A function which specifies the parameter of interest. A call \code{fpsi(theta)} should return a scalar 
    value. 
    }
\item{datagen}{
    A function which simulates a data set. A call \code{datagen(theta, data)} will generate a copy of the 
    \code{data} list, with the values of the response variable replaced by a set of values simulated from  
    the parametric statistical model assumed for the response variable.
    }
\item{R}{
    The number of Monte Carlo replicates used for computing the r* statistic. A positive integer, default is 
    \code{1000}.
    }
\item{seed}{
    Optional positive integer, the random seed for the Monte Carlo computation. Default is \code{NULL}.
    }  
\item{ronly}{
     Logical. If set to \code{TRUE} the computation of the r* statistic will be skipped, and only the value 
     of the signed likelihood ratio test statistic r will be returned by the procedure, without any Monte 
     Carlo computation. Default is \code{FALSE}.
    }          
 \item{psidesc}{
    An optional character string describing the nature of the parameter of interest. Default is \code{NULL}.
    }
  \item{constr.opt}{  
    Constrained optimizer used for maximizing the log likelihood function under the null hypothesis. Possible
    values are   \code{"solnp"} or \code{"alabama"}, with the former employing  the \code{solnp} function from
    package \code{Rsolnp} and the latter the \code{constrOptim.nl} from the package \code{alabama}. Defauls is 
    \code{"solnp"}. 
  }    
 \item{lower, upper}{
    Optional numeric values defining the lower/upper limit of a grid of points for the parameter of interest, 
    where the r* statistic will be evaluated. Default is \code{NULL}.
    }
 \item{control}{   
    A list of parameters for controlling the computation of confidence intervals. See \link{rstar.ci.control}.
   }   
 \item{...}{
 	Arguments to be used to form the default \code{control} argument if it is not supplied directly.
 	}   
}
 \value{
  The returned value is an object of class \code{"rstarci"}, containing the following components:

  \item{psivals}{
    A list of values for the parameter of interest for which the r and r* statistics have been evaluated.
  }
  \item{rvals}{
    A numerical list containing the values of the r statistic evaluated at each element of \code{psivals}.
  }
  \item{NPvals, INFvals}{
    Numerical lists containing the values of the Nuisance Parameter adjustment (NP) and the Information 
    adjustment (INF) from the decomposition of the r*-r adjustment, for each of the \code{psivals} values.  
    Not computed when \code{ronly = TRUE}.
  }
  \item{rsvals}{
    The observed value of the r* statistic at each element of \code{psivals}.  Not computed when \code{ronly = TRUE}.
  }
  \item{CIr}{
    A 3 x 2 matrix containing the 90\%, 95\% and 99\% confidence intervals for the parameter of interest (first, 
    second and third row respectively) based on the first-order r statistic.
  }
  \item{CIrs}{
    A 3 x 2 matrix containing the 90\%, 95\% and 99\% confidence intervals for the parameter of interest (first, 
    second and third row respectively) absed on the r* statistic. Not computed when \code{ronly = TRUE}.
    }
  \item{seed}{
    Random seed used for Monte Carlo replicates used for computing the r* statistic. Not returned when \code{ronly = TRUE}.
  }
  \item{psidesc}{
   A character string describing the nature of the parameter of interest. 
  }
  \item{R}{
    Number of Monte Carlo replicates used for computing the r* statistic.  Not returned when \code{ronly = TRUE}.
  }
  
  There are \code{print}, \code{summary} and \code{plot} methods for this class. 
}
\details{
  The function obtains 90\%, 95\% and 99\%  two-sided confidence intervals for the scalar function of 
  interest based on the r* statistic.
    
  The function requires the user to provide three functions defining the log likelihood function,
  the scalar parametric function of interest, and a function for generating
  a data set from the assumed statistical model. A further function returning the gradient of the 
  log likelihood is not required, but if provided it will speed up the computation. 
  
  When \code{ronly = TRUE} the function literally returns the value of the signed likelihood ratio test
  statistic r only. The function handles also one-parameter models.

  The function provides two different strategies to obtain the various confidence intervals. The default
  strategy, invoked by leaving either \code{lower} or \code{upper} to \code{NULL}, starts from the MLE 
  and moves away in a stepwise fashion, until the r* statistic crosses the standard normal quantiles 
  corresponding to the 99\% two-sided confidence interval. It is crucial to start the search a bit away 
  from the MLE, where the r* is singular, and this is regulated by the \code{away} argument of the 
  \link{rstar.ci.control} function. The first strategy may fail to cross the target normal quantiles when 
  the profile likelihood has an upper asymptote. For such cases, and for any other instances when the 
  output of the default strategy is deemed not satisfactory, it is possible to specify the 
  range of a grid of values where the r* statistic will be evaluated. The \code{lower} and \code{upper} 
  argument specify the lower and upper limit of such grid, whereas the number of points is controlled by the 
  \code{npoints} of the \link{rstar.ci.control} function.  

}

\references{
The method implemented in this function was proposed in

Skovgaard, I.M. (1996). An explicit large-deviation approximation to one-parameter tests. 
\emph{Bernoulli}, \bold{2}, 145--165.

For a general review

Severini, T.A. (2000). Likelihood Methods in Statistics. Oxford University Press.
}
\seealso{
  \code{\link{rstar}, \link{rstar.ci.control}}.
}
\examples{
# A negative binomial example, taken from Venables and Ripley (2002, MASS4 book)
library(MASS)
# The quine data are analysed in Section 7.4
data(quine)
# We fit a model with just the main effects
quine.nb1 <- glm.nb(Days ~ Eth + Sex + Age + Lrn, data = quine) 
# The data list includes the design matrix and the response vector
quinedata <- list(X=model.matrix(quine.nb1), y=quine$Days)      
# Let us define the required functions
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
genDataNbin <- function(theta,data)
{
  out <- data
  X <- data$X
  eta<- X \%*\% theta[1:ncol(X)] 
  mu <- exp(eta)
  out$y <- rnegbin(length(data$y), mu=mu, theta=theta[ncol(X)+1])
  return(out)
}	
# Confidence intervals for the coefficient of EthN 
\dontrun{
obj <- rstar.ci(quinedata, thetainit=c(coef(quine.nb1),quine.nb1$theta),  floglik=logLikNbin, 
                datagen=genDataNbin, fscore=gradLikNbin, fpsi=function(theta) theta[2], R=1000, 
                psidesc="Coefficient of EthN")
print(obj)
summary(obj)
plot(obj)
# Confidence intervals for the overdispersion parameter
obj <- rstar.ci(quinedata, thetainit=c(coef(quine.nb1),quine.nb1$theta),  floglik=logLikNbin, 
                datagen=genDataNbin, fscore=gradLikNbin, fpsi=function(theta) theta[8], R=1000,
                psidesc="Overdispersion parameter")
summary(obj)
plot(obj)
}
}
\keyword{htest}
