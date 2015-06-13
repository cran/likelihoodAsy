\name{rstar}
\alias{rstar}
\title{
Inference on a scalar function of interest by the r* statistic
}
\description{
  This function evaluates the  r* statistic for testing of a scalar function of interest.
}
\usage{
rstar(data, thetainit, floglik, fscore=NULL, fpsi, psival, datagen, R=1000, seed=NULL, 
      trace=TRUE, ronly=FALSE,  psidesc=NULL, constr.opt="solnp")
}
\arguments{
  \item{data}{
    The data as a list. All the elements required to compute the log likelihood function 
    at a given parameter value should be included in this list.
    }
  \item{thetainit}{
    A numerical vector containing the initial value for the parameter of the model. It will be used 
    as starting point in the numerical optimization of the log likelihood function.
    }
  \item{floglik}{
    A function which returns the log likelihood function at a given parameter value. 
    In particular, for a certain parameter value contained in a numerical vector \code{theta},
    a call \code{floglik(theta, data)} should return a scalar numerical value, the log likelihood
    function at \code{theta}.
    }
 \item{fscore}{
    An optional function which returns the score function at a given parameter value. It must return a numerical
    vector of the same length of \code{thetainit}. For a certain parameter value contained in a numerical vector \code{theta},
    a call \code{fscore(theta, data)} should return the gradient of the log likelihood function at \code{theta}. Default 
    is \code{NULL}, implying that numerical differentiation will be employed.
    }
\item{fpsi}{
    A function which specifies the parameter of interest. A call \code{fpsi(theta)} should return a scalar value. 
    }
\item{psival}{
    A numerical scalar value containing the value of the parameter of interest under testing.
    }  
\item{datagen}{
    A function which simulates a data set. A call \code{datagen(theta, data)} will generate a copy of the \code{data} list, with the  
    values of the response variable replaced by a set of values simulated from the parametric statistical model assumed for the response 	
    variable.
    }
\item{R}{
    The number of Monte Carlo replicates used for computing the r* statistic. A positive integer, default is \code{1000}.
    }
\item{seed}{
    Optional positive integer, the random seed for the Monte Carlo computation. Default is \code{NULL}.
    }  
\item{trace}{
     Logical. When set to \code{TRUE} will cause some information on the computation to be printed. Default is \code{FALSE}.
    }   
\item{ronly}{
     Logical. If set to \code{TRUE} the computation of the r* statistic will be skipped, and only the value of the 
     signed  likelihood ratio test statistic r will be returned by the procedure, without any Monte Carlo computation. 
     Default is \code{FALSE}.
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
}
 \value{
  The returned value is an object of class \code{"rstar"}, containing the following components:

  \item{r}{
    The observed value the signed likelihodo ratio test statistic r for testing \code{fpsi(theta)=psival}. 
  }
  \item{NP, INF}{
    The Nuisance Parameter adjustment (NP) and the Information adjustment (INF) from the decomposition of the r*-r adjustment. The former is not computed for one-parameter models. Neither one is computed when \code{ronly = TRUE}.
  }
  \item{rs}{
    The observed value of the r* statistic.  Not computed when \code{ronly = TRUE}.
  }
  \item{theta.hat}{
    The maximum likelihood estimate of the parameter theta, the argument of the \code{floglik}, \code{fscore}, \code{datagen} and 
    \code{fpsi} functions.
  }
  \item{info.hat}{
    The observed information matrix evaluated at \code{theta.hat}. Not computed when \code{ronly = TRUE}.
  } 
  
  \item{se.theta.hat}{
    The estimated standard error of \code{theta.hat}. Not computed when \code{ronly = TRUE}.  
  }
  \item{psi.hat}{
    The parameter of interest evaluated at \code{theta.hat}. 
  }
  \item{se.psi.hat}{
    The estimated standard error for the parameter of interest. Not computed when \code{ronly = TRUE}.
  }
  \item{theta.hyp}{
    The constrained estimate of the parameter, under the null hypothesis \code{fpsi(theta)=psival}. 
  }
  \item{psi.hyp}{
    The value under testing, equal to \code{psival}.
  }
  \item{seed}{
    Random seed used for Monte Carlo trials. Not returned when \code{ronly = TRUE}.
  }
  \item{psidesc}{
   A character string describing the nature of the parameter of interest. 
  }
  \item{R}{
    Number of Monte Carlo replicates used for computing the r* statistic. Not returned when \code{ronly = TRUE}.
  }   
 
  
  There are \code{print} and \code{summary} methods for this class. 
}
\details{
  The function computes the r* statistic proposed by Skovgaard (1996) for accurate 
  computation of the asymptotic distribution of the signed likelihood ratio test for
  a scalar function of interest. 
    
  The function requires the user to provide three functions defining the log likelihood function,
  the scalar parametric function of interest, and a function for generating
  a data set from the assumed statistical model. A further function returning the gradient of the 
  log likelihood is not required, but if provided it will speed up the computation. 
  
  When \code{ronly = TRUE} the function returns the value of the signed likelihood ratio test
  statistic r only. 
  
  The function handles also one-parameter models.
}

\references{
The method implemented in this function was proposed in

Skovgaard, I.M. (1996) An explicit large-deviation approximation to one-parameter tests. 
\emph{Bernoulli}, \bold{2}, 145--165.

For a general review

Severini, T.A. (2000). Likelihood Methods in Statistics. Oxford University Press.
}
\seealso{
  \code{\link{rstar.ci}}.
}
\examples{
# Autoregressive model of order 1
# We use the lh data from MASS
library(MASS)
data(lh)
dat.y <- list(y=as.numeric(lh))
# First let us define the function returning the log likelihood function
# We employ careful parameterizations for the correlation and variance to
# avoid numerical problems
likAR1 <- function(theta, data)
{ 
  y <- data$y
  mu <- theta[1]
  phi <- theta[2] ### phi is log(sigma) 
  sigma2 <- exp(phi*2)
  z <- theta[3]   ### z is Fisher'z transform for rho
  rho <- (exp(2*z)-1) / (1 + exp(2*z))
  n <- length(y)
  Gamma1 <- diag(1+c(0,rep(rho^2,n-2),0))
  for(i in 2:n)
    Gamma1[i,i-1]<- Gamma1[i-1,i] <- -rho 
  lik <- -n/2 * log(sigma2) + 0.5 * log(1-rho^2) -1/(2*sigma2) * 
        mahalanobis(y, rep(mu,n), Gamma1, inverted = TRUE)
  return(lik)
}
# We need a function for simulating a data set
genDataAR1 <- function(theta, data)  
{
  out <- data
  mu <- theta[1]
  sigma <- exp(theta[2])
  z <- theta[3]
  rho <- (exp(2*z)-1) / (1 + exp(2*z))
  n <- length(data$y)
  y <- rep(0,n)
  y[1] <- rnorm(1,mu,s=sigma*sqrt(1/(1-rho^2)))
  for(i in 2:n)
    y[i] <- mu + rho * (y[i-1]-mu) + rnorm(1) * sigma 
  out$y <- y 
  return(out)
}
# For inference on the mean parameter we need a function returning the first component of theta
psifcn.mu <- function(theta) theta[1]
# Now we can call the function
rs.mu <- rstar(dat.y, c(0,0,0), likAR1, fpsi=psifcn.mu, psival=2, datagen=genDataAR1, R=1000, 
               trace=TRUE, psidesc="mean parameter")
summary(rs.mu)
}
\keyword{htest}