\name{rstar.ci.control}
\alias{rstar.ci.control}
\title{
 Auxiliary function for controlling computation of r*-based confidence intervals	
}
\description{
  Auxiliary function for \code{rstar.ci}. 
}
\usage{
rstar.ci.control(npoints=10, away=0.3, stepsizefac=3, maxstep=50, trace=TRUE)
}
\arguments{
  \item{npoints}{
    Integer giving the number of points at which the r* and r statistics will be evaluated away from the MLE 
    in each direction when both \code{lower} and \code{upper} are not null. When either \code{lower} or \code{upper} 
    are \code{NULL}, such value is only roughly proportional to the number of evaluation points. Default is \code{10}.
    }
  \item{away}{
    Positive value indicating how far from the MLE the computation of the r* and r statistics will be started, expressed
    in units of standard error of the scalar function of interest. Default is \code{0.3}.
    }
  \item{stepsizefac}{
    Positive value used to determine the stepsize of the confidence interval algorithm when either \code{lower} or 
    \code{upper} are null. In particular, the stepsize is given by \code{stepsizefac/npoints} times the 
    standard error of the scalar function of interest. Default is \code{3}.
    }
 \item{maxstep}{
    Integer giving the maximum number of steps that will be taken for crossing the normal quantiles for a 99%
    confidence interval for the r* statistic. Default is \code{50}.
    }
\item{trace}{
     Logical indicating if output should be produced during the computation. Default is \code{TRUE}.
    }   
}
 \value{
  A list with components named as the arguments. 
  }
\seealso{
  \code{\link{rstar.ci}}.
}
\examples{
# A variation on example(rstar.ci):
\dontrun{
obj <- rstar.ci(quinedata, thetainit=c(coef(quine.nb1),quine.nb1$theta),  floglik=logLikNbin, 
                datagen=genDataNbin, fscore=gradLikNbin, fpsi=function(theta) theta[2], R=1000, 
                psidesc="Coefficient of EthN", npoints=5, away=0.1)
plot(obj)
}
}
\keyword{htest}