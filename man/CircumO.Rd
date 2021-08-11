\name{CircumO}
\alias{CircumO}
\title{Circumplex Models Ordinary Least Square Estimation}
\usage{
CircumO(rawdt, m=1, mcsc="unconstrained", type="ordinal", simulation_based = F, icorrection = 0, ncore = 8, maxit = 1000, factr = 1e7 ,pgtol = 0, lmm = NULL, N_star = 10000)
}
\arguments{
\item{rawdt}{indicates the raw data which is a n by p matrix where n is the number of participants and p is the number of manifest variables.}

\item{m}{indicates the number of cosine function coefficients. Defaults to 1.}

\item{mcsc}{minimum common score correlation value: \code{"unconstrained"} (default), \code{-1}}

\item{type}{indicates the type of data: \code{"ordinal"} (default), \code{"continuous"}}

\item{simulation_based}{determines whether an asymptotic covariance matrix estimate is a sample-based ACM estimate or a Monte Carlo ACM estimate: \code{False} (default, sample-based ACM), \code{True} (Monte Carlo ACM)}

\item{icorrection}{Methods to adjust for empty cells: a scalar where 0 is no adjustment is done (default), 1  adds 1/(nc*nr) to all cells where nc and nr is the number of columns and rows of the contingency table respectively, 2 adds 0.1 to all cells, 3 adds 0.5 to all cells, 11 adds 1/(nc*nr) to only zero cells, 12 adds 0.1 to only zero cells, and 13 adds 0.5 to only zero cells}

\item{ncore}{indicates the number of cores for parallel computing. It defaults to 8.}

\item{maxit}{The maximum number of iterations. It defaults to 1000.}

\item{factr}{controls the convergence of the "L-BFGS-B" method. Convergence occurs when the reduction in the objective is within this factor of the machine tolerance. Default is 1e7, that is a tolerance of about 1e-8.}

\item{pgtol}{helps control the convergence of the "L-BFGS-B" method. It is a tolerance on the projected gradient in the current search direction. This defaults to zero, when the check is suppressed.}

\item{lmm}{is an integer giving the number of BFGS updates retained in the "L-BFGS-B" method. It defaults to number of free parameters.}

\item{N_star}{is the sample size of a simulated sample for Monte Carlo ACM. It defaults to 10000.}
}
\description{
The function is to compute ordinary least square estimates of circumplex model parameters and their standard errors. It allows both of continuous and ordinal variables.
}
\details{
CirucmOLS

Code modified from function CircE.BFGS obtained from https://cran.r-project.org/src/contrib/Archive/CircE/CircE_1.1.tar.gz
}
\references{
Browne, M. W. (1992). Circumplex models for correlation matrices. Psychometrika, 57, 469-497. doi: 10.1007/BF02294416

Monroe, S. (2018). Contributions to estimation of polychoric correlations. Multivariate Behavioral Research, 53, 247-266. doi: 10.1080/00273171.2017.1419851

Grassi, M., Luccio, R., & Di Blas, L. (2010). Circe: An r implementation of browne’s circular stochastic process model. Behavior Research Methods, 42(1), 55-73. doi:10.3758/BRM.42.1.55

Lee, D., & Zhang, G. (in preparation). Circumplex models with ordinal data.

Zhang, G., Trichtinger, L., Lee, D., & Jiang, G. (2021). Polychoricrm: A computationally efficient r function for estimating polychoric correlations and their asymptotic covariance matrix. Structural Equation Modeling: A Multidisciplinary Journal. doi:10.1080/10705511.2021.1929996
}
\author{
Dayoung Lee \email{dlee33@nd.edu}
}

\examples{
#An example using the data set included in the package:

data("emotion") # data from (Trampe et al., 2015), N = 12108

CircumOLS(emotion, m = 1, type = "ordinal")
}
