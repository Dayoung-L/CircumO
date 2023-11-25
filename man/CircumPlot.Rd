\name{CircumPlot}
\alias{CircumPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Circumplex Model Plot
}
\description{
The function is to graphically represent ordinary least square estimates of circumplex model parameters.
}
\usage{
CircumPlot <-
  function(C.object,p,m=1,r=1pchar=NULL,bg.points="red",ef=0.4,big.points=15,big.labels=15,bg.plot="white",col.axis="black",color="black",col.text="white",labels=TRUE,reverse=FALSE,title=NULL)}
\arguments{
  \item{C.object}{indicates the result object of "CircumO".}

  \item{m}{indicates the number of cosine function coefficients. Defaults to 1.}

  \item{r}{indicates the order of the reference variable. Defaults to 1.}

  \item{labels}{ logical. If \code{TRUE} (default), variable labels are reported in the plot.}

  \item{pchar}{ plotting character, i.e., symbol to use. It can be either   a single value or a vector of length equal to numbers of variables used in correlation matrix. For further details see \code{?par}.  }

  \item{bg.points}{ background color for the symbols selected with \code{pchar}. Either a single value or a vector of color.}

  \item{ef}{ "explosion factor". It specifies the distances at which the name of the variables are drawn around the circumference. Default to 0.4 (if negative labels are drawn inside the circumference).}

  \item{big.points, big.labels}{ numerical values giving the amount by which plotting labels and symbols should be magnified relative to the default (\code{cex} = 1/(number of observed variables)). }

  \item{bg.plot}{ background color of the plot region. }

  \item{col.axis}{ the color to be used for axis annotation, axis line and labels. }

  \item{color}{ the color to be used for plot lines and symbols foreground. }

  \item{col.text}{ the color to be used for plot texts. }

  \item{label}{ indicates a vector object of the variable labels. If \code{NULL} (default), variable labels are reported according to p. }

  \item{reverse}{ logical. If \code{TRUE}, plot the figure using the 360-angular positions (just reversing the circular order).}
}

\details{
CirucmPlot

Code modified from function CircE.Plot obtained from https://cran.r-project.org/src/contrib/Archive/CircE/CircE_1.1.tar.gz
}

\author{
Dayoung Lee \email{dlee33@nd.edu}
}

\examples{
#An example using the data set included in the package:

data("emotion") # data from (Trampe et al., 2015), N = 12108

result = CircumO(emotion, m = 1, type = "O")

CircumPlot(result, p = 18)
}
