\name{prc}
\alias{prc}
\alias{coef.prc}
\alias{predict.prc}
\alias{plot.prc}
\alias{print.prc}
\alias{lines.prc}
\alias{four_pl_prc}
\alias{s.dot.f}
\alias{prcsp}
\alias{compute.A}
\alias{mixlik}
\alias{mixlik.prc}
\alias{quadratic.eiv}
\alias{coef.quad}
\alias{lines.quad}
\alias{plot.quad}
\alias{print.quad}
\alias{quad.f}
\alias{quadratic.eiv.sp}
\title{
Four Parameter Paired Response Curve
}
\description{
Estimate four-parameter paired response curve using both total least squares and semiparametric methods.
}
\usage{

prc(xvar, dil.x, yvar, dil.y, model=c("4P","3P"), method=c("M-estimator","naive"),
  init.method=c("gnls","optim"), opt.method=c("gnls","optim"), reltol=1e-3, max.iter=50,
  init=NULL, verbose=FALSE)


prcsp (xvar, dil.x, yvar, dil.y, grid.density=200, init.method=c("gnls","optim"),
    reltol=1e-3, max.iter=20, model=c("4P"), init=NULL, try.additiona.support.sets=TRUE,
    verbose=FALSE)

mixlik (object, ...) 
\method{mixlik}{prc}(object, ...) 

compute.A (logc,logd,b,f, dil.r, sigma.sq, support, xvar, yvar) 

    
\method{predict}{prc}(object, new.dilution, xvar=NULL, dil.x=NULL, ret.sd=FALSE, ...)
\method{coef}{prc}(object, ...)
\method{print}{prc}(x, ...)
\method{plot}{prc}(x, type=c("b","l","p"), add=FALSE, diag.line=TRUE, lcol=2, pcol=1,
    log.axis=TRUE, xlab=NULL, ylab=NULL, lwd=2, xlim=NULL, ylim=NULL, ...) 
\method{lines}{prc}(x, col=1, ...) 

four_pl_prc(c, d, b, f, xx, k)
s.dot.f(c,d,b,f, r, k)

quadratic.eiv (xvar, yvar, init=NULL, reltol=1e-3, opt.method=c("optim"), 
    stop.after.init=FALSE, max.iter=50, verbose=FALSE)
\method{coef}{quad}(object, ...)
\method{print}{quad}(x, ...)
\method{plot}{quad}(x, type=c("b","l"), add=FALSE, lcol=2, pcol=1, 
    xlab=NULL, ylab=NULL, lwd=2, x.range=NULL, log.axis=TRUE,...) 
\method{lines}{quad}(x, col=1, x.range=NULL, ...) 

quadratic.eiv.sp (xvar, yvar, grid.density=200, init=NULL, reltol=1e-3, 
    opt.method=c("optim"), max.iter=50, fix.sigma.sq=FALSE, verbose=FALSE)

quad.f(a,b,c,x) 



}
\arguments{
  \item{xvar}{
vector of numbers on the ln(FI) scale
}
  \item{dil.x}{
numeric.
}
  \item{yvar}{
vector of numbers on the ln(FI) scale
}
  \item{dil.y}{
numeric.
}
  \item{init.method}{gnls uses gnls to find initial parameter estimate and optim uses optim. In the initial estimation, measurement errors in the xvar are ignored.}
  \item{reltol}{relative tolerance convergence criterion used in optimization}
  \item{opt.method}{optimization method}
  \item{method}{M-estimator minimizes total least squares. Naive minmizes least squares, ignoring measurement errors in the xvar.}
  \item{model}{4P is a four-parameter model proposed in Fong et al. 3P is a three-parameter model that assumes f=1 in 4P.}
  \item{max.iter}{maximum number of iterations to perform in optimization. Each iteration comprises finding the best r given theta and finding the best theta given r}
  \item{verbose}{Boolean. If TRUE, print messages during execution}

  \item{object}{an object of type prc}
  \item{new.dilution}{a new dilution}
  \item{type}{"b" for both points and line, "l" for line only}
  \item{add}{add to an existing plot or not}
  \item{lcol}{line color}
  \item{pcol}{points color}
  \item{col}{color}
  \item{diag.line}{whether to add a diagonal line}
  \item{log.axis}{whether the axis should be on the log scale}
  \item{xlab}{xlab}
  \item{ylab}{ylab}
  \item{lwd}{line width. Does not affect points}

  \item{x}{object of class prc or quad depending on the functions}

  \item{xx}{vector of numbers on the log scale}
  \item{k}{dilution ratio: dilution of x axis/dilution of y axis, e.g. 2500/500}
  \item{c}{parameter c in the four-parameter paired response model}
  \item{d}{parameter d in the four-parameter paired response model}
  \item{b}{parameter b in the four-parameter paired response model}
  \item{f}{parameter f in the four-parameter paired response model}

  \item{a}{parameter a in the quadratic model}
  \item{r}{the x axis of a point on the prc on the log scale}
  
  \item{init}{initial parameter value
  }
  \item{ret.sd}{whether to return standard deviation
  }
  \item{xlim}{xlim
  }
  \item{ylim}{ylim
  }

  \item{logc}{log(c) in the four-parameter paired response model}
  \item{logd}{log(d) in the four-parameter paired response model}
  \item{dil.r}{dilution ratio}
  \item{sigma.sq}{sigma square}
  \item{support}{the support of the distribution}


  \item{try.additiona.support.sets}{boolean}
  \item{stop.after.init}{Boolean. If true (not recommended), a naive curve estimate which pretends there is no measurement error in xvar is returned. If false, a proper least squares estimate accounting for errors-in-variables is returned.}
  \item{x.range}{a pair of low and high}

  \item{grid.density}{number of grid points}
  \item{fix.sigma.sq}{whether to assume simgasq is fixed}

  \item{...}{more args}
}
\details{
prc() operates on the scale on which the distance is defined. For example, for MBA/Luminex readouts, this is the log(FI) scale.
s.dot.f() is created with deriv3 and it returns prc function value, as well as gradient and hessian at the point (r,s(r))

quadratic.eiv and quadratic.eiv.sp deal with a quadratic model as an illustration
}
\value{
Both prc and prcsp return an object of type prc.
  \item{coefficients}{curve estimate}
  \item{sigma.sq}{error variance estimate}
  \item{dilution.ratio}{dilution.x/dilution.y}
  \item{dilution.x}{sample dilution for the readouts plotted on the x axis}
  \item{dilution.y}{sample dilution for the readouts plotted on the y axis}
  \item{xvar}{observations on the x axis}
  \item{yvar}{observations on the y axis}
  \item{rvar}{closest point to (x,y) on the curve is (r,s(r))}
  
Both quadratic.eiv and quadratic.eiv.sp return an object of type quad, which has a similar structure as type prc.

compute.A returns a matrix.
}
\author{
Youyi Fong
}
\references{
Fong, Y., Permar, S., Tomaras, G. (2014) Four-Parameter Paired Response Curve for Serial Dilution Assay: Estimation and Prediction. Submitted
}

\examples{

# estimation
fit=prc (mtct.eg$V3_BioV3B_2500, 2500, mtct.eg$V3_BioV3B_500, 500, verbose=TRUE) 
plot(fit)

# prediction
logfi.1000 = predict(fit, new.dilution=1000)
points(exp(fit$rvar), exp(logfi.1000), col=2, cex=.5)

theta=coef(fit)
four_pl_prc(theta["c"], theta["d"], theta["b"], theta["f"], log(500), k=5)

}
