# dil.r is dilution ratio, previously named k

logcdbf=c("logc","logd","b","f")

# verbose=1: only one line per iteration
# verbose=2: several short lines per iteration
# verbose=4: print optim output
# reltol=1e-3; max.iter=50; verbose=3; model="4P"; grid.density=200
prcsp=function(xvar, dil.x, yvar, dil.y, grid.density=200, init.method=c("gnls","optim"), 
    reltol=1e-3, max.iter=20, model=c("4P"), init=NULL, try.additiona.support.sets=TRUE, verbose=FALSE) {    
    
    # removing this line requires adding package namespace to Matrix and mosek function calls
    have.Rmosek=try(require(Rmosek))
    if (!have.Rmosek) stop("have.Rmosek does not load successfully. The package is not on CRAN and needs to be installed.")
    
    
    opt.method="nlminb" # optim or nlnimb does not seem to make a big difference
    init.method=match.arg(init.method)
    if (verbose) {myprint(opt.method, init.method)}
    
    stop.when.dropping=FALSE
    
    dil.r=dil.x/dil.y
    stopifnot(length(xvar)==length(yvar))
    n=length(xvar)    
    xvar=as.vector(xvar); yvar=as.vector(yvar) # strip off potential attributes; otherwise it may cause trouble with gnls
    dat=data.frame(readout.x=unname(xvar), readout.y=unname(yvar))    
    
    res=list(dilution.ratio=dil.r, dilution.x=dil.x, dilution.y=dil.y, xvar=xvar, yvar=yvar)
    class(res)=c("prc",class(res))
    
    #### run prc for a small number of iterations to initialize
    if (is.null(init)) {
        if (verbose) cat("\nrun prc to get initial theta\n")
        fit.init=suppressWarnings( prc (xvar, dil.x, yvar, dil.y, max.iter=2, init.method=init.method, verbose=verbose) )
        theta=c(coef(fit.init), sigma.sq=fit.init$sigma.sq)
        theta=c(logc=unname(log(theta["c"])), logd=unname(log(theta["d"])), theta[-match(c("c","d"), names(theta))])
    } else {
        theta=init
        if (length(theta)!=5) stop("init need to include sigma.sq")
    }
        
    if (verbose) cat("\nInitial theta:", theta, "\n")
    iterations=0
    old.lik=-Inf
    old.support=NULL
    old.p=NULL
    histories=list()
    mix.lik.after.thetas=NULL
    
    while(TRUE) {
    
        iterations=iterations+1
        if (verbose>=2) {
            cat("=========== Iter "%+%iterations%+%" ==============\n") 
        } else {
            cat(formatC(iterations,digits=2,), ") ", sep="")
        }
        
        # Step 1: get nonparametric estimate of distribution of u
        
        # define the grid and compute A
        u=seq((theta["logc"]), (theta["logd"]), length=2+grid.density)
        u=u[-c(1,grid.density)] # length of u is K
        A=compute.A (theta["logc"], theta["logd"], theta["b"], theta["f"], dil.r, theta["sigma.sq"], u, xvar, yvar)         
        
        # run mosek
        sco1 <- list(sense = "min")
        sco1$c <- rep(0,n)
        sco1$A <- Matrix::Matrix (t(A), sparse = TRUE )
        sco1$bc <- rbind(blc = rep(-Inf,length(u)), buc = rep(n,length(u)))
        sco1$bx <- rbind(blx = rep(0,n), bux = rep(Inf,n))
        opro <- matrix(list(), nrow=5, ncol=n, dimnames=list(c("type","j","f","g","h"), NULL))
        for (i in 1:n) opro[,i] <- list("LOG", i, -1.0, 1.0, 0.0)
        sco1$scopt <- list(opro=opro)
        fit.mosek <- Rmosek::mosek(sco1, opts=list(verbose=1))
        if (fit.mosek$sol$itr$solsta!="OPTIMAL") warning("fit.mosek$sol$itr$solsta!=OPTIMAL")
        if (fit.mosek$sol$itr$prosta!="PRIMAL_AND_DUAL_FEASIBLE") warning("fit.mosek$sol$itr$prosta!=PRIMAL_AND_DUAL_FEASIBLE")
        #mosek_write (sco1, "mosekfiles/yf_ex1_run_1.task", list(scofile="yf_ex1.sco"))
        mix.lik.dual = sum(-log(fit.mosek$sol$itr$xx))
        
        if (verbose>=2) {
            myprint(mix.lik.dual)               
            plot(fit.mosek$sol$itr$xc, sco1$A %*% fit.mosek$sol$itr$xx); abline(0,1) # sanity check
            title(main="Iteration "%+%iterations, outer=T, line=-1)
            empty.plot()
            #print(fit.mosek)
            #print(table(fit.mosek$sol$itr$skc))        
#            print(table(fit.mosek$sol$itr$skx))
#            print(summary(fit.mosek$sol$itr$xx))
#            print(summary(fit.mosek$sol$itr$xc))
#            print(cbind(fit.mosek$sol$itr$skx, fit.mosek$sol$itr$xx, fit.mosek$sol$itr$slx, fit.mosek$sol$itr$sux))
            if(verbose>=4) print(cbind(fit.mosek$sol$itr$skc, fit.mosek$sol$itr$xc, fit.mosek$sol$itr$slc, fit.mosek$sol$itr$suc)[fit.mosek$sol$itr$skc=="UL",])
        }
        
        #print(sort(fit.mosek$sol$itr$xc, decreasing=TRUE))
    
        # recover primal soln
        support.set.1 = suppressWarnings(which(fit.mosek$sol$itr$skc=="UL")) # UL seems to choose less than ideal set
        if (length(support.set.1)<=1) {
            support.set.1 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:floor(n*.1)] # 10% sample size
            support.set.2 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:floor(n*.2)] # 20% sample size
            support.set.3 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:floor(n*.4)] # 40% sample size
        } else {
            support.set.2 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:min(length(support.set.1)*1.5,n*.8)] # twice as many as support.set.1, but not over 80% sample size
            support.set.3 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:min(length(support.set.1)*2.0,n*.8)] # four times as many as support.set.1, but not over 80% sample size
        }
        support.set.5 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:which(cumsum(sort(fit.mosek$sol$itr$suc,decreasing=TRUE))>=1-1e-6)[1]] # cumsum of suc just above 1
        if (verbose>=2) cat("support sizes UL|1|2|3|5: ", suppressWarnings(sum(fit.mosek$sol$itr$skc=="UL")), length(support.set.1), length(support.set.2), length(support.set.3), length(support.set.5), "\n")
        # Others tried:
        # fit.mosek$sol$itr$xc==n sometimes return empty set
        # a modification of set.1 that includes all with suc greater than any UL suc, does not make a big difference
#        support.set.2 = which(fit.mosek$sol$itr$suc>0.01) 
#        support.set.3 = which(fit.mosek$sol$itr$suc>0.001) 
        
        # try three choices of support set and use the one with the best mix.lik.primal
        support.set.list=list()
        if (length(support.set.2)<length(support.set.1)) num.sets=1 else if (length(support.set.3)<length(support.set.2)) num.sets=2 else num.sets=3
        if (!try.additiona.support.sets) num.sets=1
        mix.lik.primals=numeric(num.sets)
        for (support.idx in 1:num.sets) {
            if (verbose>=2) myprint(support.idx)
            support.set = get("support.set."%+%support.idx)            
            while (TRUE) {
                p.fit=lm.fit(x = A[,support.set,drop=FALSE], y = 1/fit.mosek$sol$itr$xx)
                p.new=coef(p.fit)
                if (any(is.na(p.new))) { 
                    #stop("some p.new are NA")
                    p.new[is.na(p.new)]=-1
                } 
#               p.new=fit.mosek$sol$itr$suc[support.set] # performance not good
#                support.1=unique(r.x); p.1=rep(1/5,5) # mix.lik not as good as support
#                support.1=support[-c(2,4)]; p.1=rep(1/5,5) # mix.lik not as good as support
#                sum(log(compute.A (theta["logc"],theta["logd"],theta["b"],theta["f"], dil.r, theta["sigma.sq"], support.1, xvar, yvar) %*%p.1))
                if (verbose>=2) cat("sum(p.new)=", round(sum(p.new),2), ", ", sep="")
                p.new=p.new/sum(p.new)
                mix.lik.primal = suppressWarnings(sum(log(compute.A (theta["logc"],theta["logd"],theta["b"],theta["f"], dil.r, theta["sigma.sq"], u[support.set], xvar, yvar) %*%p.new)))                
                if (verbose>=2) {
                    myprint(mix.lik.primal, digits=6)
                    plot(p.fit$fitted.values, 1/fit.mosek$sol$itr$xx) # making sure we are doing a good job at recovering primal soln
                    plot(0,0, type="n", ylim=range(c(0,p.new)), xlim=range(u[support.set]), xlab="r", ylab="p")#, main="Iteration "%+%iterations)
                    points(u[support.set], p.new, pch=19)
                    for (i in 1:length(p.new)) lines(rep(u[support.set][i],2), c(0,p.new[i]))            
                    if (verbose>=3) {
                        print(cbind(fit.mosek$sol$itr$skc, fit.mosek$sol$itr$xc, fit.mosek$sol$itr$slc, fit.mosek$sol$itr$suc)[order(fit.mosek$sol$itr$suc,decreasing=TRUE)[1:n],])
                        myprint(sum(fit.mosek$sol$itr$suc[support.set]))
                        myprint(sum(fit.mosek$sol$itr$suc[support.set.5]))
                        print(rbind(u[support.set], p.new))
                    }                            
                }
                
                support.set=support.set[p.new>0]
                if (all(p.new>0)) break 
            }
            support.set.list[[support.idx]]=support.set # save this
            old.par=par(no.readonly = TRUE); par(old.par) # to start a new page with the existing graphics settings, but the order of mfrow vs mfcol is not saved
            mix.lik.primals[support.idx]=mix.lik.primal            
        } # end support.idx loop
        mix.lik.primal=max(mix.lik.primals)
        if (verbose>=2) myprint(mix.lik.primal, digits=6)
        if (mix.lik.primal<mix.lik.dual/2) {
            cat("\nStop: unable to extract a good primal soln\n")
            break;
        } else {
            support.set = support.set.list[[which.max(mix.lik.primals)]]
            p.fit=lm.fit(x = A[,support.set,drop=FALSE], y = 1/fit.mosek$sol$itr$xx)
            p=coef(p.fit)
            p=p/sum(p)
            support=u[support.set]
        }
        
        ####################################################
        # loop through step 2 and 3
        
        # initialize new.theta with an approximation to the mixture likelihood 
        A=compute.A (theta["logc"], theta["logd"], theta["b"], theta["f"], dil.r, theta["sigma.sq"], support, xvar, yvar) 
        best.support=support[apply(A, 1, which.max)]
        s.sqrt.k.rvar = four_pl_prc (exp(theta["logc"]), exp(theta["logd"]), theta["b"], theta["f"], best.support, dil.r^.5) # an easy to make mistake is to use k*.5 instead of k^.5
        dat.stacked=data.frame(readout=c(unname(xvar),unname(yvar)), x=rep(s.sqrt.k.rvar,2), k=rep(c(dil.r^(-1/2),dil.r^(1/2)),each=n))
        formula.gnls = as.formula(  "(readout) ~ log(exp(logc)+(exp(logd)-exp(logc))/(1+k^b*(((exp(logd)-exp(logc))/(exp(x)-exp(logc)))^(1/f)-1))^f)"  ) 
        fit.1=try(gnls(formula.gnls, data=dat.stacked, start=theta[logcdbf], # key to take only cdbf from new.theta as start for gnls
            control=gnlsControl(
            nlsTol=1e-1,  # nlsTol seems to be important, if set to 0.01, then often does not converge
            tolerance=1e-4, 
            # msTol=1e-1, minScale=1e-1, .relStep=1e-7,
            returnObject=TRUE, # allow the return of the fit when the max iter is reached
            maxIter=5000, nlsMaxIter=50, opt="nlminb", msVerbose=T)), silent=FALSE
        )       
        new.theta=coef(fit.1)
        new.theta=c(new.theta, theta["sigma.sq"]) # has to be on two lines, otherwise names are lost
        if (verbose>=2) cat("theta after gnls:", new.theta, "\n")        
        if (any(support<new.theta["logc"]) | any(support>new.theta["logd"])) {
            if (verbose>=2) cat("some support is outside new.theta, revert new.theta to theta\n")
            new.theta=theta
        } 
        
        max.inner.iter=25
        inner.iter=0
        while (TRUE) {
            if (inner.iter>max.inner.iter) break;
            inner.iter=inner.iter+1
            new.theta.0=new.theta # used to decide convergence
            
            # Step 2: estimate sigma.sq
            #m=compute.m (theta["logc"], theta["logd"], theta["b"], theta["f"], dil.r, sigma.sq, support, xvar, yvar) 
            optimize.out = optimize(
                f = function(sigma.sq.f,...) {
                    A=compute.A (new.theta["logc"], new.theta["logd"], new.theta["b"], new.theta["f"], dil.r, sigma.sq.f, support, xvar, yvar) 
                    -suppressWarnings(mean(log(A%*%p)))
                },
                interval=c(new.theta["sigma.sq"]/100, new.theta["sigma.sq"]*100), 
            )
            if (verbose>=4) print(optimize.out)
            new.theta["sigma.sq"] = optimize.out$minimum            
            
            # Step 3: estimate curve parameters
            if (opt.method=="nlminb") {
                optim.out = suppressWarnings(nlminb(
                    new.theta[logcdbf],
                    function(theta.fn,...) {
                        c=
                        A=compute.A (theta.fn["logc"], theta.fn["logd"], theta.fn["b"], theta.fn["f"], dil.r, new.theta["sigma.sq"], support, xvar, yvar) 
                        -(mean(log(A%*%p)))
                    }, 
                    gradient=NULL, hessian = F,
                    #p, support, xvar, yvar, dil.r, 
                    control = list(trace=0)
                )) # it will warn about NA/NaN function evaluation
                if (verbose>=4) print(optim.out)
                new.theta[logcdbf] = optim.out$par
        
            } else if (opt.method=="optim") {
                optim.out = optim(
                    new.theta[logcdbf],
                    function(theta.fn,...) {
                        A=compute.A (theta.fn["logc"],theta.fn["logd"],theta.fn["b"],theta.fn["f"], dil.r, new.theta["sigma.sq"], support, xvar, yvar) 
                        -suppressWarnings(mean(log(A%*%p)))
                    }, 
                    gr=NULL,
                    #p, support, xvar, yvar, dil.r, 
                    method="BFGS", control = list(trace=0), hessian = F
                )
                if (verbose>=4) print(optim.out)
                new.theta[logcdbf] = optim.out$par
            }
            
            if (verbose>=3) myprint(new.theta)
            if (max(abs(1 - new.theta/new.theta.0)) < reltol) {
                break;
            }            
        } # end while step 2 and 3
        mix.lik.after.theta = suppressWarnings(sum(log(compute.A (new.theta["logc"],new.theta["logd"],new.theta["b"],new.theta["f"], dil.r, new.theta["sigma.sq"], support, xvar, yvar) %*%p)))
            
        if (verbose>=2) {
            myprint(mix.lik.after.theta, digits=6)
            cat("support size:", length(support), "\n")
            cat("theta:", new.theta)#, "support:", mean(support))
            cat("\n")            
        } else if (verbose) {
            cat("ll", formatC(mix.lik.after.theta,format="g",width=7))
            cat(", #support", formatC(length(support),format="g",digits=ceiling(log10(n))), sep="")
            cat(", theta", new.theta, "\n")                        
        }
        
        histories[[iterations]]=list(support=support, p=p, theta=new.theta, mix.lik.after.theta=mix.lik.after.theta)
        mix.lik.after.thetas=c(mix.lik.after.thetas, mix.lik.after.theta)
        if (max(abs(1 - new.theta/theta)) < reltol) {
            if (verbose) cat("converged\n")
            theta = new.theta # update theta 
            break;
        } else if (iterations>=max.iter) {
            cat("Stop: max iter reached\n")
            theta = new.theta # update theta 
            break;
        } else if (stop.when.dropping & mix.lik.after.theta<old.lik) {
            cat("Stop: likelihood starts decreasing.\n")
            # do not update theta 
            # revert to old values
            mix.lik.after.theta=old.lik
            support=old.support
            p=old.p
            break;
        } else {
            theta = new.theta # update theta 
            # update old values
            old.lik=mix.lik.after.theta 
            old.support=support
            old.p=p
        }        
            
    } # end while loop
        
    res$iterations=iterations
    res$histories=histories
    res$mix.lik.after.thetas=mix.lik.after.thetas
    
    best.iter=histories[[which.max(mix.lik.after.thetas)]]    
    res$support=best.iter$support
    res$p=best.iter$p
    res$coefficients=best.iter$theta[logcdbf]
    res$sigma.sq = best.iter$theta["sigma.sq"]
    res$mixlik = best.iter$mix.lik.after.theta
    
    res$A=with(res, compute.A (coefficients["logc"],coefficients["logd"],coefficients["b"],coefficients["f"], dil.r, sigma.sq, support, xvar, yvar) )
    
    #### compute asymptotic variance of theta_hat and xhat and yhat
    
    
        
    # return object
    res
}

# this should match the mixture likehood obtained directly from mosek
mixlik <- function(object, ...) UseMethod("mixlik") 
mixlik.prc=function(object, ...) {    
    with(object, sum(log(compute.A (coefficients["logc"],coefficients["logd"],coefficients["b"],coefficients["f"], dil.r, sigma.sq, support, xvar, yvar) %*%p)))
    
}

# populate A
compute.A = function(logc,logd,b,f, dil.r, sigma.sq, support, xvar, yvar) {
    c=exp(logc); d=exp(logd)
    .Call("compute_A", c,d,b,f, dil.r, sigma.sq, support, xvar, yvar)
    # .Call is a lot faster than the R implementation below
#    m.f.1=function(c,d,b,f,r,x,y,dil.r)  (y - four_pl_prc(c,d,b,f, r, dil.r))^2  +  (x - r)^2 # this is about four times as fast as m.f
#    A=matrix(NA,n,K)
#    for (i in 1:n) {
#        for (k in 1:K) {
#            m_ik = m.f.1(c, d, b, f, u[k], xvar[i], yvar[i], dil.r) 
#            A[i,k] = 1/fit$sigma.sq * exp(-m_ik/fit$sigma.sq/2)
#        }
#    }
}


compute.m = function(logc,logd,b,f, dil.r, sigma.sq, support, xvar, yvar) {
    c=exp(logc); d=exp(logd)
    n=length(xvar)
    m.f.1=function(c,d,b,f,r,x,y,dil.r)  (y - four_pl_prc(c,d,b,f, r, dil.r))^2  +  (x - r)^2 # this is about four times as fast as m.f
    K=length(support)
    m=matrix(NA,n,K)
    for (i in 1:n) {
        for (k in 1:K) {
            m[i,k] = m.f.1(c, d, b, f, support[k], xvar[i], yvar[i], dil.r) 
        }
    }
    m
}
