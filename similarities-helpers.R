# constrOptim funciton for stats package, modified to avoid feasible region problem
# caused by numerical inaccuracy of calculations
myConstrOptim = function(theta, f, grad, ui, ci, mu = 1e-04, control = list(),
          method = if (is.null(grad)) "Nelder-Mead" else "BFGS", outer.iterations = 100,
          outer.eps = 1e-05, ..., hessian = FALSE) {
    if (!is.null(control$fnscale) && control$fnscale < 0)
        mu <- -mu
    R <- function(theta, theta.old, ...) {
        ui.theta <- ui %*% theta
        gi <- ui.theta - ci
        if (any(gi < -1e-12)){
            return(NaN)
        }else{
            gi[gi<0]=0
        }
        gi.old <- ui %*% theta.old - ci
        gi.old[gi.old<0]=0

        bar <- sum(gi.old * log(gi) - ui.theta)
        if (!is.finite(bar)){
            bar <- -1e10
#             bar <- -Inf
        }
        f(theta, ...) - mu * bar
    }
    dR <- function(theta, theta.old, ...) {
        ui.theta <- ui %*% theta
        gi <- drop(ui.theta - ci)
        gi.old <- drop(ui %*% theta.old - ci)
        dbar <- colSums(ui * gi.old/gi - ui)
        grad(theta, ...) - mu * dbar
    }
    if (any(ui %*% theta - ci <= -1e-12))
        stop("initial value is not in the interior of the feasible region")
    obj <- f(theta, ...)
    r <- R(theta, theta, ...)
    fun <- function(theta, ...) R(theta, theta.old, ...)
    gradient <- if (method == "SANN") {
        if (missing(grad))
            NULL
        else grad
    }
    else function(theta, ...) dR(theta, theta.old, ...)
    totCounts <- 0
    s.mu <- sign(mu)
    for (i in seq_len(outer.iterations)) {
        obj.old <- obj
        r.old <- r
        theta.old <- theta
        a <- optim(theta.old, fun, gradient, control = control,
                   method = method, hessian = hessian, ...)
        r <- a$value
        if (is.finite(r) && is.finite(r.old) && abs(r - r.old) <
                (0.001 + abs(r)) * outer.eps)
            break
        theta <- a$par
        totCounts <- totCounts + a$counts
        obj <- f(theta, ...)
        if (s.mu * obj > s.mu * obj.old)
            break
    }
    if (i == outer.iterations) {
        a$convergence <- 7
        a$message <- gettext("Barrier algorithm ran out of iterations and did not converge")
    }
    if (mu > 0 && obj > obj.old) {
        a$convergence <- 11
        a$message <- gettextf("Objective function increased at outer iteration %d",
                              i)
    }
    if (mu < 0 && obj < obj.old) {
        a$convergence <- 11
        a$message <- gettextf("Objective function decreased at outer iteration %d",
                              i)
    }
    a$outer.iterations <- i
    a$counts <- totCounts
    a$barrier.value <- a$value
    a$value <- f(a$par, ...)
    a$barrier.value <- a$barrier.value - a$value
    a
}

RELATIVE.CARDINALITY = function(la, ua, lb, ub, sim.params, E = 1e-6){
    MAX.I = 1000
    AscIterBinSearch <- function(A, value) {
        low = 1
        high = length(A)
#         if(any(is.na(A))){
#             browser()
#         }
        if ( A[[high]] < value )
            return (high+1)

        while ( low <= high ) {
            mid <- floor((low + high)/2)
            if ( A[mid] > value )
                high <- mid - 1
            else if ( A[mid] < value )
                low <- mid + 1
            else
                return(mid)
        }
        return(max(low, mid, high))
    }
    DescIterBinSearch <- function(A, value) {
        low = 1
        high = length(A)
        if ( A[[high]] > value )
            return (high+1)

#         if(any(is.na(A))){
#             browser()
#         }
        while ( low <= high ) {
            mid <- floor((low + high)/2)
            if ( A[mid] < value )
                high <- mid - 1
            else if ( A[mid] > value )
                low <- mid + 1
            else
                return(mid)
        }
        return(max(low, mid, high))
    }

    # Init
    n = length(la);

    l.pi = apply(cbind(la, ua, lb, ub), 1, sim.params$lp)
    u.pi = apply(cbind(la, ua, lb, ub), 1, sim.params$up)

    # Step 1
    l.order = order(l.pi, decreasing=F)
    u.order = order(u.pi, decreasing=T)

    l.pi = l.pi[l.order]
    u.pi = u.pi[u.order]

    l.la = la[l.order]
    l.ua = ua[l.order]
    l.lb = lb[l.order]
    l.ub = ub[l.order]

    u.la = la[u.order]
    u.ua = ua[u.order]
    u.lb = lb[u.order]
    u.ub = ub[u.order]

    l.binded = cbind(l.la, l.ua, l.lb, l.ub)
    u.binded = cbind(u.la, u.ua, u.lb, u.ub)

    # Step 2
    l.kp = floor(n / 2.4)
    u.kp = floor(n / 2.4)

    # Step 3
    l.x = c(l.ub[1:l.kp], l.lb[(l.kp+1):n])
    u.x = c(u.ub[1:u.kp], u.lb[(u.kp+1):n])

    l.m = sum(sapply(apply(cbind(l.la, l.x), 1, sim.params$tnorm), sim.params$f))
    u.m = sum(sapply(apply(cbind(u.ua, u.x), 1, sim.params$tnorm), sim.params$f))

    # Step 4
    l.M = sum(sapply(l.x, sim.params$f))
    u.M = sum(sapply(u.x, sim.params$f))

    # Step 5
    l.yp = l.m / l.M
    u.yp = u.m / u.M

    # Step 6 for lower bound
    l.i = 1
    repeat {
        # Step 7
        l.y = l.m / l.M

        # Step 8
        repeat {
            if(l.i > MAX.I){
                warning(paste("Too many iterations, interuped by watchdog l.i=", l.i))
                break
            }
            l.i = l.i + 1

            # Step 9
            l.k=l.kp
            # Step 10
            l.kp = AscIterBinSearch(l.pi, l.y) - 1

            # Step 11 - 13
            if(l.kp > 0){
                if (l.kp < n){
                    l.x = c(
                        apply(matrix(l.binded[1:l.kp,],ncol=4), 1, sim.params$lu, l.m, l.M),
                        l.lb[(l.kp+1):n]
                    )
                } else {
                    l.x =apply(l.binded, 1, sim.params$lu, l.m, l.M)
                }
            }else{
                l.x = l.lb
            }

            l.m = sum(sapply(apply(cbind(l.la, l.x), 1, sim.params$tnorm), sim.params$f))
            l.M = sum(sapply(l.x, sim.params$f))

            # Step 14
            l.y = l.m / l.M

            # Step 15
            if( l.kp == l.k ) {
                break
            }
        }

        # Step 16
        l.error = abs(l.yp - l.y)

        # Step 17
        l.yp = l.y

        # Step 18
        if( l.error < E | l.i > MAX.I){
            break
        }
    }
    #################################################################
    # Step 6 for upper bound
    u.i = 0
    repeat {
        # Step 7
        u.y = u.m / u.M

        # Step 8
        repeat {
            if(u.i > MAX.I){
                warning(paste("Too many iterations, interuped by watchdog u.i=", u.i))
                break
            }
            u.i = u.i + 1

            # Step 9
            u.k=u.kp
            # Step 10
            u.kp = DescIterBinSearch(u.pi, u.y) - 1

            # Step 11 - 13
            if(u.kp > 0){
                if (u.kp < n){
                    u.x = c(
                        apply(matrix(u.binded[1:u.kp,],ncol=4), 1, sim.params$uu, u.m, u.M),
                        u.lb[(u.kp+1):n]
                    )
                } else {
                    u.x =apply(u.binded, 1, sim.params$uu, u.m, u.M)
                }
            }else{
                u.x = u.lb
            }

            u.m = sum(sapply(apply(cbind(u.ua, u.x), 1, sim.params$tnorm), sim.params$f))
            u.M = sum(sapply(u.x, sim.params$f))

            # Step 14
            u.y = u.m / u.M

            # Step 15
            if( u.kp == u.k ) {
                break
            }
        }

        # Step 16
        u.error = abs(u.yp - u.y)

        # Step 17
        u.yp = u.y

        # Step 18
        if( u.error < E | u.i > MAX.I){
            break
        }
    }

    # Validate the result
    l.x = l.x[order(l.order)]
    if( any(l.x - lb < -1e-6) | any(l.x - ub > 1e-6)){
        print("Incorrect lower")
        print(lb)
        print(ub)
        print(l.x)
        warning("Incorrect resulting l.x")
    }
    u.x = u.x[order(u.order)]
    if( any(u.x - lb < -1e-6) | any(u.x -ub > 1e-6)){
        print("Incorrect upper")
        print(lb)
        print(ub)
        print(u.x)
        warning("Incorrect resulting u.x")
    }

    # Step 19
    return(list(ly=l.y, uy=u.y, lx=l.x, ux = u.x))
}

# numeric optimisation based RC implementation
RELATIVE.CARDINALITY.OPT = function(la, ua, lb, ub, sim.params, E = 1e-6){
    #   tnorm, f, lu, uu, lp, up
    n=length(la)
    ui = matrix(c(sapply(1:n, function(x){c(diag(n)[,x], -diag(n)[,x])}, simplify = T)), nrow=2*n, ncol=n, byrow=T)
    ci = c(rbind(lb, -ub))

    l.results = lapply(seq(0,1,0.5), function(alpha){
        theta = alpha * lb + (1 - alpha) * ub
        myConstrOptim(theta, function(x){
            m = sum(sapply(apply(cbind(la, x), 1, sim.params$tnorm), sim.params$f))
            M = sum(sapply(x, sim.params$f))
            return(m/M)
        }, NULL, ui, ci, outer.eps = E, outer.iterations = 100)
    })
    l = which.min(lapply(l.results, '[[','value'))
    l.result=l.results[[l]]

    u.results = lapply(seq(0,1,0.5), function(alpha){
        theta = alpha * lb + (1 - alpha) * ub
        myConstrOptim(theta, function(x){
            m = sum(sapply(apply(cbind(ua, x), 1, sim.params$tnorm), sim.params$f))
            M = sum(sapply(x, sim.params$f))
            return(m/M)
        }, NULL, ui, ci, outer.eps = E, outer.iterations = 100, control = list(fnscale = -1))
    })
    u = which.max(lapply(u.results, '[[','value'))
    u.result=u.results[[u]]
    return(list(ly = l.result$value, uy = u.result$value, lx = l.result$par, ux = u.result$par))
}


GEN.SIM.JACCARD = function(sim.params, optimisation = F){
    force(sim.params)
    if(optimisation){
        RC = RELATIVE.CARDINALITY.OPT
    }else{
        RC = RELATIVE.CARDINALITY
    }
    return(function(iA, iB){
        #calculate AB and AvB
        iAiB = list(
            lower = apply(cbind(iA$lower, iB$lower), 1, sim.params$tnorm),
            upper = apply(cbind(iA$upper, iB$upper), 1, sim.params$tnorm))
#         print("iAiB")
#         print(iAiB)
        tconorm  = function(x){
            return(1-sim.params$tnorm(c(1-x[[1]],1-x[[2]])))
        }
        iAviB = list(
            lower = apply(cbind(iA$lower, iB$lower), 1, tconorm),
            upper = apply(cbind(iA$upper, iB$upper), 1, tconorm))
#         print("iAviB")
#         print(iAviB)

        result = RC(iAiB$lower, iAiB$upper, iAviB$lower, iAviB$upper, sim.params)
        return(result)
    })
}

GEN.SIM.PARAMS.OPTIMISATION = function(tnorm.p, f.p){
    compLimit = function(fn, b0){
        eps = 1e-4
        if(b0 + eps<1){
            x = seq(b0 + eps, b0+1e-6, by=-1e-6)
        }else{
            x = seq(b0 - eps, b0-1e-6, by=1e-6)
        }
        return(mean(sapply(x, fn), na.rm=T))
    }
    force(tnorm.p);force(f.p)
    return(list(
        tnorm = function(x){
            return(tnorm.p(x[[1]], x[[2]]))
        },
        f = f.p,
        lp = function(x){
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]

            fn = function(b){
                return(min(1, (f.p(tnorm.p(la, b))-f.p(tnorm.p(la, lb)))/(f.p(b)-f.p(lb))))
            }

            if(abs(ub - lb) <1e-4){
                return(compLimit(fn, mean(lb, ub)))
            }else{
                theta = (lb + ub)/2
                res = optim(theta, fn, method="Brent", lower = lb, upper = ub)
#                 if(is.nan(res$value)){
#                     browser()
#                 }
#                 if(is.na(res$value)){
#                     browser()
#                 }
                return(res$value)
            }
        },
        up = function(x){
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]

            fn = function(b){
                return(min(1, (f.p(tnorm.p(ua, b))-f.p(tnorm.p(ua, lb)))/(f.p(b)-f.p(lb))))
            }

            if(abs(ub - lb) <1e-4){
                return(compLimit(fn, mean(lb, ub)))
            }else{
                theta = (lb + ub)/2
                res = optim(theta, fn, method="Brent", lower = lb, upper = ub, control = list(fnscale = -1))
#                 if(is.nan(res$value)){
#                     browser()
#                 }
#                 if(is.na(res$value)){
#                     browser()
#                 }
                return(res$value)
            }
        },
        lu = function(x, m, M){
#             print("lu")
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
#             ui = matrix(c(1, -1), nrow=2, ncol=1, byrow=T)
#             ci = c(lb, -ub)
            if(lb == ub){
                return(lb)
            }else{
                theta = (lb + ub)/2
                res = optim(theta, function(b){
                    return((m + f.p(tnorm.p(la, b)))/(M + f.p(b)))
                }, method="Brent", lower = lb, upper = ub)
    #             res = myConstrOptim(theta, function(b){
    #                 return((m + f.p(tnorm.p(la, b)))/(M + f.p(b)))
    #             }, NULL, ui, ci, outer.eps = 1e-6, method="Brent")
                return(res$par)
            }
        },
        uu = function(x, m, M){
#             print("uu")
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
#             ui = matrix(c(1, -1), nrow=2, ncol=1, byrow=T)
#             ci = c(lb, -ub)
            if(lb == ub){
                return(lb)
            }else{
                theta = (lb + ub)/2
                res = optim(theta, function(b){
                    return((m + f.p(tnorm.p(ua, b)))/(M + f.p(b)))
                }, method="Brent", lower = lb, upper = ub, control = list(fnscale = -1))
    #             res = myConstrOptim(theta, function(b){
    #                 return((m + f.p(tnorm.p(ua, b)))/(M + f.p(b)))
    #             }, NULL, ui, ci, outer.eps = 1e-6, method="Brent", control = list(fnscale = -1))
                return(res$par)
            }
        }))
}

GEN.TNORM.SS = function(p){
    force(p)
    if(p<0){
        return(function(a,b){
            return((a^p + b^p - 1)^(1/p))
        })
    } else if (p==0){
        return(function(a,b){
            return(a * b)
        })
    } else {
        return(function(a,b){
            return(max(0, a^p + b^p - 1)^(1/p))
        })
    }
}
SIM.PARAMS.OPT.SS.05.ID = GEN.SIM.PARAMS.OPTIMISATION(GEN.TNORM.SS(0.5), identity)
SIM.PARAMS.OPT.SS.m05.ID = GEN.SIM.PARAMS.OPTIMISATION(GEN.TNORM.SS(-0.5), identity)
SIM.PARAMS.OPT.SS.m1.ID = GEN.SIM.PARAMS.OPTIMISATION(GEN.TNORM.SS(-1), identity)
SIM.PARAMS.OPT.SS.m2.ID = GEN.SIM.PARAMS.OPTIMISATION(GEN.TNORM.SS(-2), identity)
SIM.PARAMS.OPT.SS.m5.ID = GEN.SIM.PARAMS.OPTIMISATION(GEN.TNORM.SS(-5), identity)
SIM.PARAMS.OPT.SS.m25.ID = GEN.SIM.PARAMS.OPTIMISATION(GEN.TNORM.SS(-25), identity)

SIM.PARAMS.EXACT.MINIMUM.ID = list(
    tnorm = function(x){return(min(x[[1]], x[[2]]))},
    f=identity,
    lp = function(x) {
        la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
        if(la<lb) {
            return(0)
        } else if (lb <= la & la < ub){
            return((la-lb)/(ub - lb))
        } else {
            return(1)
        }
    },
    up = function(x) {
        la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
        if(ua<lb) {
            return(0)
        } else if (lb <= ua & ua < ub){
            return(1)
        } else {
            return(1)
        }
    },
    lu = function(x, m, M) {
        return(x[[4]])
    },
    uu = function(x, m, M) {
        la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
        if(ua<lb) {
            return(lb)
        } else if (lb <= ua & ua < ub){
            return(ua)
        } else {
            return(ub)
        }
    }
)

SIM.PARAMS.EXACT.PRODUCT.ID = list(
    tnorm = function(x){return(x[[1]] * x[[2]])},
    f=identity,
    lp = function(x) {
        return(x[[1]])
    },
    up = function(x) {
        return(x[[2]])
    },
    lu = function(x, m, M) {
        return(x[[4]])
    },
    uu = function(x, m, M) {
        return(x[[4]])
    }
)

GEN.SIM.PARAMS.EXACT.SW.ID = function(lambda){
    if(!(lambda > -1 & lambda<=0)){
        stop("SW tnorm has not u-constant property for this lambda")
    }
    return(list(
        tnorm = function(x){
            return(max(0, (x[[1]]+x[[2]]-1+lambda*x[[1]]*x[[2]])/(1+lambda)))
        },
        f = identity,
        lp = function(x, tnorm, f) {
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
            if( (1-la)/(1+lambda*la) < lb ){
                return((1+lambda*la)/(1+lambda))
            } else if (lb <= (1-la)/(1+lambda*la) & (1-la)/(1+lambda*la) <ub) {
                return(0)
            }else{
                return(0)
            }
        },
        up = function(x, tnorm, f) {
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
            if( (1-ua)/(1+lambda*ua) < lb ){
                return((1+lambda*ua)/(1+lambda))
            } else if (lb <= (1-ua)/(1+lambda*ua) & (1-ua)/(1+lambda*ua) <ub) {
                return((max(0, (ua+ub-1+lambda*ua*ub)/(1+lambda)))/(ub-lb))
            }else{
                return(0)
            }
        },
        lu = function(x, m, M, tnorm, f) {
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
            if( (1-la)/(1+lambda*la) < lb ){
                return(ub)
            } else if (lb <= (1-la)/(1+lambda*la) & (1-la)/(1+lambda*la) <ub) {
                val = (1-la)/(1+lambda*la)
                return(min(ub, max(lb, val)))
            }else{
                return(ub)
            }
        },
        uu = function(x, m, M, tnorm, f) {
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
            if( (1-ua)/(1+lambda*ua) < lb ){
                return(ub)
            } else if (lb <= (1-ua)/(1+lambda*ua) & (1-ua)/(1+lambda*ua) <ub) {
                return(ub)
            }else{
                return(lb)
            }
        }
    ))
}
SIM.PARAMS.EXACT.SW.m0.5.ID = GEN.SIM.PARAMS.EXACT.SW.ID(-0.5)

GEN.SIM.PARAMS.EXACT.SS.ID = function(lambda){
    if(!(lambda >= 1)){
        stop("SS tnorm has not u-constant property for this lambda")
    }
    return(list(
        tnorm = function(x){
            return(max(0, (x[[1]]^lambda+x[[2]]^lambda-1))^(1/lambda))
        },
        f = identity,
        lp = function(x, tnorm, f) {
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
            if( (1-la^lambda)^(1/lambda) < lb ){
                return(1)
            } else if (lb <= (1-la^lambda)^(1/lambda) & (1-la^lambda)^(1/lambda) <ub) {
                return(0)
            }else{
                return(0)
            }
        },
        up = function(x, tnorm, f) {
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
            if( (1-ua^lambda)^(1/lambda) < lb ){
                return(1)
            } else if (lb <= (1-ua^lambda)^(1/lambda) & (1-ua^lambda)^(1/lambda) <ub) {
                return(max(0, ((ua^lambda + ub^lambda -1)^(1/lambda))/(ub - lb) ))
            }else{
                return(0)
            }
        },
        lu = function(x, m, M, tnorm, f) {
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
            if( (1-la^lambda)^(1/lambda) < lb ){
                return(lb)
            } else if (lb <= (1-la^lambda)^(1/lambda) & (1-la^lambda)^(1/lambda) <ub) {
                return((1-la^lambda)^(1/lambda))
            }else{
                return(ub)
            }
        },
        uu = function(x, m, M, tnorm, f) {
            la = x[[1]]; ua = x[[2]]; lb = x[[3]]; ub = x[[4]]
            if( (1-ua^lambda)^(1/lambda) < lb ){
                return(ub)
            } else if (lb <= (1-ua^lambda)^(1/lambda) & (1-ua^lambda)^(1/lambda) <ub) {
                return(ub)
            }else{
                return(lb)
            }
        }
    ))
}
SIM.PARAMS.EXACT.LUK.ID = GEN.SIM.PARAMS.EXACT.SS.ID(1)
SIM.PARAMS.EXACT.SS.2.ID = GEN.SIM.PARAMS.EXACT.SS.ID(2)
SIM.PARAMS.EXACT.SS.5.ID = GEN.SIM.PARAMS.EXACT.SS.ID(5)
SIM.PARAMS.EXACT.SS.25.ID = GEN.SIM.PARAMS.EXACT.SS.ID(25)

SIM.PARAMS = list(SIM.PARAMS.EXACT.MINIMUM.ID,
               SIM.PARAMS.OPT.SS.m25.ID,
               SIM.PARAMS.OPT.SS.m5.ID,
               SIM.PARAMS.OPT.SS.m2.ID,
               SIM.PARAMS.OPT.SS.m1.ID,
               SIM.PARAMS.OPT.SS.m05.ID,
               SIM.PARAMS.EXACT.PRODUCT.ID,
               SIM.PARAMS.OPT.SS.05.ID,
               SIM.PARAMS.EXACT.LUK.ID,
               SIM.PARAMS.EXACT.SS.2.ID,
               SIM.PARAMS.EXACT.SS.5.ID,
               SIM.PARAMS.EXACT.SS.25.ID)