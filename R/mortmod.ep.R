mortmod.ep <- function(recal, entry.par, prev, child.art=NULL, adult.art=NULL, region=1, sex=1, opt=TRUE, entry.nax=NULL){
  for(f in recal) load(file=f)
  which.user.par <- c(which(ages==user.par[1]),which(ages==user.par[2]))
  
  lt.mx <- function(nmx, age = ages, nax=entry.nax){
    n <- c(diff(age), 999)
    if(!is.null(nax)){
      nax <- nax
    } else {
      nax <- 0.5*n
    }   
    nqx <- (n * nmx)/(1 + (n - nax) * nmx)
    nqx <- c(nqx[-(length(nqx))], 1)
    for (i in 1:length(nqx)) {
      if (nqx[i] > 1) 
        nqx[i] <- 1
    }
    nage <- length(age)
    nqx <- round(nqx, 4)
    npx <- 1 - nqx
    l0 = 1e+05
    lx <- round(cumprod(c(l0, npx)))
    ndx <- -diff(lx)
    lxpn <- lx[-1]
    nLx <- n * lxpn + ndx * nax
    Tx <- c(rev(cumsum(rev(nLx[-length(nLx)]))),0)
    lx <- lx[1:length(age)]
    ex <- Tx/lx
    lt <- cbind(Age = age, nax = c(round(nax[-length(nax)], 3),NA), nmx = round(nmx,4), nqx = round(nqx, 4), npx = round(npx, 4), ndx = ndx, 
                lx = lx, nLx = c(round(nLx[-length(nLx)]),NA), Tx = c(round(Tx[-length(Tx)]),NA), ex = c(round(ex[-length(ex)],2),NA))
    lt <- lt[lt[, 6] != 0, ]
    e0 <- lt[1, 10]
    if(n[4]>1){
      lt.45q15 <- 1 - (lx[14]/lx[5])
      lt.5q0 <- 1 - (lx[3]/lx[1])
    }
    if(n[4]==1){
      lt.45q15 <- 1 - (lx[61]/lx[16])
      lt.5q0 <- 1 - (lx[6]/lx[1])
    }
    return(list(e0 = e0, lt.5q0 = lt.5q0, lt.45q15 = lt.45q15, 
                lt = lt))
  }  
  
if(region==1 & sex==0){ # Africa, male
	intercept <- median(svd.coeffs.xp.m[africa.nums,1])
	b1.m <- predict.lm(co1.epmod.a.m, newdata=data.frame(ep.a.m=entry.par, prev.a=prev))
	b2.m <- predict.lm(co2.epmod.a.m, newdata=data.frame(ep.a.m=entry.par, prev.a=prev))
	b3.m <- predict.lm(co3.epmod.a.m, newdata=data.frame(ep.a.m=entry.par, prev.a=prev))
  
	if(!is.null(adult.art)){
	  b2.m <- predict.lm(co2.epmod2.a.m, newdata=data.frame(ep.a.m=entry.par, prev.a=prev, aart.a=adult.art))
	  b3.m <- predict.lm(co3.epmod2.a.m, newdata=data.frame(ep.a.m=entry.par, prev.a=prev, aart.a=adult.art)) 
	}
	
	if(!is.null(adult.art)&!is.null(child.art)){
	  b2.m <- predict.lm(co2.epmod3.a.m, newdata=data.frame(ep.a.m=entry.par, prev.a=prev, aart.a=adult.art, cart.a=child.art))
	  b3.m <- predict.lm(co3.epmod3.a.m, newdata=data.frame(ep.a.m=entry.par, prev.a=prev, aart.a=adult.art, cart.a=child.art)) 
	}
	
	if(opt==FALSE){
	out.mort <- intercept + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]} # if
	
	
## optimize reduction or addition to 1q0 and 4q1 to make it match
if(opt==TRUE){
  out.mort.func <- function(intercept.alter){
    out.mort <- intercept.alter + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]
    lt.out <- lt.mx(nmx=exp(out.mort), age=ages)
    lt.ep <- unname(1-(lt.out$lt[which.user.par[2],7]/lt.out$lt[which.user.par[1],7]))
    ep.diff <- abs(entry.par-lt.ep)
    return(ep.diff)	
  } # function to be optimized
  
  intercept.opt <- optimize(f=out.mort.func, interval=c(-5,5))$minimum
  
  out.mort <- intercept.opt + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]} # if # if

return(exp(out.mort))
}

if(region==1 & sex==1){ # Africa, female
	intercept <- median(svd.coeffs.xp.f[africa.nums,1])
	b1.f <- predict.lm(co1.epmod.a.f, newdata=data.frame(ep.a.f=entry.par, prev.a=prev))
	b2.f <- predict.lm(co2.epmod.a.f, newdata=data.frame(ep.a.f=entry.par, prev.a=prev))
	b3.f <- predict.lm(co3.epmod.a.f, newdata=data.frame(ep.a.f=entry.par, prev.a=prev))
  
	if(!is.null(adult.art)){
	  b2.f <- predict.lm(co2.epmod2.a.f, newdata=data.frame(ep.a.f=entry.par, prev.a=prev, aart.a=adult.art))
	  b3.f <- predict.lm(co3.epmod2.a.f, newdata=data.frame(ep.a.f=entry.par, prev.a=prev, aart.a=adult.art)) 
	}
	
	if(!is.null(adult.art)&!is.null(child.art)){
	  b2.f <- predict.lm(co2.epmod3.a.f, newdata=data.frame(ep.a.f=entry.par, prev.a=prev, aart.a=adult.art, cart.a=child.art))
	  b3.f <- predict.lm(co3.epmod3.a.f, newdata=data.frame(ep.a.f=entry.par, prev.a=prev, aart.a=adult.art, cart.a=child.art)) 
	}
	
	if(opt==FALSE){
	out.mort <- intercept + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]} # if
		
	## optimize reduction or addition to 1q0 and 4q1 to make it match
	if(opt==TRUE){
	  out.mort.func <- function(intercept.alter){
	    out.mort <- intercept.alter + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]
	    lt.out <- lt.mx(nmx=exp(out.mort), age=ages)
	    lt.ep <- unname(1-(lt.out$lt[which.user.par[2],7]/lt.out$lt[which.user.par[1],7]))
	    ep.diff <- abs(entry.par-lt.ep)
	    return(ep.diff)	
	  } # function to be optimized
	  
	  intercept.opt <- optimize(f=out.mort.func, interval=c(-5,5))$minimum
	  
	  out.mort <- intercept.opt + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]} # if # if

return(exp(out.mort))
}

if(region==0 & sex==0){ # Non-Africa, male
	intercept <- median(svd.coeffs.xp.m[la.nums,1])
	b1.m <- predict.lm(co1.epmod.na.m, newdata=data.frame(ep.na.m=entry.par))
	b2.m <- predict.lm(co2.epmod.na.m, newdata=data.frame(ep.na.m=entry.par, prev.na=prev))
	b3.m <- predict.lm(co3.epmod.na.m, newdata=data.frame(ep.na.m=entry.par, prev.na=prev))
  
	if(!is.null(adult.art)){
	  b2.m <- predict.lm(co2.epmod2.na.m, newdata=data.frame(ep.na.m=entry.par, prev.na=prev, aart.na=adult.art))
	  b3.m <- predict.lm(co3.epmod2.na.m, newdata=data.frame(ep.na.m=entry.par, prev.na=prev, aart.na=adult.art)) 
	}
	
	if(!is.null(adult.art)&!is.null(child.art)){
	  b2.m <- predict.lm(co2.epmod3.na.m, newdata=data.frame(ep.na.m=entry.par, prev.na=prev, aart.na=adult.art, cart.na=child.art))
	  b3.m <- predict.lm(co3.epmod3.na.m, newdata=data.frame(ep.na.m=entry.par, prev.na=prev, aart.na=adult.art, cart.na=child.art)) 
	}
	
	if(opt==FALSE){
	out.mort <- intercept + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]} # if
	
	## optimize reduction or addition to 1q0 and 4q1 to make it match
	if(opt==TRUE){
	  out.mort.func <- function(intercept.alter){
	    out.mort <- intercept.alter + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]
	    lt.out <- lt.mx(nmx=exp(out.mort), age=ages)
	    lt.ep <- unname(1-(lt.out$lt[which.user.par[2],7]/lt.out$lt[which.user.par[1],7]))
	    ep.diff <- abs(entry.par-lt.ep)
	    return(ep.diff)	
	  } # function to be optimized
	  
	  intercept.opt <- optimize(f=out.mort.func, interval=c(-5,5))$minimum
	  
	  out.mort <- intercept.opt + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]}# if
  
return(exp(out.mort))
}

if(region==0 & sex==1){ # Non-Africa, female
	intercept <- median(svd.coeffs.xp.f[la.nums,1])
	b1.f <- predict.lm(co1.epmod.na.f, newdata=data.frame(ep.na.f=entry.par))
	b2.f <- predict.lm(co2.epmod.na.f, newdata=data.frame(ep.na.f=entry.par, prev.na=prev))
	b3.f <- predict.lm(co3.epmod.na.f, newdata=data.frame(ep.na.f=entry.par, prev.na=prev))
  
	if(!is.null(adult.art)){
	  b2.f <- predict.lm(co2.epmod2.na.f, newdata=data.frame(ep.na.f=entry.par, prev.na=prev, aart.na=adult.art))
	  b3.f <- predict.lm(co3.epmod2.na.f, newdata=data.frame(ep.na.f=entry.par, prev.na=prev, aart.na=adult.art)) 
	}
	
	if(!is.null(adult.art)&!is.null(child.art)){
	  b2.f <- predict.lm(co2.epmod3.na.f, newdata=data.frame(ep.na.f=entry.par, prev.na=prev, aart.na=adult.art, cart.na=child.art))
	  b3.f <- predict.lm(co3.epmod3.na.f, newdata=data.frame(ep.na.f=entry.par, prev.na=prev, aart.na=adult.art, cart.na=child.art)) 
	}
	
	if(opt==FALSE){
	  out.mort <- intercept + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]} # if
	
	## optimize reduction or addition to 1q0 and 4q1 to make it match
	if(opt==TRUE){
	  out.mort.func <- function(intercept.alter){
	    out.mort <- intercept.alter + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]
	    lt.out <- lt.mx(nmx=exp(out.mort), age=ages)
	    lt.ep <- unname(1-(lt.out$lt[which.user.par[2],7]/lt.out$lt[which.user.par[1],7]))
	    ep.diff <- abs(entry.par-lt.ep)
	    return(ep.diff)	
	  } # function to be optimized
	  
	  intercept.opt <- optimize(f=out.mort.func, interval=c(-5,5))$minimum
	  
	  out.mort <- intercept.opt + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]} # if # if
	
	return(exp(out.mort))
}
}
