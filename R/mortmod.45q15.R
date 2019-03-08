mortmod.45q15 <- function(child.mort, adult.mort, prev, child.art=NULL, adult.art=NULL, region=1, sex=1, opt=TRUE, recal=NULL){
  if(!is.null(recal)){for(f in recal) load(file=f)} 
  
  lt.mx <- function (nmx, sex="female", age = ages, nax=NULL){
    n <- c(diff(age), 999)
    if(is.null(nax)){
      nax <- 0.5*n
      if(n[2]==4){
        if(sex == "male"){
          if (nmx[1] >= 0.107) {
            nax[1] <- 0.33
            nax[2] <- 1.352
          } else {
            nax[1] <- 0.045 + 2.684 * nmx[1]
            nax[2] <- 1.651 - 2.816 * nmx[1]
          } # else
        } # if (sex == "male")
        if(sex == "female"){
          if (nmx[1] >= 0.107) {
            nax[1] <- 0.35
            nax[2] <- 1.361
          } else {
            nax[1] <- 0.053 + 2.8 * nmx[1]
            nax[2] <- 1.522 - 1.518 * nmx[1]
          } # else
        } #if (sex == "female")
      } #if(n[2]==4)
    } #if(is.null(nax))
    
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
    lt.45q15 <- 1 - (lx[which(age==60)]/lx[which(age==15)])
    lt.5q0 <- 1 - (lx[which(age==5)]/lx[which(age==0)])
    
    return(list(e0 = e0, lt.5q0 = lt.5q0, lt.45q15 = lt.45q15, 
                lt = lt))
  }
  
if(region==1 & sex==0){ # Africa, male
	intercept <- median(svd.coeffs.xp.m[africa.nums,1])
	b1.m <- predict.lm(co1.45q15mod.a.m, newdata=data.frame(cmort.a.m=child.mort, amort.a.m=adult.mort))
	b2.m <- predict.lm(co2.45q15mod.a.m, newdata=data.frame(cmort.a.m=child.mort, amort.a.m=adult.mort, prev.a=prev))
	b3.m <- predict.lm(co3.45q15mod.a.m, newdata=data.frame(cmort.a.m=child.mort, amort.a.m=adult.mort, prev.a=prev))
  
	if(!is.null(adult.art)){
	  b2.m <- predict.lm(co2.45q15mod2.a.m, newdata=data.frame(cmort.a.m=child.mort, amort.a.m=adult.mort, prev.a=prev, aart.a=adult.art))
	  b3.m <- predict.lm(co3.45q15mod2.a.m, newdata=data.frame(cmort.a.m=child.mort, amort.a.m=adult.mort, prev.a=prev, aart.a=adult.art)) 
	}
  
	if(!is.null(adult.art)&!is.null(child.art)){
	  b2.m <- predict.lm(co2.45q15mod3.a.m, newdata=data.frame(cmort.a.m=child.mort, amort.a.m=adult.mort, prev.a=prev, aart.a=adult.art, cart.a=child.art))
	  b3.m <- predict.lm(co3.45q15mod3.a.m, newdata=data.frame(cmort.a.m=child.mort, amort.a.m=adult.mort, prev.a=prev, aart.a=adult.art, cart.a=child.art)) 
	}
	
	if(opt==FALSE){
	out.mort <- intercept + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]} # if
	
# optimize the intercept when predicting the mortality rates from the weights
	if(opt==TRUE){
		out.mort.func <- function(intercept.alter){
		out.mort <- intercept.alter + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]
    #lt.out <- LifeTableMx(mx=exp(out.mort), sex="Male")
    lt.out <- lt.mx(nmx=exp(out.mort), sex="male", age=ages)
    #lt.out.45q15 <- 1-lt.out$lx[14]/lt.out$lx[5]
    lt.out.45q15 <- lt.out$lt.45q15
    amort.diff <- abs(adult.mort-lt.out.45q15)     
	return(amort.diff)	
	} # function to be optimized
	
	intercept.opt <- optimize(f=out.mort.func, interval=c(-2,2))$minimum
	
	out.mort <- intercept.opt + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]} # if

	
## optimize reduction or addition to 1q0 and 4q1 to make it match
if(opt==TRUE){
	out.mort.start <- out.mort
	#lt.out <- LifeTableMx(mx=exp(out.mort.start), sex="Male")
	lt.out <- lt.mx(nmx=exp(out.mort.start), sex="male", age=ages)
	
	diff.weight.opt <- log(1-child.mort)/(log(1-lt.out$lt[1,4])+log(1-lt.out$lt[2,4]))
	qx.child.new.opt <- 1-(c((1-lt.out$lt[1,4])^diff.weight.opt, (1-lt.out$lt[2,4])^diff.weight.opt))
	cmort.opt <- 1-((1-qx.child.new.opt[1])*(1-qx.child.new.opt[2]))
	
	axs <- c(1-lt.out$lt[1,2], 4-lt.out$lt[2,2])
	mx1 <- qx.child.new.opt[1]/(1-axs[1]*qx.child.new.opt[1])
	mx2 <- qx.child.new.opt[2]/(4-axs[2]*qx.child.new.opt[2])
	
	mx.child.new.opt <- c(mx1, mx2)
	out.mort.new.opt <- unname(c(log(mx.child.new.opt), out.mort.start[3:length(out.mort.start)]))
	if(out.mort.new.opt[3]>out.mort.new.opt[2]){
	  # y <- out.mort.new.opt[1:2]
	  # x <- c(1,2)
	  # child.mod.lm <- lm(y~x)
	  # new.5q5 <- predict.lm(child.mod.lm, newdata=data.frame(x = 3))
	  new.5m5 <- approx(x=c(2,4),y=out.mort.new.opt[c(2,4)],xout=3)$y
	  out.mort.new.opt[3] <- new.5m5
	}
	#lt.out.new.opt <- LifeTableMx(mx=exp(out.mort.new.opt), sex="Male")
	lt.out.new.opt <- lt.mx(nmx=exp(out.mort.new.opt), sex="male", age=ages, nax=NULL)
	out.mort <- out.mort.new.opt
} # if


	return(exp(out.mort))
}

if(region==1 & sex==1){ # Africa, female
  intercept <- median(svd.coeffs.xp.f[africa.nums,1])
  b1.f <- predict.lm(co1.45q15mod.a.f, newdata=data.frame(cmort.a.f=child.mort, amort.a.f=adult.mort))
  b2.f <- predict.lm(co2.45q15mod.a.f, newdata=data.frame(cmort.a.f=child.mort, amort.a.f=adult.mort, prev.a=prev))
  b3.f <- predict.lm(co3.45q15mod.a.f, newdata=data.frame(cmort.a.f=child.mort, amort.a.f=adult.mort, prev.a=prev))

  if(!is.null(adult.art)){
    b2.f <- predict.lm(co2.45q15mod2.a.f, newdata=data.frame(cmort.a.f=child.mort, amort.a.f=adult.mort, prev.a=prev, aart.a=adult.art))
    b3.f <- predict.lm(co3.45q15mod2.a.f, newdata=data.frame(cmort.a.f=child.mort, amort.a.f=adult.mort, prev.a=prev, aart.a=adult.art)) 
  }
  
  if(!is.null(adult.art)&!is.null(child.art)){
    b2.f <- predict.lm(co2.45q15mod3.a.f, newdata=data.frame(cmort.a.f=child.mort, amort.a.f=adult.mort, prev.a=prev, aart.a=adult.art, cart.a=child.art))
    b3.f <- predict.lm(co3.45q15mod3.a.f, newdata=data.frame(cmort.a.f=child.mort, amort.a.f=adult.mort, prev.a=prev, aart.a=adult.art, cart.a=child.art)) 
  }
  
  if(opt==FALSE){
    out.mort <- intercept + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]} # if
  
  # optimize the intercept when predicting the mortality rates from the weights
  if(opt==TRUE){
    out.mort.func <- function(intercept.alter){
      out.mort <- intercept.alter + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]
      #lt.out <- LifeTableMx(mx=exp(out.mort), sex="Female")
      lt.out <- lt.mx(nmx=exp(out.mort), sex="female", age=ages)
      #lt.out.45q15 <- 1-lt.out$lx[14]/lt.out$lx[5]
      lt.out.45q15 <- lt.out$lt.45q15
      amort.diff <- abs(adult.mort-lt.out.45q15)     
      return(amort.diff)	
    } # function to be optimized
    
    intercept.opt <- optimize(f=out.mort.func, interval=c(-2,2))$minimum
    
    out.mort <- intercept.opt + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]} # if
  
  
  ## optimize reduction or addition to 1q0 and 4q1 to make it match
  if(opt==TRUE){
    out.mort.start <- out.mort
    #lt.out <- LifeTableMx(mx=exp(out.mort.start), sex="Female")
    lt.out <- lt.mx(nmx=exp(out.mort.start), sex="female", age=ages)
    
    diff.weight.opt <- log(1-child.mort)/(log(1-lt.out$lt[1,4])+log(1-lt.out$lt[2,4]))
    qx.child.new.opt <- 1-(c((1-lt.out$lt[1,4])^diff.weight.opt, (1-lt.out$lt[2,4])^diff.weight.opt))
    cmort.opt <- 1-((1-qx.child.new.opt[1])*(1-qx.child.new.opt[2]))
    
    axs <- c(1-lt.out$lt[1,2], 4-lt.out$lt[2,2])
    mx1 <- qx.child.new.opt[1]/(1-axs[1]*qx.child.new.opt[1])
    mx2 <- qx.child.new.opt[2]/(4-axs[2]*qx.child.new.opt[2])
    
    mx.child.new.opt <- c(mx1, mx2)
    out.mort.new.opt <- unname(c(log(mx.child.new.opt), out.mort.start[3:length(out.mort.start)]))
    if(out.mort.new.opt[3]>out.mort.new.opt[2]){
      # y <- out.mort.new.opt[1:2]
      # x <- c(1,2)
      # child.mod.lm <- lm(y~x)
      # new.5q5 <- predict.lm(child.mod.lm, newdata=data.frame(x = 3))
      new.5m5 <- approx(x=c(2,4),y=out.mort.new.opt[c(2,4)],xout=3)$y
      out.mort.new.opt[3] <- new.5m5
    }
    #lt.out.new.opt <- LifeTableMx(mx=exp(out.mort.new.opt), sex="Female")
    lt.out.new.opt <- lt.mx(nmx=exp(out.mort.new.opt), sex="female", age=ages, nax=NULL)
    out.mort <- out.mort.new.opt
  } # if
  
  
  return(exp(out.mort))
}

if(region==0 & sex==0){ # Non-Africa, male
  intercept <- median(svd.coeffs.xp.m[la.nums,1])
  b1.m <- predict.lm(co1.45q15mod.na.m, newdata=data.frame(cmort.na.m=child.mort, amort.na.m=adult.mort))
  b2.m <- predict.lm(co2.45q15mod.na.m, newdata=data.frame(cmort.na.m=child.mort, amort.na.m=adult.mort, prev.na=prev))
  b3.m <- predict.lm(co3.45q15mod.na.m, newdata=data.frame(cmort.na.m=child.mort, amort.na.m=adult.mort, prev.na=prev))
  
  if(!is.null(adult.art)){
    b2.m <- predict.lm(co2.45q15mod2.na.m, newdata=data.frame(cmort.na.m=child.mort, amort.na.m=adult.mort, prev.na=prev, aart.na=adult.art))
    b3.m <- predict.lm(co3.45q15mod2.na.m, newdata=data.frame(cmort.na.m=child.mort, amort.na.m=adult.mort, prev.na=prev, aart.na=adult.art)) 
  }
  
  if(!is.null(adult.art)&!is.null(child.art)){
    b2.m <- predict.lm(co2.45q15mod3.na.m, newdata=data.frame(cmort.na.m=child.mort, amort.na.m=adult.mort, prev.na=prev, aart.na=adult.art, cart.na=child.art))
    b3.m <- predict.lm(co3.45q15mod3.na.m, newdata=data.frame(cmort.na.m=child.mort, amort.na.m=adult.mort, prev.na=prev, aart.na=adult.art, cart.na=child.art)) 
  }
  
  if(opt==FALSE){
    out.mort <- intercept + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]} # if
  
  # optimize the intercept when predicting the mortality rates from the weights
  if(opt==TRUE){
    out.mort.func <- function(intercept.alter){
      out.mort <- intercept.alter + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]
      #lt.out <- LifeTableMx(mx=exp(out.mort), sex="Male")
      lt.out <- lt.mx(nmx=exp(out.mort), sex="male", age=ages)
      #lt.out.45q15 <- 1-lt.out$lx[14]/lt.out$lx[5]
      lt.out.45q15 <- lt.out$lt.45q15
      amort.diff <- abs(adult.mort-lt.out.45q15)     
      return(amort.diff)	
    } # function to be optimized
    
    intercept.opt <- optimize(f=out.mort.func, interval=c(-2,2))$minimum
    
    out.mort <- intercept.opt + b1.m*Mx.svd.scores.m[,1] + b2.m*Mx.svd.scores.m[,2] + b3.m*Mx.svd.scores.m[,3]} # if
  
  
  ## optimize reduction or addition to 1q0 and 4q1 to make it match
  if(opt==TRUE){
    out.mort.start <- out.mort
    #lt.out <- LifeTableMx(mx=exp(out.mort.start), sex="Male")
    lt.out <- lt.mx(nmx=exp(out.mort.start), sex="male", age=ages)
    
    diff.weight.opt <- log(1-child.mort)/(log(1-lt.out$lt[1,4])+log(1-lt.out$lt[2,4]))
    qx.child.new.opt <- 1-(c((1-lt.out$lt[1,4])^diff.weight.opt, (1-lt.out$lt[2,4])^diff.weight.opt))
    cmort.opt <- 1-((1-qx.child.new.opt[1])*(1-qx.child.new.opt[2]))
    
    axs <- c(1-lt.out$lt[1,2], 4-lt.out$lt[2,2])
    mx1 <- qx.child.new.opt[1]/(1-axs[1]*qx.child.new.opt[1])
    mx2 <- qx.child.new.opt[2]/(4-axs[2]*qx.child.new.opt[2])
    
    mx.child.new.opt <- c(mx1, mx2)
    out.mort.new.opt <- unname(c(log(mx.child.new.opt), out.mort.start[3:length(out.mort.start)]))
    if(out.mort.new.opt[3]>out.mort.new.opt[2]){
      # y <- out.mort.new.opt[1:2]
      # x <- c(1,2)
      # child.mod.lm <- lm(y~x)
      # new.5q5 <- predict.lm(child.mod.lm, newdata=data.frame(x = 3))
      new.5m5 <- approx(x=c(2,4),y=out.mort.new.opt[c(2,4)],xout=3)$y
      out.mort.new.opt[3] <- new.5m5
    }
    #lt.out.new.opt <- LifeTableMx(mx=exp(out.mort.new.opt), sex="Male")
    lt.out.new.opt <- lt.mx(nmx=exp(out.mort.new.opt), sex="male", age=ages, nax=NULL)
    out.mort <- out.mort.new.opt
  } # if
  
  
  return(exp(out.mort))
}

if(region==0 & sex==1){ # Non-Africa, female
  intercept <- median(svd.coeffs.xp.f[la.nums,1])
  b1.f <- predict.lm(co1.45q15mod.na.f, newdata=data.frame(cmort.na.f=child.mort, amort.na.f=adult.mort))
  b2.f <- predict.lm(co2.45q15mod.na.f, newdata=data.frame(cmort.na.f=child.mort, amort.na.f=adult.mort, prev.na=prev))
  b3.f <- predict.lm(co3.45q15mod.na.f, newdata=data.frame(cmort.na.f=child.mort, amort.na.f=adult.mort, prev.na=prev))
  
  if(!is.null(adult.art)){
    b2.f <- predict.lm(co2.45q15mod2.na.f, newdata=data.frame(cmort.na.f=child.mort, amort.na.f=adult.mort, prev.na=prev, aart.na=adult.art))
    b3.f <- predict.lm(co3.45q15mod2.na.f, newdata=data.frame(cmort.na.f=child.mort, amort.na.f=adult.mort, prev.na=prev, aart.na=adult.art)) 
  }
  
  if(!is.null(adult.art)&!is.null(child.art)){
    b2.f <- predict.lm(co2.45q15mod3.na.f, newdata=data.frame(cmort.na.f=child.mort, amort.na.f=adult.mort, prev.na=prev, aart.na=adult.art, cart.na=child.art))
    b3.f <- predict.lm(co3.45q15mod3.na.f, newdata=data.frame(cmort.na.f=child.mort, amort.na.f=adult.mort, prev.na=prev, aart.na=adult.art, cart.na=child.art)) 
  }
  
  if(opt==FALSE){
    out.mort <- intercept + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]} # if
  
  # optimize the intercept when predicting the mortality rates from the weights
  if(opt==TRUE){
    out.mort.func <- function(intercept.alter){
      out.mort <- intercept.alter + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]
      #lt.out <- LifeTableMx(mx=exp(out.mort), sex="Female")
      lt.out <- lt.mx(nmx=exp(out.mort), sex="female", age=ages)
      #lt.out.45q15 <- 1-lt.out$lx[14]/lt.out$lx[5]
      lt.out.45q15 <- lt.out$lt.45q15
      amort.diff <- abs(adult.mort-lt.out.45q15)     
      return(amort.diff)  
    } # function to be optimized
    
    intercept.opt <- optimize(f=out.mort.func, interval=c(-2,2))$minimum
    
    out.mort <- intercept.opt + b1.f*Mx.svd.scores.f[,1] + b2.f*Mx.svd.scores.f[,2] + b3.f*Mx.svd.scores.f[,3]} # if
  
  
  ## optimize reduction or addition to 1q0 and 4q1 to make it match
  if(opt==TRUE){
    out.mort.start <- out.mort
    #lt.out <- LifeTableMx(mx=exp(out.mort.start), sex="Female")
    lt.out <- lt.mx(nmx=exp(out.mort.start), sex="female", age=ages)
    
    diff.weight.opt <- log(1-child.mort)/(log(1-lt.out$lt[1,4])+log(1-lt.out$lt[2,4]))
    qx.child.new.opt <- 1-(c((1-lt.out$lt[1,4])^diff.weight.opt, (1-lt.out$lt[2,4])^diff.weight.opt))
    cmort.opt <- 1-((1-qx.child.new.opt[1])*(1-qx.child.new.opt[2]))
    
    axs <- c(1-lt.out$lt[1,2], 4-lt.out$lt[2,2])
    mx1 <- qx.child.new.opt[1]/(1-axs[1]*qx.child.new.opt[1])
    mx2 <- qx.child.new.opt[2]/(4-axs[2]*qx.child.new.opt[2])
    
    mx.child.new.opt <- c(mx1, mx2)
    out.mort.new.opt <- unname(c(log(mx.child.new.opt), out.mort.start[3:length(out.mort.start)]))
    if(out.mort.new.opt[3]>out.mort.new.opt[2]){
      # y <- out.mort.new.opt[1:2]
      # x <- c(1,2)
      # child.mod.lm <- lm(y~x)
      # new.5q5 <- predict.lm(child.mod.lm, newdata=data.frame(x = 3))
      new.5m5 <- approx(x=c(2,4),y=out.mort.new.opt[c(2,4)],xout=3)$y
      out.mort.new.opt[3] <- new.5m5
    }
    #lt.out.new.opt <- LifeTableMx(mx=exp(out.mort.new.opt), sex="Female")
    lt.out.new.opt <- lt.mx(nmx=exp(out.mort.new.opt), sex="female", age=ages, nax=NULL)
    out.mort <- out.mort.new.opt
  } # if
  
  
  return(exp(out.mort))
}
}
