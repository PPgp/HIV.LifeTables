hiv.mlt.recal <- function(mx.m, mx.f, prev, c.art, a.art, ages=c(0,1,seq(5,100,5)), la.countries=c("Bahamas", "Belize", "Cambodia", "Dominican Republic", "Estonia", "Guyana", "Haiti", "Honduras", "Jamaica", "Panama", "Russian Federation", "Suriname", "Thailand", "Trinidad and Tobago", "Ukraine"), midperiod=seq(1973,2008,5), determ=TRUE, user.par=NULL, nax=NULL, save.output="Ws-models.Rdata"){
  ###########################################################
  ## Set up objects to be used in function from input data ##
  ###########################################################
  
  ## Check if mx and epi data have same countries and periods -- otherwise error and stop ##
  #   if(nrow(mx.m)==nrow(mx.f)|ncol(mx.m)!=ncol(mx.f)){
  #     stop("mx input matrixes do not match columns or rows")
  #   }
  # 
  #   if(prev$country!=c.art$country|prev$country!=a.art$country|a.art$country!=c.art$country){
  #     stop("prevalence and ART input matrixes have discrepancy in country")
  #   }
  
  ## make a data frame with the prevalence and art info in one object
  print("Making data frame for epidemiological data...")
  midper.X <- paste("X",midperiod,sep="")
  cols.midper <- rep(NA,length(midper.X))
  for(i in 1:length(midper.X)){
    cols.midper[i] <- which(names(prev)==midper.X[i])
  }
  order.nums.rows <- order(rep(prev$country,length(midperiod)))
  prev.art.frame <- data.frame(sort(rep(prev$country,length(midperiod))), rep(prev$country_code,length(midperiod))[order.nums.rows], rep(midperiod,length(prev$country)), c(as.matrix(prev[,cols.midper]))[order.nums.rows], c(as.matrix(a.art[,cols.midper]))[order.nums.rows], c(as.matrix(c.art[,cols.midper]))[order.nums.rows])
  names(prev.art.frame) <- c("country", "country_code", "midperiod", "prev", "aart", "cart")
  
  Intervals <- length(ages)  
  
  # make matrix with mortality rates
  print("Making data frame for mortality data...")
  country.names.mx <- as.character(unique(mx.m$country))  
  for(i in order(country.names.mx)){
    if(i==order(country.names.mx)[1]){
      in.mx.m <- as.matrix(mx.m[mx.m$country==country.names.mx[i],-c(1:3)])
      in.mx.f <- as.matrix(mx.f[mx.f$country==country.names.mx[i],-c(1:3)])
    } else {
      in.mx.m <- cbind(in.mx.m, as.matrix(mx.m[mx.m$country==country.names.mx[i],-c(1:3)]))
      in.mx.f <- cbind(in.mx.f, as.matrix(mx.f[mx.f$country==country.names.mx[i],-c(1:3)]))
    } ## else
  } ## i loop
  
  #   # make sure prev.art.frame and mx have same number of cases
  #   if(ncol(in.mx.m)!=nrow(prev.art.frame)){
  #     stop("number of country-periods in mx input does not match number of country-periods in epi input")
  #   }
  
  # get input life tables to derive 5q0 and 45q15
  print("Calculating input life tables...")
  which.opt <- c(which(ages==user.par[1]),which(ages==user.par[2]))
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
  
  #library(bayesPop) ## may not need to call bayesPop package if we build a package out of this function
  in.lt.m <- array(NA, dim=c(Intervals, 14, ncol(in.mx.m)))
  in.lt.f <- array(NA, dim=c(Intervals, 14, ncol(in.mx.f)))
  for(i in 1:ncol(in.mx.m)){
    if(!is.null(nax)){
      lt.i.m <- lt.mx(nmx=unname(in.mx.m[,i]), sex="male", age=ages, nax=nax[i,])
      lt.i.f <- lt.mx(nmx=unname(in.mx.f[,i]), sex="female", age=ages, nax=nax[i,])      
    } else {
      lt.i.m <- lt.mx(nmx=unname(in.mx.m[,i]), sex="male", age=ages)
      lt.i.f <- lt.mx(nmx=unname(in.mx.f[,i]), sex="female", age=ages)   
    }
  
    in.lt.m[1:dim(lt.i.m$lt)[1],1:dim(lt.i.m$lt)[2],i] <- lt.i.m$lt
    in.lt.f[1:dim(lt.i.f$lt)[1],1:dim(lt.i.f$lt)[2],i] <- lt.i.f$lt
    in.lt.m[1,11,i] <- lt.i.m$e0
    in.lt.f[1,11,i] <- lt.i.f$e0
    in.lt.m[1,12,i] <- lt.i.m$lt.5q0
    in.lt.f[1,12,i] <- lt.i.f$lt.5q0
    in.lt.m[1,13,i] <- lt.i.m$lt.45q15
    in.lt.f[1,13,i] <- lt.i.f$lt.45q15
    if(!is.null(user.par)){
      in.lt.m[1,14,i] <- 1-(lt.i.m$lt[which(ages==user.par[2]),7]/lt.i.m$lt[which(ages==user.par[1]),7])
      in.lt.f[1,14,i] <- 1-(lt.i.f$lt[which(ages==user.par[2]),7]/lt.i.f$lt[which(ages==user.par[1]),7])  
    }
  } # i loop
  
  in.e0.m <- in.lt.m[1,11,]
  in.e0.f <- in.lt.f[1,11,]
  in.5q0.m <- in.lt.m[1,12,]
  in.5q0.f <- in.lt.f[1,12,]
  in.45q15.m <- in.lt.m[1,13,]
  in.45q15.f <- in.lt.f[1,13,]
  in.ep.m <- in.lt.m[1,14,]
  in.ep.f <- in.lt.f[1,14,]
  
  # which rows are African countries so we can make non-African and African models later
  for(i in 1:length(la.countries)){
    la.country.i <- unique(la.countries)[i]
    if(i==1){
      la.nums <- which(prev.art.frame$country==la.country.i)
    } else {
      la.nums <- c(la.nums, which(prev.art.frame$country==la.country.i))
    } # else
  } # i loop
  
  africa.nums <- c(1:nrow(prev.art.frame))[-la.nums]
  africa <- rep(1,nrow(prev.art.frame))
  africa[la.nums] <- 0
  
  #########
  ## SVD ##
  #########
  ## : perform Singular Value Decompositon on mortality rate schedules
  print("SVD...")
  Mx.xp.ln.m <- log(in.mx.m)
  Mx.xp.ln.f <- log(in.mx.f)

  svd.m <- svd(Mx.xp.ln.m) # SVD
  svd.f <- svd(Mx.xp.ln.f) # SVD

  Mx.svd.scores.m <- svd.m$u[,1:10]
  Mx.svd.scores.f <- svd.f$u[,1:10]
  
  #############
  ## Regress ##
  #############
  ## : regress the mortaity rate schdeules on the SVD components 
  print("Regress nmx on SVD...")
  num.components <- 10

  svd.coeffs.m <- lm(Mx.xp.ln.m ~ Mx.svd.scores.m[,1:num.components])$coefficients 
  svd.coeffs.f <- lm(Mx.xp.ln.f ~ Mx.svd.scores.f[,1:num.components])$coefficients 

  # transpose coefficients
  svd.coeffs.xp.m <- matrix(nrow=dim(svd.coeffs.m)[2],ncol=dim(svd.coeffs.m)[1])
  for(i in 1:dim(svd.coeffs.m)[2]){
    for(j in 1:dim(svd.coeffs.m)[1]){
      svd.coeffs.xp.m[i,j] <- svd.coeffs.m[j,i]
    }  # j loop
  } # i loop

  svd.coeffs.xp.f <- matrix(nrow=dim(svd.coeffs.f)[2],ncol=dim(svd.coeffs.f)[1])
  for(i in 1:dim(svd.coeffs.f)[2]){
    for(j in 1:dim(svd.coeffs.f)[1]){
      svd.coeffs.xp.f[i,j] <- svd.coeffs.f[j,i]
    }  # j loop
  } # i loop
  
  
  
  #################
  ## Make models ##
  #################
  ### used Bayesian Model Averaging for model selection -- code not included here #### 
  print("Calibrating weight models...")
  ## set up independent variables for modeling weights (svd.coeffs)
  in.prev <- prev.art.frame$prev
  in.aart <- prev.art.frame$aart
  in.cart <- prev.art.frame$cart
  
  # independent variables for African countries
  prev.a <- in.prev[africa.nums]
  aart.a <- in.aart[africa.nums]
  cart.a <- in.cart[africa.nums]
  le.a.m <- in.e0.m[africa.nums]
  le.a.f <- in.e0.f[africa.nums]
  cmort.a.m <- in.5q0.m[africa.nums]
  cmort.a.f <- in.5q0.f[africa.nums]
  amort.a.m <- in.45q15.m[africa.nums]
  amort.a.f <- in.45q15.f[africa.nums]
  
  # independent variables for non-African countries
  prev.na <- in.prev[la.nums]
  aart.na <- in.aart[la.nums]
  cart.na <- in.cart[la.nums]
  le.na.m <- in.e0.m[la.nums]
  le.na.f <- in.e0.f[la.nums]
  cmort.na.m <- in.5q0.m[la.nums]
  cmort.na.f <- in.5q0.f[la.nums]
  amort.na.m <- in.45q15.m[la.nums]
  amort.na.f <- in.45q15.f[la.nums]
  
  in.e0.average <- (in.e0.m+in.e0.f)/2
  le.total.a <- in.e0.average[africa.nums]
  le.total.na <- in.e0.average[la.nums]

  # additional entry parameters for African countries
  if(!is.null(user.par)){
    ep.a.m <- in.ep.m[africa.nums]
    ep.a.f <- in.ep.f[africa.nums]
    ep.na.m <- in.ep.m[la.nums]
    ep.na.f <- in.ep.f[la.nums]
  }


  # using three components is sufficient about 99.5% of age variation 
  ## e0 and prevalence
  ## sex-specific e0 models
  # # Africa: male
  co1.e0mod.a.m <- lm(svd.coeffs.xp.m[africa.nums,2] ~ le.a.m)
  co2.e0mod.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ le.a.m + prev.a)
  co3.e0mod.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ le.a.m + prev.a)
  
  # Africa: female
  co1.e0mod.a.f <- lm(svd.coeffs.xp.f[africa.nums,2] ~ le.a.f)
  co2.e0mod.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ le.a.f + prev.a)
  co3.e0mod.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ le.a.f + prev.a)
  
  # non-Africa: male
  co1.e0mod.na.m <- lm(svd.coeffs.xp.m[la.nums,2] ~ le.na.m)
  co2.e0mod.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ le.na.m + prev.na)
  co3.e0mod.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ le.na.m + prev.na)
  
  # non-Africa: female
  co1.e0mod.na.f <- lm(svd.coeffs.xp.f[la.nums,2] ~ le.na.f)
  co2.e0mod.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ le.na.f + prev.na)
  co3.e0mod.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ le.na.f + prev.na)
  
  ## e0 and prevalence and adult ART
  ## sex-specific e0 models
  # # Africa: male
  co2.e0mod2.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ le.a.m + prev.a + aart.a)
  co3.e0mod2.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ le.a.m + prev.a + aart.a)
  
  # Africa: female
  co2.e0mod2.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ le.a.f + prev.a + aart.a)
  co3.e0mod2.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ le.a.f + prev.a + aart.a)
  
  # non-Africa: male
  co2.e0mod2.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ le.na.m + prev.na + aart.na)
  co3.e0mod2.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ le.na.m + prev.na + aart.na)
  
  # non-Africa: female
  co2.e0mod2.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ le.na.f + prev.na + aart.na)
  co3.e0mod2.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ le.na.f + prev.na + aart.na)
  
  ## e0 and prevalence and adult ART and child ART
  ## sex-specific e0 models
  # # Africa: male
  co2.e0mod3.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ le.a.m + prev.a + aart.a + cart.a)
  co3.e0mod3.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ le.a.m + prev.a + aart.a + cart.a)
  
  # Africa: female
  co2.e0mod3.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ le.a.f + prev.a + aart.a + cart.a)
  co3.e0mod3.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ le.a.f + prev.a + aart.a + cart.a)
  
  # non-Africa: male
  co2.e0mod3.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ le.na.m + prev.na + aart.na + cart.na)
  co3.e0mod3.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ le.na.m + prev.na + aart.na + cart.na)
  
  # non-Africa: female
  co2.e0mod3.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ le.na.f + prev.na + aart.na + cart.na)
  co3.e0mod3.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ le.na.f + prev.na + aart.na + cart.na)
  
  
  ## 5q0 and prevalence 
  # Africa: male
  co1.5q0mod.a.m <- lm(svd.coeffs.xp.m[africa.nums,2] ~ cmort.a.m)
  co2.5q0mod.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ cmort.a.m + prev.a)
  co3.5q0mod.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ cmort.a.m + prev.a)
  
  # Africa: female
  co1.5q0mod.a.f <- lm(svd.coeffs.xp.f[africa.nums,2] ~ cmort.a.f)
  co2.5q0mod.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ cmort.a.f + prev.a)
  co3.5q0mod.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ cmort.a.f + prev.a)
  
  # non-Africa: male
  co1.5q0mod.na.m <- lm(svd.coeffs.xp.m[la.nums,2] ~ cmort.na.m)
  co2.5q0mod.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ cmort.na.m + prev.na)
  co3.5q0mod.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ cmort.na.m + prev.na)
  
  # non-Africa: female
  co1.5q0mod.na.f <- lm(svd.coeffs.xp.f[la.nums,2] ~ cmort.na.f)
  co2.5q0mod.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ cmort.na.f + prev.na)
  co3.5q0mod.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ cmort.na.f + prev.na)
  
  ## 5q0 and prevalence and adult ART
  # Africa: male
  co2.5q0mod2.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ cmort.a.m + prev.a + aart.a)
  co3.5q0mod2.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ cmort.a.m + prev.a + aart.a)
  
  # Africa: female
  co2.5q0mod2.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ cmort.a.f + prev.a + aart.a)
  co3.5q0mod2.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ cmort.a.f + prev.a + aart.a)
  
  # non-Africa: male
  co2.5q0mod2.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ cmort.na.m + prev.na + aart.na)
  co3.5q0mod2.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ cmort.na.m + prev.na + aart.na)
  
  # non-Africa: female
  co2.5q0mod2.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ cmort.na.f + prev.na + aart.na)
  co3.5q0mod2.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ cmort.na.f + prev.na + aart.na)
  
  ## 5q0 and prevalence and adult ART and child ART
  # Africa: male
  co2.5q0mod3.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ cmort.a.m + prev.a + aart.a + cart.a)
  co3.5q0mod3.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ cmort.a.m + prev.a + aart.a + cart.a)
  
  # Africa: female
  co2.5q0mod3.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ cmort.a.f + prev.a + aart.a + cart.a)
  co3.5q0mod3.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ cmort.a.f + prev.a + aart.a + cart.a)
  
  # non-Africa: male
  co2.5q0mod3.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ cmort.na.m + prev.na + aart.na + cart.na)
  co3.5q0mod3.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ cmort.na.m + prev.na + aart.na + cart.na)
  
  # non-Africa: female
  co2.5q0mod3.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ cmort.na.f + prev.na + aart.na + cart.na)
  co3.5q0mod3.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ cmort.na.f + prev.na + aart.na + cart.na)
  
  ## 5q0 and 45q15 and prevalence
  # when modeling the first weight for African countries just adding 45q15 to the 5q0 model causes the sign on prevalence to change in the wrong direction so prev is dropped from the regressions for the first weight
  # # Africa: male
  co1.45q15mod.a.m <- lm(svd.coeffs.xp.m[africa.nums,2] ~ cmort.a.m + amort.a.m)
  co2.45q15mod.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ cmort.a.m + amort.a.m + prev.a)
  co3.45q15mod.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ cmort.a.m + amort.a.m + prev.a)
  
  # Africa: female
  co1.45q15mod.a.f <- lm(svd.coeffs.xp.f[africa.nums,2] ~ cmort.a.f + amort.a.f)
  co2.45q15mod.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ cmort.a.f + amort.a.f + prev.a)
  co3.45q15mod.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ cmort.a.f + amort.a.f + prev.a)
  
  # non-Africa: male 
  co1.45q15mod.na.m <- lm(svd.coeffs.xp.m[la.nums,2] ~ cmort.na.m + amort.na.m)
  co2.45q15mod.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ cmort.na.m + amort.na.m + prev.na)
  co3.45q15mod.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ cmort.na.m + amort.na.m + prev.na)
  
  # non-Africa: female
  co1.45q15mod.na.f <- lm(svd.coeffs.xp.f[la.nums,2] ~ cmort.na.f + amort.na.f)
  co2.45q15mod.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ cmort.na.f + amort.na.f + prev.na)
  co3.45q15mod.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ cmort.na.f + amort.na.f + prev.na)
  
  ## 5q0 and 45q15 and prevalence and adult ART
  # when modeling the first weight for African countries just adding 45q15 to the 5q0 model causes the sign on prevalence to change in the wrong direction so prev is dropped from the regressions for the first weight
  # # Africa: male
  co2.45q15mod2.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ cmort.a.m + amort.a.m + prev.a + aart.a)
  co3.45q15mod2.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ cmort.a.m + amort.a.m + prev.a + aart.a)
  
  # Africa: female
  co2.45q15mod2.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ cmort.a.f + amort.a.f + prev.a + aart.a)
  co3.45q15mod2.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ cmort.a.f + amort.a.f + prev.a + aart.a)
  
  # non-Africa: male 
  co2.45q15mod2.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ cmort.na.m + amort.na.m + prev.na + aart.na)
  co3.45q15mod2.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ cmort.na.m + amort.na.m + prev.na + aart.na)
  
  # non-Africa: female
  co2.45q15mod2.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ cmort.na.f + amort.na.f + prev.na + aart.na)
  co3.45q15mod2.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ cmort.na.f + amort.na.f + prev.na + aart.na)
  
  ## 5q0 and 45q15 and prevalence and adult ART and child ART
  # when modeling the first weight for African countries just adding 45q15 to the 5q0 model causes the sign on prevalence to change in the wrong direction so prev is dropped from the regressions for the first weight
  # # Africa: male
  co2.45q15mod3.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ cmort.a.m + amort.a.m + prev.a + aart.a + cart.a)
  co3.45q15mod3.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ cmort.a.m + amort.a.m + prev.a + aart.a + cart.a)
  
  # Africa: female
  co2.45q15mod3.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ cmort.a.f + amort.a.f + prev.a + aart.a + cart.a)
  co3.45q15mod3.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ cmort.a.f + amort.a.f + prev.a + aart.a + cart.a)
  
  # non-Africa: male 
  co2.45q15mod3.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ cmort.na.m + amort.na.m + prev.na + aart.na + cart.na)
  co3.45q15mod3.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ cmort.na.m + amort.na.m + prev.na + aart.na + cart.na)
  
  # non-Africa: female
  co2.45q15mod3.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ cmort.na.f + amort.na.f + prev.na + aart.na + cart.na)
  co3.45q15mod3.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ cmort.na.f + amort.na.f + prev.na + aart.na + cart.na)
  
## user-supplied entry parameter and prevalence
if(!is.null(user.par)){
# # Africa: male
co1.epmod.a.m <- lm(svd.coeffs.xp.m[africa.nums,2] ~ ep.a.m)
co2.epmod.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ ep.a.m + prev.a)
co3.epmod.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ ep.a.m + prev.a)

# Africa: female
co1.epmod.a.f <- lm(svd.coeffs.xp.f[africa.nums,2] ~ ep.a.f)
co2.epmod.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ ep.a.f + prev.a)
co3.epmod.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ ep.a.f + prev.a)

# non-Africa: male
co1.epmod.na.m <- lm(svd.coeffs.xp.m[la.nums,2] ~ ep.na.m)
co2.epmod.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ ep.na.m + prev.na)
co3.epmod.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ ep.na.m + prev.na)

# non-Africa: female
co1.epmod.na.f <- lm(svd.coeffs.xp.f[la.nums,2] ~ ep.na.f)
co2.epmod.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ ep.na.f + prev.na)
co3.epmod.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ ep.na.f + prev.na)

## user-supplied entry parameter and prevalence and adult ART
# # Africa: male
co2.epmod2.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ ep.a.m + prev.a + aart.a)
co3.epmod2.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ ep.a.m + prev.a + aart.a)

# Africa: female
co2.epmod2.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ ep.a.f + prev.a + aart.a)
co3.epmod2.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ ep.a.f + prev.a + aart.a)

# non-Africa: male
co2.epmod2.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ ep.na.m + prev.na + aart.na)
co3.epmod2.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ ep.na.m + prev.na + aart.na)

# non-Africa: female
co2.epmod2.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ ep.na.f + prev.na + aart.na)
co3.epmod2.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ ep.na.f + prev.na + aart.na)

## user-supplied entry parameter and prevalence and adult ART and child ART
# # Africa: male
co2.epmod3.a.m <- lm(svd.coeffs.xp.m[africa.nums,3] ~ ep.a.m + prev.a + aart.a + cart.a)
co3.epmod3.a.m <- lm(svd.coeffs.xp.m[africa.nums,4] ~ ep.a.m + prev.a + aart.a + cart.a)

# Africa: female
co2.epmod3.a.f <- lm(svd.coeffs.xp.f[africa.nums,3] ~ ep.a.f + prev.a + aart.a + cart.a)
co3.epmod3.a.f <- lm(svd.coeffs.xp.f[africa.nums,4] ~ ep.a.f + prev.a + aart.a + cart.a)

# non-Africa: male
co2.epmod3.na.m <- lm(svd.coeffs.xp.m[la.nums,3] ~ ep.na.m + prev.na + aart.na + cart.na)
co3.epmod3.na.m <- lm(svd.coeffs.xp.m[la.nums,4] ~ ep.na.m + prev.na + aart.na + cart.na)

# non-Africa: female
co2.epmod3.na.f <- lm(svd.coeffs.xp.f[la.nums,3] ~ ep.na.f + prev.na + aart.na + cart.na)
co3.epmod3.na.f <- lm(svd.coeffs.xp.f[la.nums,4] ~ ep.na.f + prev.na + aart.na + cart.na)
}

  # put the results in arrays for displaying function output quickly
if(is.null(user.par)){
  # Male, African
  coeff.array.a.m <- array(NA, dim=c(7,9,3), dimnames=list(c("Intercept", "e0", "5q0", "45q15", "prev", "ART (adult)", "ART (child)"),c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9"),c("co1", "co2", "co3")))
  coeff.array.a.m[1:2,1,1] <- co1.e0mod.a.m$coeff
  coeff.array.a.m[c(1,3),4,1] <- co1.5q0mod.a.m$coeff
  coeff.array.a.m[c(1,3:4),7,1] <- co1.45q15mod.a.m$coeff

  coeff.array.a.m[c(1:2,5),1,2]   <- co2.e0mod.a.m$coeff
  coeff.array.a.m[c(1:2,5:6),2,2] <- co2.e0mod2.a.m$coeff
  coeff.array.a.m[c(1:2,5:7),3,2] <- co2.e0mod3.a.m$coeff
  coeff.array.a.m[c(1,3,5),4,2]   <- co2.5q0mod.a.m$coeff
  coeff.array.a.m[c(1,3,5:6),5,2] <- co2.5q0mod2.a.m$coeff
  coeff.array.a.m[c(1,3,5:7),6,2] <- co2.5q0mod3.a.m$coeff
  coeff.array.a.m[c(1,3:5),7,2]   <- co2.45q15mod.a.m$coeff
  coeff.array.a.m[c(1,3:6),8,2]   <- co2.45q15mod2.a.m$coeff
  coeff.array.a.m[c(1,3:7),9,2]   <- co2.45q15mod3.a.m$coeff
  
  coeff.array.a.m[c(1:2,5),1,3]   <- co3.e0mod.a.m$coeff
  coeff.array.a.m[c(1:2,5:6),2,3] <- co3.e0mod2.a.m$coeff
  coeff.array.a.m[c(1:2,5:7),3,3] <- co3.e0mod3.a.m$coeff
  coeff.array.a.m[c(1,3,5),4,3]   <- co3.5q0mod.a.m$coeff
  coeff.array.a.m[c(1,3,5:6),5,3] <- co3.5q0mod2.a.m$coeff
  coeff.array.a.m[c(1,3,5:7),6,3] <- co3.5q0mod3.a.m$coeff
  coeff.array.a.m[c(1,3:5),7,3]   <- co3.45q15mod.a.m$coeff
  coeff.array.a.m[c(1,3:6),8,3]   <- co3.45q15mod2.a.m$coeff
  coeff.array.a.m[c(1,3:7),9,3]   <- co3.45q15mod3.a.m$coeff
  
  # Male, non-African
  coeff.array.na.m <- array(NA, dim=c(7,9,3), dimnames=list(c("Intercept", "e0", "5q0", "45q15", "prev", "ART (adult)", "ART (child)"),c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9"),c("co1", "co2", "co3")))
  coeff.array.na.m[1:2,1,1] <- co1.e0mod.na.m$coeff
  coeff.array.na.m[c(1,3),4,1] <- co1.5q0mod.na.m$coeff
  coeff.array.na.m[c(1,3:4),7,1] <- co1.45q15mod.na.m$coeff
  
  coeff.array.na.m[c(1:2,5),1,2]   <- co2.e0mod.na.m$coeff
  coeff.array.na.m[c(1:2,5:6),2,2] <- co2.e0mod2.na.m$coeff
  coeff.array.na.m[c(1:2,5:7),3,2] <- co2.e0mod3.na.m$coeff
  coeff.array.na.m[c(1,3,5),4,2]   <- co2.5q0mod.na.m$coeff
  coeff.array.na.m[c(1,3,5:6),5,2] <- co2.5q0mod2.na.m$coeff
  coeff.array.na.m[c(1,3,5:7),6,2] <- co2.5q0mod3.na.m$coeff
  coeff.array.na.m[c(1,3:5),7,2]   <- co2.45q15mod.na.m$coeff
  coeff.array.na.m[c(1,3:6),8,2]   <- co2.45q15mod2.na.m$coeff
  coeff.array.na.m[c(1,3:7),9,2]   <- co2.45q15mod3.na.m$coeff
  
  coeff.array.na.m[c(1:2,5),1,3]   <- co3.e0mod.na.m$coeff
  coeff.array.na.m[c(1:2,5:6),2,3] <- co3.e0mod2.na.m$coeff
  coeff.array.na.m[c(1:2,5:7),3,3] <- co3.e0mod3.na.m$coeff
  coeff.array.na.m[c(1,3,5),4,3]   <- co3.5q0mod.na.m$coeff
  coeff.array.na.m[c(1,3,5:6),5,3] <- co3.5q0mod2.na.m$coeff
  coeff.array.na.m[c(1,3,5:7),6,3] <- co3.5q0mod3.na.m$coeff
  coeff.array.na.m[c(1,3:5),7,3]   <- co3.45q15mod.na.m$coeff
  coeff.array.na.m[c(1,3:6),8,3]   <- co3.45q15mod2.na.m$coeff
  coeff.array.na.m[c(1,3:7),9,3]   <- co3.45q15mod3.na.m$coeff
  
  # Female, African
  coeff.array.a.f <- array(NA, dim=c(7,9,3), dimnames=list(c("Intercept", "e0", "5q0", "45q15", "prev", "ART (adult)", "ART (child)"),c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9"),c("co1", "co2", "co3")))
  coeff.array.a.f[1:2,1,1] <- co1.e0mod.a.f$coeff
  coeff.array.a.f[c(1,3),4,1] <- co1.5q0mod.a.f$coeff
  coeff.array.a.f[c(1,3:4),7,1] <- co1.45q15mod.a.f$coeff
  
  coeff.array.a.f[c(1:2,5),1,2]   <- co2.e0mod.a.f$coeff
  coeff.array.a.f[c(1:2,5:6),2,2] <- co2.e0mod2.a.f$coeff
  coeff.array.a.f[c(1:2,5:7),3,2] <- co2.e0mod3.a.f$coeff
  coeff.array.a.f[c(1,3,5),4,2]   <- co2.5q0mod.a.f$coeff
  coeff.array.a.f[c(1,3,5:6),5,2] <- co2.5q0mod2.a.f$coeff
  coeff.array.a.f[c(1,3,5:7),6,2] <- co2.5q0mod3.a.f$coeff
  coeff.array.a.f[c(1,3:5),7,2]   <- co2.45q15mod.a.f$coeff
  coeff.array.a.f[c(1,3:6),8,2]   <- co2.45q15mod2.a.f$coeff
  coeff.array.a.f[c(1,3:7),9,2]   <- co2.45q15mod3.a.f$coeff
  
  coeff.array.a.f[c(1:2,5),1,3]   <- co3.e0mod.a.f$coeff
  coeff.array.a.f[c(1:2,5:6),2,3] <- co3.e0mod2.a.f$coeff
  coeff.array.a.f[c(1:2,5:7),3,3] <- co3.e0mod3.a.f$coeff
  coeff.array.a.f[c(1,3,5),4,3]   <- co3.5q0mod.a.f$coeff
  coeff.array.a.f[c(1,3,5:6),5,3] <- co3.5q0mod2.a.f$coeff
  coeff.array.a.f[c(1,3,5:7),6,3] <- co3.5q0mod3.a.f$coeff
  coeff.array.a.f[c(1,3:5),7,3]   <- co3.45q15mod.a.f$coeff
  coeff.array.a.f[c(1,3:6),8,3]   <- co3.45q15mod2.a.f$coeff
  coeff.array.a.f[c(1,3:7),9,3]   <- co3.45q15mod3.a.f$coeff
  
  # Female, non-African
  coeff.array.na.f <- array(NA, dim=c(7,9,3), dimnames=list(c("Intercept", "e0", "5q0", "45q15", "prev", "ART (adult)", "ART (child)"),c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9"),c("co1", "co2", "co3")))
  coeff.array.na.f[1:2,1,1] <- co1.e0mod.na.f$coeff
  coeff.array.na.f[c(1,3),4,1] <- co1.5q0mod.na.f$coeff
  coeff.array.na.f[c(1,3:4),7,1] <- co1.45q15mod.na.f$coeff
  
  coeff.array.na.f[c(1:2,5),1,2]   <- co2.e0mod.na.f$coeff
  coeff.array.na.f[c(1:2,5:6),2,2] <- co2.e0mod2.na.f$coeff
  coeff.array.na.f[c(1:2,5:7),3,2] <- co2.e0mod3.na.f$coeff
  coeff.array.na.f[c(1,3,5),4,2]   <- co2.5q0mod.na.f$coeff
  coeff.array.na.f[c(1,3,5:6),5,2] <- co2.5q0mod2.na.f$coeff
  coeff.array.na.f[c(1,3,5:7),6,2] <- co2.5q0mod3.na.f$coeff
  coeff.array.na.f[c(1,3:5),7,2]   <- co2.45q15mod.na.f$coeff
  coeff.array.na.f[c(1,3:6),8,2]   <- co2.45q15mod2.na.f$coeff
  coeff.array.na.f[c(1,3:7),9,2]   <- co2.45q15mod3.na.f$coeff
  
  coeff.array.na.f[c(1:2,5),1,3]   <- co3.e0mod.na.f$coeff
  coeff.array.na.f[c(1:2,5:6),2,3] <- co3.e0mod2.na.f$coeff
  coeff.array.na.f[c(1:2,5:7),3,3] <- co3.e0mod3.na.f$coeff
  coeff.array.na.f[c(1,3,5),4,3]   <- co3.5q0mod.na.f$coeff
  coeff.array.na.f[c(1,3,5:6),5,3] <- co3.5q0mod2.na.f$coeff
  coeff.array.na.f[c(1,3,5:7),6,3] <- co3.5q0mod3.na.f$coeff
  coeff.array.na.f[c(1,3:5),7,3]   <- co3.45q15mod.na.f$coeff
  coeff.array.na.f[c(1,3:6),8,3]   <- co3.45q15mod2.na.f$coeff
  coeff.array.na.f[c(1,3:7),9,3]   <- co3.45q15mod3.na.f$coeff
  } else {
  # Male, African
  coeff.array.a.m <- array(NA, dim=c(8,12,3), dimnames=list(c("Intercept", "e0", "5q0", "45q15", "prev", "ART (adult)", "ART (child)", "user entry par"),c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12"),c("co1", "co2", "co3")))
  coeff.array.a.m[1:2,1,1] <- co1.e0mod.a.m$coeff
  coeff.array.a.m[c(1,3),4,1] <- co1.5q0mod.a.m$coeff
  coeff.array.a.m[c(1,3:4),7,1] <- co1.45q15mod.a.m$coeff
  coeff.array.a.m[c(1,8),10,1] <- co1.epmod.a.m$coeff
  
  coeff.array.a.m[c(1:2,5),1,2]   <- co2.e0mod.a.m$coeff
  coeff.array.a.m[c(1:2,5:6),2,2] <- co2.e0mod2.a.m$coeff
  coeff.array.a.m[c(1:2,5:7),3,2] <- co2.e0mod3.a.m$coeff
  coeff.array.a.m[c(1,3,5),4,2]   <- co2.5q0mod.a.m$coeff
  coeff.array.a.m[c(1,3,5:6),5,2] <- co2.5q0mod2.a.m$coeff
  coeff.array.a.m[c(1,3,5:7),6,2] <- co2.5q0mod3.a.m$coeff
  coeff.array.a.m[c(1,3:5),7,2]   <- co2.45q15mod.a.m$coeff
  coeff.array.a.m[c(1,3:6),8,2]   <- co2.45q15mod2.a.m$coeff
  coeff.array.a.m[c(1,3:7),9,2]   <- co2.45q15mod3.a.m$coeff
  coeff.array.a.m[c(1,8,5),10,2]    <- co2.epmod.a.m$coeff
  coeff.array.a.m[c(1,8,5:6),11,2]    <- co2.epmod2.a.m$coeff
  coeff.array.a.m[c(1,8,5:7),12,2]    <- co2.epmod3.a.m$coeff
  
  coeff.array.a.m[c(1:2,5),1,3]   <- co3.e0mod.a.m$coeff
  coeff.array.a.m[c(1:2,5:6),2,3] <- co3.e0mod2.a.m$coeff
  coeff.array.a.m[c(1:2,5:7),3,3] <- co3.e0mod3.a.m$coeff
  coeff.array.a.m[c(1,3,5),4,3]   <- co3.5q0mod.a.m$coeff
  coeff.array.a.m[c(1,3,5:6),5,3] <- co3.5q0mod2.a.m$coeff
  coeff.array.a.m[c(1,3,5:7),6,3] <- co3.5q0mod3.a.m$coeff
  coeff.array.a.m[c(1,3:5),7,3]   <- co3.45q15mod.a.m$coeff
  coeff.array.a.m[c(1,3:6),8,3]   <- co3.45q15mod2.a.m$coeff
  coeff.array.a.m[c(1,3:7),9,3]   <- co3.45q15mod3.a.m$coeff
  coeff.array.a.m[c(1,8,5),10,3]    <- co3.epmod.a.m$coeff
  coeff.array.a.m[c(1,8,5:6),11,3]    <- co3.epmod2.a.m$coeff
  coeff.array.a.m[c(1,8,5:7),12,3]    <- co3.epmod3.a.m$coeff
  
  # Male, non-African
  coeff.array.na.m <- array(NA, dim=c(8,12,3), dimnames=list(c("Intercept", "e0", "5q0", "45q15", "prev", "ART (adult)", "ART (child)", "user entry par"),c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12"),c("co1", "co2", "co3")))
  coeff.array.na.m[1:2,1,1] <- co1.e0mod.na.m$coeff
  coeff.array.na.m[c(1,3),4,1] <- co1.5q0mod.na.m$coeff
  coeff.array.na.m[c(1,3:4),7,1] <- co1.45q15mod.na.m$coeff
  coeff.array.na.m[c(1,8),10,1] <- co1.epmod.na.m$coeff
  
  coeff.array.na.m[c(1:2,5),1,2]   <- co2.e0mod.na.m$coeff
  coeff.array.na.m[c(1:2,5:6),2,2] <- co2.e0mod2.na.m$coeff
  coeff.array.na.m[c(1:2,5:7),3,2] <- co2.e0mod3.na.m$coeff
  coeff.array.na.m[c(1,3,5),4,2]   <- co2.5q0mod.na.m$coeff
  coeff.array.na.m[c(1,3,5:6),5,2] <- co2.5q0mod2.na.m$coeff
  coeff.array.na.m[c(1,3,5:7),6,2] <- co2.5q0mod3.na.m$coeff
  coeff.array.na.m[c(1,3:5),7,2]   <- co2.45q15mod.na.m$coeff
  coeff.array.na.m[c(1,3:6),8,2]   <- co2.45q15mod2.na.m$coeff
  coeff.array.na.m[c(1,3:7),9,2]   <- co2.45q15mod3.na.m$coeff
  coeff.array.na.m[c(1,8,5),10,2]    <- co2.epmod.na.m$coeff
  coeff.array.na.m[c(1,8,5:6),11,2]    <- co2.epmod2.na.m$coeff
  coeff.array.na.m[c(1,8,5:7),12,2]    <- co2.epmod3.na.m$coeff
  
  coeff.array.na.m[c(1:2,5),1,3]   <- co3.e0mod.na.m$coeff
  coeff.array.na.m[c(1:2,5:6),2,3] <- co3.e0mod2.na.m$coeff
  coeff.array.na.m[c(1:2,5:7),3,3] <- co3.e0mod3.na.m$coeff
  coeff.array.na.m[c(1,3,5),4,3]   <- co3.5q0mod.na.m$coeff
  coeff.array.na.m[c(1,3,5:6),5,3] <- co3.5q0mod2.na.m$coeff
  coeff.array.na.m[c(1,3,5:7),6,3] <- co3.5q0mod3.na.m$coeff
  coeff.array.na.m[c(1,3:5),7,3]   <- co3.45q15mod.na.m$coeff
  coeff.array.na.m[c(1,3:6),8,3]   <- co3.45q15mod2.na.m$coeff
  coeff.array.na.m[c(1,3:7),9,3]   <- co3.45q15mod3.na.m$coeff
  coeff.array.na.m[c(1,8,5),10,3]    <- co3.epmod.na.m$coeff
  coeff.array.na.m[c(1,8,5:6),11,3]    <- co3.epmod2.na.m$coeff
  coeff.array.na.m[c(1,8,5:7),12,3]    <- co3.epmod3.na.m$coeff
  
  # Female, African
  coeff.array.a.f <- array(NA, dim=c(8,12,3), dimnames=list(c("Intercept", "e0", "5q0", "45q15", "prev", "ART (adult)", "ART (child)", "user entry par"),c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12"),c("co1", "co2", "co3")))
  coeff.array.a.f[1:2,1,1] <- co1.e0mod.a.f$coeff
  coeff.array.a.f[c(1,3),4,1] <- co1.5q0mod.a.f$coeff
  coeff.array.a.f[c(1,3:4),7,1] <- co1.45q15mod.a.f$coeff
  coeff.array.a.f[c(1,8),10,1] <- co1.epmod.a.f$coeff
  
  coeff.array.a.f[c(1:2,5),1,2]   <- co2.e0mod.a.f$coeff
  coeff.array.a.f[c(1:2,5:6),2,2] <- co2.e0mod2.a.f$coeff
  coeff.array.a.f[c(1:2,5:7),3,2] <- co2.e0mod3.a.f$coeff
  coeff.array.a.f[c(1,3,5),4,2]   <- co2.5q0mod.a.f$coeff
  coeff.array.a.f[c(1,3,5:6),5,2] <- co2.5q0mod2.a.f$coeff
  coeff.array.a.f[c(1,3,5:7),6,2] <- co2.5q0mod3.a.f$coeff
  coeff.array.a.f[c(1,3:5),7,2]   <- co2.45q15mod.a.f$coeff
  coeff.array.a.f[c(1,3:6),8,2]   <- co2.45q15mod2.a.f$coeff
  coeff.array.a.f[c(1,3:7),9,2]   <- co2.45q15mod3.a.f$coeff
  coeff.array.a.f[c(1,8,5),10,2]    <- co2.epmod.a.f$coeff
  coeff.array.a.f[c(1,8,5:6),11,2]    <- co2.epmod2.a.f$coeff
  coeff.array.a.f[c(1,8,5:7),12,2]    <- co2.epmod3.a.f$coeff
  
  coeff.array.a.f[c(1:2,5),1,3]   <- co3.e0mod.a.f$coeff
  coeff.array.a.f[c(1:2,5:6),2,3] <- co3.e0mod2.a.f$coeff
  coeff.array.a.f[c(1:2,5:7),3,3] <- co3.e0mod3.a.f$coeff
  coeff.array.a.f[c(1,3,5),4,3]   <- co3.5q0mod.a.f$coeff
  coeff.array.a.f[c(1,3,5:6),5,3] <- co3.5q0mod2.a.f$coeff
  coeff.array.a.f[c(1,3,5:7),6,3] <- co3.5q0mod3.a.f$coeff
  coeff.array.a.f[c(1,3:5),7,3]   <- co3.45q15mod.a.f$coeff
  coeff.array.a.f[c(1,3:6),8,3]   <- co3.45q15mod2.a.f$coeff
  coeff.array.a.f[c(1,3:7),9,3]   <- co3.45q15mod3.a.f$coeff
  coeff.array.a.f[c(1,8,5),10,3]    <- co3.epmod.a.f$coeff
  coeff.array.a.f[c(1,8,5:6),11,3]    <- co3.epmod2.a.f$coeff
  coeff.array.a.f[c(1,8,5:7),12,3]    <- co3.epmod3.a.f$coeff
  
  # Female, non-African
  coeff.array.na.f <- array(NA, dim=c(8,12,3), dimnames=list(c("Intercept", "e0", "5q0", "45q15", "prev", "ART (adult)", "ART (child)", "user entry par"),c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7", "mod8", "mod9", "mod10", "mod11", "mod12"),c("co1", "co2", "co3")))
  coeff.array.na.f[1:2,1,1] <- co1.e0mod.na.f$coeff
  coeff.array.na.f[c(1,3),4,1] <- co1.5q0mod.na.f$coeff
  coeff.array.na.f[c(1,3:4),7,1] <- co1.45q15mod.na.f$coeff
  coeff.array.na.f[c(1,8),10,1] <- co1.epmod.na.f$coeff
  
  coeff.array.na.f[c(1:2,5),1,2]   <- co2.e0mod.na.f$coeff
  coeff.array.na.f[c(1:2,5:6),2,2] <- co2.e0mod2.na.f$coeff
  coeff.array.na.f[c(1:2,5:7),3,2] <- co2.e0mod3.na.f$coeff
  coeff.array.na.f[c(1,3,5),4,2]   <- co2.5q0mod.na.f$coeff
  coeff.array.na.f[c(1,3,5:6),5,2] <- co2.5q0mod2.na.f$coeff
  coeff.array.na.f[c(1,3,5:7),6,2] <- co2.5q0mod3.na.f$coeff
  coeff.array.na.f[c(1,3:5),7,2]   <- co2.45q15mod.na.f$coeff
  coeff.array.na.f[c(1,3:6),8,2]   <- co2.45q15mod2.na.f$coeff
  coeff.array.na.f[c(1,3:7),9,2]   <- co2.45q15mod3.na.f$coeff
  coeff.array.na.f[c(1,8,5),10,2]    <- co2.epmod.na.f$coeff
  coeff.array.na.f[c(1,8,5:6),11,2]    <- co2.epmod2.na.f$coeff
  coeff.array.na.f[c(1,8,5:7),12,2]    <- co2.epmod3.na.f$coeff
  
  coeff.array.na.f[c(1:2,5),1,3]   <- co3.e0mod.na.f$coeff
  coeff.array.na.f[c(1:2,5:6),2,3] <- co3.e0mod2.na.f$coeff
  coeff.array.na.f[c(1:2,5:7),3,3] <- co3.e0mod3.na.f$coeff
  coeff.array.na.f[c(1,3,5),4,3]   <- co3.5q0mod.na.f$coeff
  coeff.array.na.f[c(1,3,5:6),5,3] <- co3.5q0mod2.na.f$coeff
  coeff.array.na.f[c(1,3,5:7),6,3] <- co3.5q0mod3.na.f$coeff
  coeff.array.na.f[c(1,3:5),7,3]   <- co3.45q15mod.na.f$coeff
  coeff.array.na.f[c(1,3:6),8,3]   <- co3.45q15mod2.na.f$coeff
  coeff.array.na.f[c(1,3:7),9,3]   <- co3.45q15mod3.na.f$coeff
  coeff.array.na.f[c(1,8,5),10,3]    <- co3.epmod.na.f$coeff
  coeff.array.na.f[c(1,8,5:6),11,3]    <- co3.epmod2.na.f$coeff
  coeff.array.na.f[c(1,8,5:7),12,3]    <- co3.epmod3.na.f$coeff
}

if(is.null(user.par)){
  save.list <- c("co1.e0mod.a.m", "co2.e0mod.a.m", "co3.e0mod.a.m", 
                  "co2.e0mod2.a.m", "co3.e0mod2.a.m",
                  "co2.e0mod3.a.m", "co3.e0mod3.a.m",
                  "co1.e0mod.na.m", "co2.e0mod.na.m", "co3.e0mod.na.m", 
                  "co2.e0mod2.na.m", "co3.e0mod2.na.m",
                  "co2.e0mod3.na.m", "co3.e0mod3.na.m",
                  "co1.e0mod.a.f", "co2.e0mod.a.f", "co3.e0mod.a.f",
                  "co2.e0mod2.a.f", "co3.e0mod2.a.f",
                  "co2.e0mod3.a.f", "co3.e0mod3.a.f",
                  "co1.e0mod.na.f", "co2.e0mod.na.f", "co3.e0mod.na.f",
                  "co2.e0mod2.na.f", "co3.e0mod2.na.f",
                  "co2.e0mod3.na.f", "co3.e0mod3.na.f",
                  "co1.5q0mod.a.m", "co2.5q0mod.a.m", "co3.5q0mod.a.m",
                  "co2.5q0mod2.a.m", "co3.5q0mod2.a.m",
                  "co2.5q0mod3.a.m", "co3.5q0mod3.a.m",
                  "co1.5q0mod.na.m", "co2.5q0mod.na.m", "co3.5q0mod.na.m", 
                  "co2.5q0mod2.na.m", "co3.5q0mod2.na.m",
                  "co2.5q0mod3.na.m", "co3.5q0mod3.na.m",
                  "co1.5q0mod.a.f", "co2.5q0mod.a.f", "co3.5q0mod.a.f",
                  "co2.5q0mod2.a.f", "co3.5q0mod2.a.f",
                  "co2.5q0mod3.a.f", "co3.5q0mod3.a.f",
                  "co1.5q0mod.na.f", "co2.5q0mod.na.f", "co3.5q0mod.na.f",
                  "co2.5q0mod2.na.f", "co3.5q0mod2.na.f",
                  "co2.5q0mod3.na.f", "co3.5q0mod3.na.f",
                  "co1.45q15mod.a.m", "co2.45q15mod.a.m", "co3.45q15mod.a.m",
                  "co2.45q15mod2.a.m", "co3.45q15mod2.a.m",
                  "co2.45q15mod3.a.m", "co3.45q15mod3.a.m",
                  "co1.45q15mod.na.m", "co2.45q15mod.na.m", "co3.45q15mod.na.m",
                  "co2.45q15mod2.na.m", "co3.45q15mod2.na.m",
                  "co2.45q15mod3.na.m", "co3.45q15mod3.na.m",
                  "co1.45q15mod.a.f", "co2.45q15mod.a.f", "co3.45q15mod.a.f",
                  "co2.45q15mod2.a.f", "co3.45q15mod2.a.f",
                  "co2.45q15mod3.a.f", "co3.45q15mod3.a.f",
                  "co1.45q15mod.na.f", "co2.45q15mod.na.f", "co3.45q15mod.na.f",
                  "co2.45q15mod2.na.f", "co3.45q15mod2.na.f",
                  "co2.45q15mod3.na.f", "co3.45q15mod3.na.f",
                  "mu.age.m", "mu.age.f", "sd.age.m", "sd.age.f",
                  "africa.nums", "la.nums", "Mx.svd.scores.m", "Mx.svd.scores.f", "svd.coeffs.xp.m", "svd.coeffs.xp.f",
                  "ages", "user.par")
} else {
  save.list <- c("co1.e0mod.a.m", "co2.e0mod.a.m", "co3.e0mod.a.m", 
                 "co2.e0mod2.a.m", "co3.e0mod2.a.m",
                 "co2.e0mod3.a.m", "co3.e0mod3.a.m",
                 "co1.e0mod.na.m", "co2.e0mod.na.m", "co3.e0mod.na.m", 
                 "co2.e0mod2.na.m", "co3.e0mod2.na.m",
                 "co2.e0mod3.na.m", "co3.e0mod3.na.m",
                 "co1.e0mod.a.f", "co2.e0mod.a.f", "co3.e0mod.a.f",
                 "co2.e0mod2.a.f", "co3.e0mod2.a.f",
                 "co2.e0mod3.a.f", "co3.e0mod3.a.f",
                 "co1.e0mod.na.f", "co2.e0mod.na.f", "co3.e0mod.na.f",
                 "co2.e0mod2.na.f", "co3.e0mod2.na.f",
                 "co2.e0mod3.na.f", "co3.e0mod3.na.f",
                 "co1.5q0mod.a.m", "co2.5q0mod.a.m", "co3.5q0mod.a.m",
                 "co2.5q0mod2.a.m", "co3.5q0mod2.a.m",
                 "co2.5q0mod3.a.m", "co3.5q0mod3.a.m",
                 "co1.5q0mod.na.m", "co2.5q0mod.na.m", "co3.5q0mod.na.m", 
                 "co2.5q0mod2.na.m", "co3.5q0mod2.na.m",
                 "co2.5q0mod3.na.m", "co3.5q0mod3.na.m",
                 "co1.5q0mod.a.f", "co2.5q0mod.a.f", "co3.5q0mod.a.f",
                 "co2.5q0mod2.a.f", "co3.5q0mod2.a.f",
                 "co2.5q0mod3.a.f", "co3.5q0mod3.a.f",
                 "co1.5q0mod.na.f", "co2.5q0mod.na.f", "co3.5q0mod.na.f",
                 "co2.5q0mod2.na.f", "co3.5q0mod2.na.f",
                 "co2.5q0mod3.na.f", "co3.5q0mod3.na.f",
                 "co1.45q15mod.a.m", "co2.45q15mod.a.m", "co3.45q15mod.a.m",
                 "co2.45q15mod2.a.m", "co3.45q15mod2.a.m",
                 "co2.45q15mod3.a.m", "co3.45q15mod3.a.m",
                 "co1.45q15mod.na.m", "co2.45q15mod.na.m", "co3.45q15mod.na.m",
                 "co2.45q15mod2.na.m", "co3.45q15mod2.na.m",
                 "co2.45q15mod3.na.m", "co3.45q15mod3.na.m",
                 "co1.45q15mod.a.f", "co2.45q15mod.a.f", "co3.45q15mod.a.f",
                 "co2.45q15mod2.a.f", "co3.45q15mod2.a.f",
                 "co2.45q15mod3.a.f", "co3.45q15mod3.a.f",
                 "co1.45q15mod.na.f", "co2.45q15mod.na.f", "co3.45q15mod.na.f",
                 "co2.45q15mod2.na.f", "co3.45q15mod2.na.f",
                 "co2.45q15mod3.na.f", "co3.45q15mod3.na.f",
                 "co1.epmod.a.m", "co2.epmod.a.m", "co3.epmod.a.m", 
                 "co2.epmod2.a.m", "co3.epmod2.a.m",
                 "co2.epmod3.a.m", "co3.epmod3.a.m",
                 "co1.epmod.na.m", "co2.epmod.na.m", "co3.epmod.na.m", 
                 "co2.epmod2.na.m", "co3.epmod2.na.m",
                 "co2.epmod3.na.m", "co3.epmod3.na.m",
                 "co1.epmod.a.f", "co2.epmod.a.f", "co3.epmod.a.f",
                 "co2.epmod2.a.f", "co3.epmod2.a.f",
                 "co2.epmod3.a.f", "co3.epmod3.a.f",
                 "co1.epmod.na.f", "co2.epmod.na.f", "co3.epmod.na.f",
                 "co2.epmod2.na.f", "co3.epmod2.na.f",
                 "co2.epmod3.na.f", "co3.epmod3.na.f",
                 "mu.age.m", "mu.age.f", "sd.age.m", "sd.age.f",
                 "africa.nums", "la.nums", "Mx.svd.scores.m", "Mx.svd.scores.f", "svd.coeffs.xp.m", "svd.coeffs.xp.f",
                 "ages", "user.par")
}

  if(determ==FALSE){
  ## calculate residuals from fitting data with the various models and getting residuals
  print("Calculating residuals by age for determ=FALSE...")
  
  mu.age.m <- rep(0, Intervals)
  mu.age.f <- rep(0, Intervals)
  # e0 model
  mx.out.m <- matrix(NA, length(in.e0.m), Intervals)
  mx.out.f <- matrix(NA, length(in.e0.f), Intervals)
  for(i in 1:length(in.e0.m)){
    mx.out.m[i,] <- mortmod.e0(e0=in.e0.m[i], prev=in.prev[i], region=africa[i], sex=0, opt=TRUE, determ=TRUE, recal=save.output)
    mx.out.f[i,] <- mortmod.e0(e0=in.e0.f[i], prev=in.prev[i], region=africa[i], sex=1, opt=TRUE, determ=TRUE, recal=save.output)
  }
  resid.m <- unname(log(t(in.mx.m))-log(mx.out.m))
  resid.f <- unname(log(t(in.mx.f))-log(mx.out.f))
  sd.age.m <- sqrt(apply(resid.m^2,2,mean))
  sd.age.f <- sqrt(apply(resid.f^2,2,mean))
  } # if
  
  if(determ==TRUE){
    mu.age.m <- rep(0, Intervals)
    mu.age.f <- rep(0, Intervals)
    sd.age.m <- rep(0, Intervals)
    sd.age.f <- rep(0, Intervals)
  }

  if(!is.null(save.output)){
    print("Saving output parameters...")
    save(list=save.list,
         file=save.output, compress = "xz")
  }

  if(is.null(user.par)){
    mort.frame.m <- data.frame(prev.art.frame[,1:3], in.e0.m, in.5q0.m, in.45q15.m)
    names(mort.frame.m) <- c("country", "country_code", "midperiod", "input e0", "input 5q0", "input 45q15")
    mort.frame.f <- data.frame(prev.art.frame[,1:3], in.e0.f, in.5q0.f, in.45q15.f)
    names(mort.frame.f) <- c("country", "country_code", "midperiod", "input e0", "input 5q0", "input 45q15")
  } else {
    mort.frame.m <- data.frame(prev.art.frame[,1:3], in.e0.m, in.5q0.m, in.45q15.m, in.ep.m)
    names(mort.frame.m) <- c("country", "country_code", "midperiod", "input e0", "input 5q0", "input 45q15", "input user-defined")
    mort.frame.f <- data.frame(prev.art.frame[,1:3], in.e0.f, in.5q0.f, in.45q15.f, in.ep.f)
    names(mort.frame.f) <- c("country", "country_code", "midperiod", "input e0", "input 5q0", "input 45q15", "input user-defined")
  }
  
  print("Calibration Complete (See directory in save.output for .RData file.)") 
  return(list(svd.comps.m=Mx.svd.scores.m[,1:3], svd.comps.f=Mx.svd.scores.f[,1:3], weights.m=svd.coeffs.xp.m, weights.f=svd.coeffs.xp.f, params.a.m=coeff.array.a.m, params.a.f=coeff.array.a.f, params.na.m=coeff.array.na.m, params.na.f=coeff.array.na.f, input.epi=prev.art.frame, input.mort.m=mort.frame.m, input.mort.f=mort.frame.f))
} ## end of function