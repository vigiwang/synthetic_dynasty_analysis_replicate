#### Tom's Copula Functions #######

calcDiscreteCopulaBasic <- function(porig) {
  xdim <- dim(porig)
  r <- xdim[1]
  c <- xdim[2]
  row_targets <- rep(1/r,r)
  col_targets <- rep(1/c,c)
  oddsratio <- porig[1,1] * porig[2:r,2:c]       # the p00*pxy part
  x1 <- porig[2:r,1] %o% porig[1,2:c]            # the p0x*p0Y part
  oddsratio <- oddsratio / x1                    # oddsratio matrix
  if (all(dim(porig) == 2)) {    # If this is 2x2, simple formulae
    oddsratio <- (porig[1,1]*porig[2,2]) / (porig[1,2]*porig[2,1])
    x1 <- sqrt(oddsratio)
    x2 <- 2*(1+x1)
    pbar <- matrix(c(x1/x2,1/x2,1/x2,x1/x2),2,2,byrow=TRUE)
    YuleY <- (x1-1)/(x1+1)
  }
  else {   # Not 2x2, so must do ipfr
    pbar <- ipu_matrix(porig, row_targets, col_targets)  # create pbar via iterative proportional fittinng, per Geenens p 430
    #    pbar <- matrix(unlist(pbar),ncol=c)
    rown <- seq(r) - 1    # vector of row and columns index numbers (from zero)
    coln <- seq(c) - 1
    YuleY <- (4/((r-1)*(c-1)))*rown %*% (pbar %*% coln) - 1   # This should be sum(u*v*puv)
    YuleY <- YuleY * 3 * sqrt(((r-1)*(c-1))/((r+1)*(c+1)))
  }
  return(list(pbar=pbar,YuleY=YuleY))
}  

calcDiscreteCopula <- function(porig) {
  xdim <- dim(porig)
  r <- xdim[1]
  c <- xdim[2]
  row_targets <- rep(1/r,r)
  col_targets <- rep(1/c,c)
  oddsratio <- porig[1,1] * porig[2:r,2:c]       # the p00*pxy part
  x1 <- porig[2:r,1] %o% porig[1,2:c]            # the p0x*p0Y part
  oddsratio <- oddsratio / x1                    # oddsratio matrix
  if (all(dim(porig) == 2)) {    # If this is 2x2, simple formulae
    oddsratio <- (porig[1,1]*porig[2,2]) / (porig[1,2]*porig[2,1])
    x1 <- sqrt(oddsratio)
    x2 <- 2*(1+x1)
    pbar <- matrix(c(x1/x2,1/x2,1/x2,x1/x2),2,2,byrow=TRUE)
    YuleY <- (x1-1)/(x1+1)
  }
  else {   # Not 2x2, so must do ipfr
    pbar <- ipu_matrix(porig, row_targets, col_targets)  # create pbar via iterative proportional fittinng, per Geenens p 430
    #    pbar <- matrix(unlist(pbar),ncol=c)
    rown <- seq(r) - 1    # vector of row and columns index numbers (from zero)
    coln <- seq(c) - 1
    YuleY <- (4/((r-1)*(c-1)))*rown %*% (pbar %*% coln) - 1   # This should be sum(u*v*puv)
    YuleY <- YuleY * 3 * sqrt(((r-1)*(c-1))/((r+1)*(c+1)))
  }
  rowsums = rowSums(porig)
  colsums = colSums(porig)
  torig <- porig / matrix(rowsums,r,c)   # Transition matrix
  # Calculate long-term proportions as first eigenvector
  # The vector will be first column
  # Does this work for non-square?
  #  x1 <- eigen(t(torig))
  #  if (is.complex(x1$values)){   # Seems that the eigen doesn't always work ??
  if (r ==c) {
    x1 <- torig %*% torig %*% torig %*% torig %*% torig
    x1 <- x1 %*% x1
    x1 <- x1 %*% x1
    x1 <- x1 %*% x1   # I think this is 32 powers
    ssorig <- x1[1,]
  }
  else {ssorig <- 0}
  #  }
  #  else{
  #    x2 <- x1$vectors / matrix(colSums(x1$vectors),r,c,byrow=TRUE)
  #    ssorig <- x2[,1]
  #  }
  
  tbar <- pbar / matrix(rowSums(pbar),r,c)
  ssbar <- col_targets  
  if (r == c) {
    x1 <- eigen(t(tbar))
    eigbar <- x1$values
    x1 <- eigen(t(torig))
    eigorig <- x1$values
  }
  else{
    eigbar = 0
    eigorig = 0
  }
  
  if (r == c) {  
    x1 <- ssorig %*% tbar    # This is population fractions transformed by pbar from ss
    # The following calculates a "structural transition matrix" but take this out  
    #  pstruct <- ipu_matrix(porig, x1[1,],ssorig)  # create prob matrix to go from pbar-transformed fractions back to ss fraction
    #  tstruct <- pstruct / matrix(rowSums(pstruct),r,c)
    
    # Calculate a 'exchange' and 'structural' state distribution, for both original and steady-state:
    #   1. Transform from original state distribution to 'exchange' by transforming with tbar
    #   2. Calculate distance original -> exchange -> final. Weighted by starting state distribution
    #   3. Do this for both original and steady-state
    stdist_act <- matrix(0,3,c)
    stdist_act[1,] <- rowsums
    stdist_act[2,] <- colsums
    stdist_act[3,] <- rowsums %*% tbar
    rownames(stdist_act) <- c('start state dist','end state dist','exchange state dist')
    dist_act <- c(0,0)
    dist_act[1] <- ((stdist_act[3,]-stdist_act[1,])^2) %*% rowsums  # Sum of squared start -> exchange
    dist_act[2] <- ((stdist_act[3,]-stdist_act[2,])^2) %*% rowsums  # Sum of squared exchange -> final
    dist_act[3] <- ((stdist_act[1,]-stdist_act[2,])^2) %*% rowsums  # Sum of squared start -> final
    dist_act <- sqrt(dist_act)
    names(dist_act) <- c('start->exch','exch->end','start->end')
    stdist_ss <- matrix(0,3,c)
    stdist_ss[1,] <- ssorig
    stdist_ss[2,] <- ssorig
    stdist_ss[3,] <- ssorig %*% tbar
    rownames(stdist_ss) <- c('start state dist','end state dist','exchange state dist')
    dist_ss <- c(0,0,0)
    dist_ss[1] <- sum((stdist_ss[3,]-stdist_ss[1,])^2)  # Sum of squared start -> exchange
    dist_ss[2] <- sum((stdist_ss[3,]-stdist_ss[2,])^2)  # Sum of squared exchange -> final
    dist_ss[3] <- sum((stdist_ss[1,]-stdist_ss[2,])^2)  # Sum of squared start -> final
    dist_ss <- sqrt(dist_ss / r)   # Here the marginals are 1/r so dividing through by r is fine
    names(dist_ss) <- c('start->exch','exch->end','start->end')
  }
  else{
    dist_act = 0
    dist_ss = 0
    stdist_ss = 0
    stdist_act = 0
  }
  
  
  
  C1<-rep(0:(r-1), each = r)
  C2 <- c(rep(0:(c-1),c))
  x1 <- as.data.frame(pbar)
  Prob <- unlist(as.list(pbar))
  
  #  C1l<-rep(c(as.character(0:(r-1)),'T'), each = (r+1))
  #  C2l <- c(rep(c(as.character(0:(c-1)),'T'),c+1))
  C1l<-as.character(c(0:(r-1),'M'))
  C2l <- as.character(c(0:(c-1),'M'))
  C1<-rep(0:(r), each = (r+1))
  C2 <- c(rep(0:(c),c+1))
  x1 <- as.data.frame(pbar)
  x1$tc <- rowsums
  x2 <- c(colsums,0)
  x1[nrow(x1)+1,] <- x2
  Prob <- unlist(as.list(t(x1)))
  mydata<-data.frame(C1,C2,Prob)
  plotpbar <- ggplot(mydata, aes(x=C1, y=C2, size = Prob)) +
    geom_point(aes(colour=Prob),show.legend=FALSE) +
    theme_bw() + 
    scale_x_continuous(breaks = 0:(c),  # Specify where the ticks should be
                       labels = C2l) +   # Provide new labels
    scale_y_continuous(breaks = 0:r,labels=C1l) +
    scale_y_reverse(breaks = 0:r,labels=C1l) 
  
  Prob <- unlist(as.list(porig))
  x1 <- as.data.frame(porig)
  x1$tc <- rowsums
  x2 <- c(colsums,0)
  x1[nrow(x1)+1,] <- x2
  Prob <- unlist(as.list(t(x1)))
  C2 <- as.integer(C2)
  mydata<-data.frame(C1,C2,Prob)
  
  plotporig <- ggplot(mydata, aes(x=C1,y=C2, size = Prob)) +
    geom_point(aes(colour=Prob),show.legend=FALSE) +
    theme_bw() + 
    scale_x_continuous(breaks = 0:(c),  # Specify where the ticks should be
                       labels = C2l) +   # Provide new labels
    scale_y_continuous(breaks = 0:r,labels=C1l) +
    scale_y_reverse(breaks = 0:r,labels=C1l) 
  #    scale_y_continuous(breaks = seq(0, (r-1), by = 1)) +
  #    scale_x_continuous(breaks = seq(0, (c-1), by = 1)) +
  #        scale_y_reverse() 
  
  # Pearson correlation from original matrix, standard scores
  xx1 <- matrix(1:c,c,r)
  xx2 <- as.vector(xx1)
  xx1 <- t(matrix(1:r,r,c))
  xx3 <- as.vector(xx1)
  xxvalues <- matrix(c(xx2,xx3),nrow=length(xx2),ncol=2)
  xxProb <- as.vector(t(porig))   # Need to transpose to get in same order as xxvalues
  Pearson = cov.wt(xxvalues,wt=xxProb,cor=TRUE)$cor[1,2]
  # Pearson correlation from original matrix, midrank
  xr <- rowSums(porig)
  xc <- colSums(porig)
  x1 <- rep(0,c)
  x1[2:c] <- cumsum(xc)[1:(c-1)]
  x2 <- xc / 2    # middle
  x3 <- x2 + x1   # midranks
  xx1 <- matrix(x3,c,r)
  xx2 <- as.vector(xx1)
  x1 <- rep(0,r)
  x1[2:r] <- cumsum(xr)[1:(r-1)]
  x2 <- xr / 2    # middle
  x3 <- x2 + x1   # midranks
  xx1 <- t(matrix(x3,r,c))
  xx3 <- as.vector(xx1)
  xxvalues <- matrix(c(xx2,xx3),nrow=length(xx2),ncol=2)
  PearsonRank = cov.wt(xxvalues,wt=xxProb,cor=TRUE)$cor[1,2]
  
  
  return(list(porig = porig,pbar = pbar,torig=torig,tbar = tbar,ssorig = ssorig, ssbar = ssbar,
              rowsums=rowsums,colsums=colsums,YuleY = YuleY, plotpbar=plotpbar,
              plotporig=plotporig,oddsratio = oddsratio, stdist_act=stdist_act,
              stdist_ss=stdist_ss,dist_act=dist_act,dist_ss=dist_ss,eigbar=eigbar,eigorig=eigorig,
              Pearson=Pearson,PearsonRank=PearsonRank))
  #  return(list(pbar = pbar, gamma = gamma))
}  

calcDiscreteFreq <- function(forig) {
  porig <- forig / sum(forig)
  retlist <- calcDiscreteCopula(porig)
  return(retlist)
}  

calcGeenensAllPos <- function(porig) {
  xdim <- dim(porig)
  r <- xdim[1]
  c <- xdim[2]
  
  
  # Make the lower-case es (RxS matrixes whose all entries are zero except the i,jth)
  # Geenens p20
  elist <- list()
  for (i in 1:r) {
    xlist <- list()
    for (j in 1:c) {
      x1 = matrix(0,r,c,byrow=TRUE)
      x1[i,j] = 1
      xlist[[j]] <- x1
    }
    elist[[i]] <- xlist
  }
  
  # Make the upper-case Eloc basis functions
  # Geenens p20
  # Eloclist <- list()
  # Eortholist <- list()
  ElocMatrix <- matrix(0,r*c,(r-1)*(c-1))
  delta_alt <- rep(0,(r-1)*(c-1))
  for (i in 2:r) {
    # xlist1 <- list()
    # xlist2 <- list()
    for (j in 2:c) {
      x1 = elist[[i-1]][[j-1]] - elist[[i]][[j-1]] - elist[[i-1]][[j]] + elist[[i]][[j]] # 2 x 2; positive on the diagonal
      delta_alt[(i-2)*(r-1) + j-1] <- log((porig[i-1,j-1]*porig[i,j])/(porig[i-1,j]*porig[i,j-1])) # for this 16 values, store their coefficients: diag/non-diag
      #      xlist1[[j]] <- x1
      ElocMatrix[,(i-2)*(r-1) + j-1] <- as.vector(t(x1))
    }
    #   Eloclist[[i]] <- xlist1
  }
  
  ElocOrtho <- gramSchmidt(ElocMatrix)[[1]]   # make bases orthonomal
  
  logps = log(as.vector(t(porig)))
  delta = logps %*% ElocOrtho        # This should be the 'delta' of Geenens p 33
  DELTA = sqrt(sum(delta*delta))     # Frobenius norm, Eq 6.2
  x1 = sqrt(r*c)     # For nxn this will be n
  GeenensD = tanh(DELTA/x1)         # Speculate that this is the right scaling. Correct for 2x2 (Geenens p 34)
  GeenensD_alt = tanh(DELTA/(x1*((x1/2)^.25)))  # This adjusts the "n" factor by something like 25% for 5x5. Purely 
  # speculative, based on comparing vs Yule's rho for some cases from 2x2 to 5x5
  #  DELTA_alt = sqrt(sum(delta_alt*delta_alt))
  
  x1 <- calcDiscreteHellinger(porig)
  if (r == c)  {      # Square case. Create GeenensRho: signed depending on distance counter to co monotonic copula
    if (x1$HellingAssocNeg > x1$HellingAssocPos) {GeenensRho <- -GeenensD} else {GeenensRho <- GeenensD}
  }
  
  return(list(delta=delta,DELTA=DELTA,GeenensD=GeenensD,GeenensD_alt=GeenensD_alt,HellingGeen=x1$HellingGeen,
              HellingTSC=x1$HellingTSC,HellingAssoc=x1$HellingAssoc,GeenensRho=GeenensRho))
  #delta_alt=delta_alt,DELTA_alt=DELTA_alt))
}


calcGeenensStructZero <- function(porig) {
  izeros <- porig == 0
  if (sum(izeros) == 0 ) {       # There are no zeros, so do 'all positive' version
    x1 <- calcGeenensAllPos(porig)
    # return(list(delta=delta,DELTA=DELTA,GeenensD=GeenensD,GeenensD_alt=GeenensD_alt,HellingGeen=x1$HellingGeen,
    #           HellingTSC=x1$HellingTSC,HellingAssoc=x1$HellingAssoc,GeenensRho=GeenensRho))
    #                       #delta_alt=delta_alt,DELTA_alt=DELTA_alt))
    return(list(delta=x1$delta,DELTA=x1$DELTA,GeenensD=x1$GeenensD,R=NA,Q=NA,GeenensD_alt=x1$GeenensD_alt,
                ndim_red = NA,HellingGeen=x1$HellingGeen,HellingTSC=x1$HellingTSC,
                HellingAssoc=x1$HellingAssoc,GeenensRho=NA,
                R_2nd=NA,GeenD_2nd=NA,delta_2nd=NA))
  }
  
  else{
    xdim <- dim(porig)
    r <- xdim[1]
    c <- xdim[2]
    
    #  izeros <- which(porig == 0,arr.ind=TRUE)  # matrix with rows the indexes of the zeros
    #  izeros <- as.vector(t(izeros))    # make it into a long vector
    #  izeros <- porig != 0
    
    # Make the lower-case es (RxS matrixes whose all entries are zero except the i,jth)
    # Geenens p20
    elist <- list()
    for (i in 1:r) {
      xlist <- list()
      for (j in 1:c) {
        x1 = matrix(0,r,c,byrow=TRUE)
        x1[i,j] = 1
        xlist[[j]] <- x1
      }
      elist[[i]] <- xlist
    }
    
    # Make the upper-case E basis functions with pivot the first of lowest value for 'zero intersections'
    # Geenens p20
    EpivMatrix <- matrix(0,r*c,(r-1)*(c-1))  # Make the matrix full rank, then trim it later
    nonzeroi <- c()
    # Figure out the pivot point to use. The first of the lowest value for 'zero intersections'
    izeros <- which(porig == 0,arr.ind=TRUE)  # matrix with rows the indexes of the zeros
    czeros <- matrix(0,r,c)
    for (i in 1:dim(izeros)[1]) {   # add for each row and column corresponding to a zero in the matrix
      czeros[izeros[i,1],] <- czeros[izeros[i,1],] + 1
      czeros[,izeros[i,2]] <- czeros[,izeros[i,2]] + 1
    }
    izeros <- porig != 0   # Now make it just a matrix with "true" if zero
    x1 <- min(czeros)
    imin <- which(czeros == x1,arr.ind=TRUE)  # matrix with rows the indexes of the minimum values of 'zero intersections'
    iind <- 1:r
    ipiv <- imin[1,1]   # pivot using the first of the lowest 'zero intersections'
    iind <- iind[iind != ipiv]   # skip the index that will be the pivot
    jind <- 1:c
    jpiv <- imin[1,2]   # pivot using the first of the lowest 'zero intersections'
    jind <- jind[jind != jpiv]   # skip the index that will be the pivot
    xcol <- 0      # counter for columns that are non-zeroed
    for (i in iind) {
      # xlist1 <- list()
      # xlist2 <- list()
      for (j in jind) {
        xcol <- xcol + 1
        if (izeros[i,j] & izeros[i,jpiv] & izeros[ipiv,j]) {   # is this the right condition for excluding zeros??
          x1 = elist[[i]][[j]] - elist[[i]][[jpiv]] - elist[[ipiv]][[j]] + elist[[ipiv]][[jpiv]]
          #        x1 = elist[[i-1]][[j-1]] - elist[[i]][[j-1]] - elist[[i-1]][[j]] + elist[[i]][[j]]
          EpivMatrix[,xcol] <- as.vector(t(x1))  # put the e into appropriate column
          nonzeroi <- append(nonzeroi,xcol)  # Record the columns (basis vectors which are not 'knocked-out' by zeros
        }
      }
    }
    
    x1 <- (apply(EpivMatrix,1,max) == 0) & (apply(EpivMatrix,1,min) == 0)    # This gets the max of each row. 
    # Rows correspond to probabilities, and some probabilities are either zero or are 'knocked out'
    # of the odds ratio matrix / vector. a knocked-out probability will have zero for max, and can 
    # be (should be) removed 
    
    x2 <- as.vector(t(porig))
    logps = log(x2[!x1])       # Keep only the non-zero 
    EpivMatrix <- EpivMatrix[!x1,]   # Remove the rows corresponding to the zeros
    EpivMatrix <- EpivMatrix[,nonzeroi]    # This trims off the trailing (zero) columns (basis vectors removed by zeros)
    if (is.matrix(EpivMatrix)) {EpivOrtho <- gramSchmidt(EpivMatrix)[[1]]}   # make bases orthonomal}
    else {EpivOrtho <- EpivMatrix}       # the basis may have been reduced to a single element (vector) in which case already orthonormal
    
    
    # ------------------------
    # Now re-do everything for the second pivot basis, if there is more than 1
    if (dim(imin)[1]<2) {            # Only 1 good pivot, so set the 2nd versions to NA
      R_2nd <- NA
      GeenD_2nd <- NA
      delta_2nd <- NA
    }
    else {
      iind <- 1:r
      ipiv <- imin[2,1]   # pivot using the second of the lowest 'zero intersections'
      iind <- iind[iind != ipiv]   # skip the index that will be the pivot
      jind <- 1:c
      jpiv <- imin[2,2]   # pivot using the first of the lowest 'zero intersections'
      jind <- jind[jind != jpiv]   # skip the index that will be the pivot
      xcol <- 0      # counter for columns that are non-zeroed
      EzerMatrix <- matrix(0,r*c,(r-1)*(c-1))  # Make the matrix full rank, then trim it later
      nonzeroi <- c()
      ncol <- 0
      for (i in iind) {
        # xlist1 <- list()
        # xlist2 <- list()
        for (j in jind) {
          xcol <- xcol + 1
          if (izeros[i,j] & izeros[i,jpiv] & izeros[ipiv,j]) {   # is this the right condition for excluding zeros??
            x1 = elist[[i]][[j]] - elist[[i]][[jpiv]] - elist[[ipiv]][[j]] + elist[[ipiv]][[jpiv]]
            #        x1 = elist[[i-1]][[j-1]] - elist[[i]][[j-1]] - elist[[i-1]][[j]] + elist[[i]][[j]]
            EzerMatrix[,xcol] <- as.vector(t(x1))  # put the e into appropriate column
            nonzeroi <- append(nonzeroi,xcol)  # Record the columns (basis vectors which are not 'knocked-out' by zeros
          }
        }
      }
      
      EzerMatrix <- EzerMatrix[,nonzeroi]    # This trims off the trailing (zero) columns (basis vectors removed by zeros)
      x1 <- (apply(EzerMatrix,1,max) == 0) & (apply(EzerMatrix,1,min) == 0)    # This gets the max of each row. 
      # Rows correspond to probabilities, and some probabilities are either zero or are 'knocked out'
      # of the odds ratio matrix / vector. a knocked-out probability will have zero for max, and can 
      # be (should be) removed 
      x2 <- as.vector(t(porig))
      logps_2nd = log(x2[!x1])       # Keep only the non-zero 
      EzerMatrix <- EzerMatrix[!x1,]   # Remove the rows corresponding to the zeros
      EzerOrtho <- gramSchmidt(EzerMatrix)[[1]]   # make bases orthonomal
    }
    # ------------------------
    
    
    
    ndim_red =  (r-1)*(c-1) - dim(EpivMatrix)[2]
    R = (ndim_red) / ((r-1)*(c-1))    # Geenens's R(X,Y) p 33   
    delta = logps %*% EpivOrtho        # This should be the 'delta' of Geenens p 33
    DELTA = sqrt(sum(delta*delta))     # Frobenius norm, Eq 6.2
    x1 <- sqrt(r*c - ndim_red)        # This is also speculative
    Q = tanh(DELTA/x1)          # Geenens's Eq 6.3, Speculate that this is the right scaling. Correct for 2x2 (Geenens p 34)
    Q_alt = tanh(DELTA/(x1*((x1/2)^.25)))  # This adjusts the "n" factor by something like 25% for 5x5. Purely 
    # speculative, based on comparing vs Yule's rho for some cases from 2x2 to 5x5
    GeenensD = R + (1-R)*Q           # Eq 6.4 
    GeenensD_alt = R + (1-R)*Q_alt
    
    if (dim(imin)[1] >= 2) {   # There are 2 (or more) good pivot points, and the 2nd version uses 2nd pivot point
      ndim_red_2nd =  (r-1)*(c-1) - dim(EzerMatrix)[2]
      R_2nd = (ndim_red_2nd) / ((r-1)*(c-1))    # Geenens's R(X,Y) p 33   
      delta_2nd = logps_2nd %*% EzerOrtho        # This should be the 'delta' of Geenens p 33
      DELTA_2nd = sqrt(sum(delta_2nd*delta_2nd))     # Frobenius norm, Eq 6.2
      x1 <- sqrt(r*c - ndim_red)        # This is also speculative
      Q_2nd = tanh(DELTA_2nd/x1)          # Geenens's Eq 6.3, Speculate that this is the right scaling. Correct for 2x2 (Geenens p 34)
      GeenD_2nd = R_2nd + (1-R_2nd)*Q_2nd           # Eq 6.4 
    }
    
    x1 <- calcDiscreteHellinger(porig)
    if (r == c)  {      # Square case. Create GeenensRho: signed depending on distance counter to co monotonic copula
      if (x1$HellingAssocNeg > x1$HellingAssocPos) {GeenensRho <- -GeenensD} else {GeenensRho <- GeenensD}
    }
    
    return(list(delta=delta,DELTA=DELTA,GeenensD=GeenensD,R=R,Q=Q,GeenensD_alt=GeenensD_alt,
                ndim_red = ndim_red,HellingGeen=x1$HellingGeen,HellingTSC=x1$HellingTSC,
                HellingAssoc=x1$HellingAssoc,GeenensRho=GeenensRho,
                R_2nd=R_2nd,GeenD_2nd=GeenD_2nd,delta_2nd=delta_2nd))
  }
}


calcGeenensAllposFreq <- function(forig) {
  porig <- forig / sum(forig)
  retlist <- calcGeenensAllPos(porig)
  return(retlist)
}  


calcDiscreteHellinger <- function(porig) {
  xdim <- dim(porig)
  r <- xdim[1]
  c <- xdim[2]
  row_targets <- rep(1/r,r)
  col_targets <- rep(1/c,c)
  oddsratio <- porig[1,1] * porig[2:r,2:c]       # the p00*pxy part
  x1 <- porig[2:r,1] %o% porig[1,2:c]            # the p0x*p0Y part
  oddsratio <- oddsratio / x1                    # oddsratio matrix
  if (all(dim(porig) == 2)) {    # If this is 2x2, simple formulae
    oddsratio <- (porig[1,1]*porig[2,2]) / (porig[1,2]*porig[2,1])
    x1 <- sqrt(oddsratio)
    x2 <- 2*(1+x1)
    pbar <- matrix(c(x1/x2,1/x2,1/x2,x1/x2),2,2,byrow=TRUE)
    YuleY <- (x1-1)/(x1+1)
  }
  else {   # Not 2x2, so must do ipfr
    pbar <- ipu_matrix(porig, row_targets, col_targets)  # create pbar via iterative proportional fittinng, per Geenens p 430
    #    pbar <- matrix(unlist(pbar),ncol=c)
    rown <- seq(r) - 1    # vector of row and columns index numbers (from zero)
    coln <- seq(c) - 1
    YuleY <- (4/((r-1)*(c-1)))*rown %*% (pbar %*% coln) - 1   # This should be sum(u*v*puv)
    YuleY <- YuleY * 3 * sqrt(((r-1)*(c-1))/((r+1)*(c+1)))
  }
  
  pbarsqrt <- sqrt(pbar)
  margsqrt <- matrix(sqrt(1/(r*c)),r,c)
  Hsq_Dep <- sum((pbarsqrt - margsqrt)^2)/2
  BhatCoeff <- 1 - Hsq_Dep
  HellingGeen <- (2/BhatCoeff^2)*sqrt(BhatCoeff^4 + sqrt(4 - 3*BhatCoeff^4) - 2)
  # TSC 26-dec-24. Speculate that scaling Hellinger's distance to match the odds ratio and then Yule's rho for 2x2
  # Seems like an idea to try
  # But it turns out that it does not work beyond 2x2. 
  omegaHelling <- (2/(1-sqrt(1-4*((1-Hsq_Dep)^2-0.5)^2)  ) - 1)            # Convert from Hellinger distance to "odds ratio" (actually sqr root)
  # This is the appropriate scaling from Hellinger distance to odds ratio for the 2x2
  HellingTSC <- (omegaHelling - 1) / (omegaHelling + 1)  # This is the appropriate conversion to Yule's colligation for 2x2
  
  # Try Hellinger distance from counter- and co-montonicity copulas. 
  # Try the difference as a "correlation" BUT only for square matrix
  if (r == c) {
    hsqpos <- diag(r) / sqrt(r)     # square root of comonotonic
    hsqpos <- sum((pbarsqrt - hsqpos)^2)/2
    Hsq_Assoc <- hsqpos
    #   hsqpos <- hsqpos ^ (sqrt(r/2))    # Totally ad-hoc scaling to make fit better vs Yule's rho beyond 2x2
    HellingAssocPos <- 1 - 4*hsqpos + 2*hsqpos^2                   # Formula for scaling "negative" Hellinger distance to 2x2 Yule's rho
    #    omegapos <- (2/(1-sqrt(1-4*((1-hsqpos)^2-0.5)^2)  ) - 1)       # This is "sqrt(omega)" (OR) from the 2x2 case
    #    if (omegapos == 0) {HellingAssocPos <- 0} 
    #    else if (omegapos == Inf) {HellingAssocPos <- 1}
    #    else if (omegapos < 1) {HellingAssocPos <- (1/omegapos - 1) / 1/(omegapos + 1)}
    #    else {HellingAssocPos <- (1/omegapos - 1) / (1/omegapos + 1)}
    #    HellingAssocPos <- 1 - HellingAssocPos
    #      # Must invert for less than 1 (to get right sign)
    #    HellingAssocPos <- 1 - (2/xpos^2)*sqrt(xpos^4 + sqrt(4 - 3*xpos^4) - 2)
    hsqneg <- diag(r)[r:1,] / sqrt(r)     # square root of countermonotonic
    hsqneg <- sum((pbarsqrt - hsqneg)^2)/2
    #    hsqneg <- hsqneg ^ (sqrt(r/2))    # Totally ad-hoc scaling to make fit better vs Yule's rho beyond 2x2
    HellingAssocNeg <- -1 + 4*hsqneg - 2*hsqneg^2                   # Formula for scaling "negative" Hellinger distance to 2x2 Yule's rho
    
    # ------
    # For Helinger distance from "Exchange mobility copula" which is 1/n*(n-1) for off-diagonal, zero diagonal
    hsqExch <- matrix((1/(r*(r-1))),r,r)    # matrix with 1/n*(n-1) everywhere
    hsqExch <- hsqExch - diag(r)/(r*(r-1))   # subtracts out the diagonal. Now it's the "Exchange Mobility" copula
    hsqExch <- sqrt(hsqExch)         # make square root
    hsqExch <- sum((pbarsqrt - hsqExch)^2)/2     # Calculate the Hellinger distance from 'Exchange' to our copula
    # ------
    
    #    xneg <- (2/(1-sqrt(1-4*((1-xneg)^2-0.5)^2)  ) - 1)       # This is "sqrt(omega)" (OR) from the 2x2 case
    #    if (xneg == 0) {HellingAssocNeg <- 0} 
    #    else if (xneg == Inf) {HellingAssocNeg <- 1}
    #    else if (xneg < 1) {HellingAssocNeg <- (1/xneg - 1) / (1/xneg + 1)}
    #    else {HellingAssocNeg <- (xneg - 1) / (xneg + 1)}
    #    HellingAssocNeg <- HellingAssocNeg - 1
    # Must invert for less than 1 (to get right sign)
    #HellingAssocNeg <- 1 - (2/xneg^2)*sqrt(xneg^4 + sqrt(4 - 3*xneg^4) - 2)
    HellingAssoc_alt <- (HellingAssocPos + HellingAssocNeg)/2
    #    hsqpos <- hsqpos ^ (sqrt(r/2))    # Totally ad-hoc scaling to make fit better vs Yule's rho beyond 2x2
    #    hsqpos <- hsqpos ^ (-1.2279/log(1-(1/sqrt(r))))   # This scaling takes the independence copula back to rho = 0.0
    # (hsqpos = 0.2929 which is the value for rho=0 for 2x2)
    x1 <- log(1/sqrt(2))/log(1/sqrt(r))
    #    hsqpos <- 1 - (1-hsqpos)^x1                # This does the same as just above, but on the "reverese" scale from 1 to 0
    x2 <- 1 - 1/sqrt(2)                    # the zero corr for 2x2
    xr <- 1 - 1/sqrt(r)                    # the zero corr for rxr
    xexp <- 1 + .2*(r-2)/r                          # Totally ad-hoc adjustment to give more curvature 
    xexp <- 1
    xmax <- (r-mod(r,2))/r                     # This is the max Hillinger distance between co- and counter-monotonic copulas
    # When r is odd the diagonals overlap and the max is only (r-1)/r instead of 1
    xmax <- 1
    if (hsqpos < xr) {hsqpos = x2*(hsqpos/xr)^(1/xexp)} else {hsqpos = 1 - (1-x2)*((xmax - hsqpos)/(1 - xr))^xexp}  # This separately scales
    # For below zero and above zero
    # This scaling appears to work best
    HellingAssoc <-  1 - 4*hsqpos + 2*hsqpos^2  
    x1 = 1 - 1/sqrt(r)              # Scaling for quadratic that matches hsq=0 -> +1, hsq(1 / 1/sqrt(n)) = 0, hsq=1 -> -1
    xc = (2*x1-1)/(x1^2 - x1)
    xb = -2 - xc                     # Scaling for quadratic that matches hsq=0 -> +1, hsq(1 / 1/sqrt(n)) = 0, hsq=1 -> -1
    #    HellingAssoc_alt <-  1 + xb*hsqpos + xc*hsqpos^2  
    x2 <- 1 - 1/sqrt(2)                    # the zero corr for 2x2
    xr <- 1 - sqrt(r-1)/sqrt(r)                    # the zero corr for rxr
    HellingExchInd <- hsqExch           # Make an index which is 1 for no mobility (furthest from 'Exchange' copula) and 0 for closest
    if (hsqExch < xr) {hsqExch = x2*(hsqExch/xr)} else {hsqExch = 1 - (1-x2)*((1 - hsqExch)/(1 - xr))}  # This separately scales
    # For below zero and above zero
    # This scaling appears to work best
    HellingExchRho <-  1 - 2*(1-hsqExch)^2  
    #   HellingExchInd <- 1 - hsqExch           # Make an index which is 0 for no mobility (furthest from 'Exchange' copula) and 1 for closest
    
  }
  
  return(list(YuleY=YuleY,Hsq_Dep=Hsq_Dep,HellingGeen=HellingGeen,HellingTSC=HellingTSC,HellingAssocPos=HellingAssocPos,
              HellingAssocNeg=HellingAssocNeg,HellingAssoc=HellingAssoc,HellingAssoc_alt=HellingAssoc_alt,
              HellingExchRho=HellingExchRho,HellingExchInd=HellingExchInd,Hsq_Assoc=Hsq_Assoc))
}  


calcKullbackLeibler <- function(pA,pB){
  # define simple function for KullbackLeibler
  kl <- function(p1,p2){
    x1 <- sum(p1 == 0)
    xp1 <- p1   # Need to use xp1 because we need to maintain original zeros in p1
    xp2 <- p2
    if (x1 > 0){   # For any zeros in p1 replace by 1 (so log works ok - will be zero'd by zero prob)
      x2 <- p1 == 0
      xp1[x2] <- 1
      xp2[x2] <- 1  # replace any "duplicate" zeros since they will be knocked out (but don't want log(.) to blow up)
    }
    x1 <- sum(p2 == 0)  # Check only for zeros in 2nd matrix, because zeros for 1st are 'knocked out'
    if (x1 == 0) {return(sum(p1 * log(xp1/xp2)))}
    else {return(NaN)}
  }
  # Calculate matrixes IPF'd to opposite marginals
  x1 <- dim(pA)
  rdA <- x1[1]
  cdA <- x1[2]
  x1 <- dim(pB)
  rdB <- x1[1]
  cdB <- x1[1]
  KL_matrix <- matrix(NaN,3,5)
  rownames(KL_matrix) <- c('KL_AB','KL_BA','KL_Avg')
  colnames(KL_matrix) <- c('Total','Dependence','Marginals','%Depend','%Margin')
  if ((rdA == rdB) & (cdA == cdB)) {    # columns and rows match
    rsA <- rowSums(pA)
    csA <- colSums(pA)
    rsB <- rowSums(pB)
    csB <- colSums(pB)
    pstarA <- ipu_matrix(pA,rsB,csB)
    pstarB <- ipu_matrix(pB,rsA,csA)
    KL_matrix[1,1] <- kl(pA,pB)
    KL_matrix[1,2] <- kl(pA,pstarB)
    KL_matrix[1,3] <- kl(pstarB,pB)
    KL_matrix[2,1] <- kl(pB,pA)
    KL_matrix[2,2] <- kl(pB,pstarA)
    KL_matrix[2,3] <- kl(pstarA,pA)
    KL_matrix[3,] <- (KL_matrix[1,] + KL_matrix[2,])/2
    KL_matrix[,4] <- KL_matrix[,2] / KL_matrix[,1]
    KL_matrix[,5] <- KL_matrix[,3] / KL_matrix[,1]
  }
  return(KL_matrix)
}

sort_parm <- function(pmat){
  xn <- dim(pmat)
  if (xn[1] != xn[2]) {return(NA)}  # not a square matrix
  else {
    #xn <- xn[1]
    #xr <- matrix(rowSums(pmat),xn,xn)      # Don't need to do the full matrix, only the diagonal. Leaving this here for reference
    #xc <- matrix(colSums(pmat),xn,xn,byrow=TRUE)
    #xsort <- pmat / (xr*xc)     # Matrix of "sorting parameters" 
    x1 <- rowSums(pmat)*colSums(pmat)     # product of marginals (probabilities for independent, conditional on the marginals)
    #x2 <- diag(pmat) / x1             # diagonal 'sorting paramters' (ratio of actual to 'independent' probabilities)
    #x2 <- sum(x2*x1) / sum(x1)
    x2 <- sum(diag(pmat)) / sum(x1)    # this is the way Greenwwod et al. does it, and same (see 2 lines above) as Eika, Mogsted, Zafar method
    return(x2)
  }
}

or_matrix <- function(pmat){  # Function that returns matrix (r-1)(s-1) of odds ratios
  xdim <- dim(pmat)
  xoratio <- matrix(0,xdim[1]-1,xdim[2]-1)
  for (i in 2:xdim[1]){
    for (j in 2:xdim[2]){
      xoratio[i-1,j-1] <- (pmat[1,1]*pmat[i-1,j-1]) / (pmat[1,j-1]*pmat[i-1,1])
    }
  }
  return(xoratio)
}

altham_metric <- function(pmat,qmat){   # Function that returns Altham metric: 
  #   sum(i,j,l,m) sqrt[ (log(pij*plm*qim*qlj)/log(pim*plj*qij*qlm))^2]
  np = dim(pmat)
  nq = dim(qmat)
  alth_met <- 0
  if ((np[1] == nq[1]) & (np[2] == nq[2])) {
    for (i in 1:np[1]) {
      for (j in 1:np[2]) {
        for (l in 1:nq[1]) {
          for (m in 1:nq[2]) {
            alth_met = alth_met + log( (pmat[i,j]*pmat[l,m]*qmat[i,m]*qmat[l,j]) / (pmat[i,m]*pmat[l,j]*qmat[i,j]*qmat[l,m]) )^2
          }
        }
      }
    }
    alth_met <- sqrt(alth_met)
    return(alth_met)
  }
  else {return(NA)}
}

altham_index <- function(pmat){   # Altham index is Altham matric with independence copula as qmat
  np <- dim(pmat)
  qmat <- matrix(1/(np[1]*np[2]),np[1],np[2])
  return(altham_metric(pmat,qmat))
}





