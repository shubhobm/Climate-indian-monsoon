## misc_functions: miscellaneous functions for ease of use
## function to generate from multivariate normal
my.mvrnorm = function(n, mu, Sigma){
  p = length(mu)
  # compute square root of covariance matrix
  eo=eigen(Sigma, symmetric=TRUE)
  sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
  
  # generate random normals from runif by box-muller transform
  rnorm.vec = sqrt(-2*log(runif(n*p)))*cos(2*pi*runif(n*p))
  
  # generate sample matrix
  sample.matrix = matrix(rep(mu, n), nrow=n, byrow=T) +
    matrix(rnorm.vec, nrow=n, ncol=p)%*%sigma.sqrt
  return(sample.matrix)
}

# function giving matrix of ones
ones = function(m,n){
  matrix(1, nrow=m, ncol=n)
}

## calculate weighted projection quantile depth
EPQD = function(X, grid, nu=1e3){
  
  p = ncol(X)
  b = apply(X, 2, median)
  X0 = X - ones(nrow(X),1) %*% b
  grid0 = grid - ones(nrow(grid),1) %*% b
  
  ## get matrix of weighted PQDs for all points
  npt = dim(grid)[1]
  Fuxu.mat = matrix(0, nrow=npt, ncol=nu)
  
  # loop over nu pts on unit circle then take max
  for(iu in 1:nu){
    u = as.matrix(rnorm(p)); u = u/sqrt(sum(u^2))
    uecdf = ecdf(X0%*%u)
    Fuxu.mat[,iu] = uecdf(grid0%*%u)
  }
  #EPQD.vec = 1/(1+apply(abs(Fuxu.mat-.5), 1, max))
  EPQD.vec = apply(1-Fuxu.mat, 1, min)
  
  return(cbind(grid,EPQD.vec))
  
}

# compute Tyler's shape matrix, optionally with depth weights
TylerSig = function(X, tol=1e-5, maxit=100, weight=NULL){
  n = nrow(X); p = ncol(X)
  iSig = diag(rep(1,p))
  
  # whether to use depth weights
  if(is.null(weight)){
    weight = rep(1,n)
  }
  
  for(i in 1:maxit){
    iiSig = matrix(0,p,p)
    inv.iSig = solve(iSig)
    for(j in 1:n){
      xj = as.matrix(X[j,])
      iiSig = iiSig + weight[j]^2 * (xj %*% t(xj))/ as.numeric(t(xj) %*% inv.iSig %*% xj)
    }
    
    iiSig = iiSig/det(iiSig)^(1/p)
    if(norm(iSig - iiSig, type="F") < tol){
      break
    }
    else{
      iSig = iiSig
    }
  }
  
  iSig
}

PcaRank <- function(x, k=0, kmax=ncol(x), delta = 0.001,na.action = na.fail,
                    scale=FALSE, signflip=TRUE, trace=FALSE,
                    proj=50, ...)
{
  
  cl <- match.call()
  
  if(missing(x)){
    stop("You have to provide at least some data")
  }
  y <- data <- x <- as.matrix(x)
  n <- nrow(y)
  p <- ncol(y)
  
  ##
  ## verify and set the input parameters: k and kmax
  ##
  kmax <- max(min(floor(kmax), rankMM(x)),1)
  if((k <- floor(k)) < 0)
    k <- 0
  else if(k > kmax) {
    warning(paste("The number of principal components k = ", k, " is larger then kmax = ", kmax, "; k is set to ", kmax,".", sep=""))
    k <- kmax
  }
  if(k != 0)
    k <- min(k, ncol(data))
  else {
    k <- min(kmax, ncol(data))
    if(trace)
      cat("The number of principal components is defined by the algorithm. It is set to ", k,".\n", sep="")
  }
  ######################################################################
  
  ## VT::15.06.2010: introduce 'scale' parameter (instead of 'corr' in this case)
  ##  return the scale in the value object
  ##
  sc = vector('numeric', p) + 1
  if(scale == TRUE)
  {
    sc = apply(data, 2, "mad")
    for(i in 1:p) {
      data[, i] = data[, i]/sc[i]
    }
  }
  
  spa = spatial.median(data, delta)
  mu = spa$mu
  ep = spa$ep
  tt = matrix(mu, n, p, byrow=TRUE)
  data = data-tt
  depth = mdepth.RP(data, data, proj=proj)$dep
  depth = max(depth) - depth
  for(i in 1:n)
  {
    z = sqrt(sum((data[i,  ])^2))
    y[i,  ] = 0 * data[i,  ]
    if(z > ep)
    {
      y[i,  ] = depth[i] * (data[i,  ]  )/z
    }
  }
  
  ##out = princomp(y, scores = TRUE, cor = FALSE, na.action=na.action, subset = TRUE)
  ## no scaling - we have already scaled with MAD
  out = PcaClassic(y, k=k, kmax=kmax, scale=FALSE, signflip=signflip, ...)
  
  k <- out@k
  scores = data %*% out@loadings
  sdev = apply(scores, 2, "mad")
  names2 = names(sdev)
  orsdev = order(sdev)
  orsdev = rev(orsdev)
  sdev  = sdev[orsdev]
  scores  = scores[,orsdev, drop=FALSE]
  loadings = out@loadings[,orsdev, drop=FALSE]
  
  names(sdev)=names2
  dimnames(scores)[[2]]=names2
  dimnames(loadings)[[2]]=names2
  
  scale       <- sc
  center      <- as.vector(mu)
  scores      <- scores[, 1:k, drop=FALSE]
  loadings    <- loadings[, 1:k, drop=FALSE]
  eigenvalues <- (sdev^2)[1:k]
  
  ######################################################################
  names(eigenvalues) <- NULL
  if(is.list(dimnames(data)))
  {
    ##dimnames(scores)[[1]] <- dimnames(data)[[1]]
    rownames(scores) <- rownames(data)
  }
  dimnames(scores)[[2]] <- as.list(paste("PC", seq_len(ncol(scores)), sep = ""))
  dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))
  
  ## fix up call to refer to the generic, but leave arg name as <formula>
  cl[[1]] <- as.name("PcaLocantore")
  res <- new("PcaLocantore", call=cl,
             loadings=loadings,
             eigenvalues=eigenvalues,
             center=center,
             scale=scale,
             scores=scores,
             k=k,
             n.obs=n)
  
  ## Compute distances and flags
  res <- rrcov:::.distances(x, p, res)
  return(res)
}

## computes the spatial median
spatial.median <- function(x, delta)
{
  dime = dim(x)
  n=dime[1]
  p=dime[2]
  delta1=delta*sqrt(p)
  mu0=apply(x,2,median)
  h=delta1+1
  tt=0
  while(h>delta1)
  {
    tt=tt+1
    TT=matrix(mu0,n,p,byrow=TRUE)
    U=(x-TT)^2
    w=sqrt(apply(U,1,sum))
    w0=median(w)
    ep=delta*w0
    
    z=(w<=ep)
    w[z]=ep
    w[!z]=1/w[!z]
    w=w/sum(w)
    x1=x
    for(i in 1:n)
      x1[i,]=w[i]*x[i,]
    mu=apply(x1,2,sum)
    h=sqrt(sum((mu-mu0)^2))
    mu0=mu
  }
  out=list(mu=mu0,ep=ep)
  out
}