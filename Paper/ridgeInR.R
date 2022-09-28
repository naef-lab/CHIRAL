
###############################
# Collection of R functions for ridge regression
# in the context of the MARA model E = NA
###############################

read.table.handlerows = function(infile){
  # Read table and handle potential duplicate rows
  dat = tryCatch({
    dat <- read.table(infile, header = T, row.names = 1)
  }, error = function(e){
    dat <- read.table(infile, header = T)
      rnames <- make.names(dat[, 1], unique = TRUE)
      rownames(dat) <- rnames
      dat <- dat[, 2:ncol(dat)]  # remove first column, it is in rownames
      return(dat)
  })
}

ridge.regression = function(N,E,lambda) {

  # Bayesian linear regression of E = N A
  # The input arguments N and E must be already centered
  
  Ns = fast.svd(N) # a list with entries: u, d, v
  rhs = crossprod(Ns$u,E)
  dia = Ns$d/(Ns$d^2 + nrow(N)*lambda)
  Ahat = sweep(Ns$v,2,dia,FUN='*') %*% rhs
  
  dimnames(Ahat) = list(colnames(N),colnames(E)) # so all derived matrices have the right row/column names
 
  Chi2 = colSums((E - N %*% Ahat)^2)    
  fov = 1-Chi2/colSums(E^2)
    
  C =  tcrossprod(sweep(Ns$v,2,1/(Ns$d^2 + nrow(N)*lambda),FUN='*'),Ns$v)
  AhatSE = sqrt(diag(C) %x% t(Chi2/nrow(E)))
  dimnames(AhatSE) = list(colnames(N), colnames(E))  # so all have right row/col names
  Zscore = Ahat/AhatSE
  combined.Zscore = sqrt(rowMeans(Zscore^2,na.rm=T))
  
  fit = list(Ahat=Ahat,
		Zscore=Zscore,
		combined.Zscore=combined.Zscore,
		fov=fov,
		AhatSE=AhatSE,lambda=lambda)
  return(fit)
}

###############################

optimize.lambda = function(N,E) {
  # optimize lambda by generalized cross-validation
  # The input arguments N and E must be already centered
  
  Ns = fast.svd(N) # a list with entries: u, d, v  
  rhs = crossprod(Ns$u,E)
  
  lambda.bnd = 10^c(-12,-6) * nrow(N) * ncol(N)

  gcv.error = function(lambda,E,Ns,rhs) {
    D = Ns$d^2/(Ns$d^2 + nrow(N)*lambda) # Hat matrix: H = VDV^t
    resid = E - sweep(Ns$u,2,D,FUN='*') %*% rhs
    GCV = sum((resid/(nrow(E)-sum(D)))^2)
    return(GCV)
  }

  opt = optimize(gcv.error,lambda.bnd,E,Ns,rhs)
  lambda.opt = opt$minimum
  gcv.opt = opt$objective
  
  return(list(lambda.opt = lambda.opt,
              gcv.opt = gcv.opt))
}

###############################

target.prediction = function(Ahat,centered.N,centered.E,N.colMeans) {
    # For every promotor p, motif m and N_pm > 0, calculate the
    # likelihood ratio of E_p with and without site m 
    # assuming Ahat stays un-changed by droping site m (checked, OK)   

    Ehat = centered.N %*% Ahat # predicted expression using all motif activities
    chi2 = rowSums((centered.E-Ehat)^2)
    chi2.mean = sum(chi2)/prod(dim(centered.E))
    # likelihood ratio: LR < 0 drop site, LR > 0 keep site
    LLR = matrix(NA,nrow(centered.E),ncol(centered.N),dimnames=list(rownames(centered.E),colnames(centered.N)))
    for (m in 1:ncol(centered.N)) {
      Nm = centered.N[,m]+N.colMeans[m]
      tidx = which(Nm>0) # target genes of motif m
      chi2.drop = rowSums((centered.E[tidx,] - (Ehat[tidx,] - Nm[tidx] %o% Ahat[m,]))^2) # the same using outer product
      LLR[tidx,m] = (chi2.drop-chi2[tidx])/chi2.mean
    }
    return(LLR) 
}

###############################
# helper functions
###############################

fast.svd = function(M,tol) {
  if (nrow(M) > 2*ncol(M)) {
    # more rows than cols
    s = svd(crossprod(M)) # svd(M^tM) = VD^2V^t 
    s$d = sqrt(s$d)
    s$u = M %*% sweep(s$v,2,s$d,FUN='/') 
   }
  else if (ncol(M) > 2*nrow(M)) {
    # more cols than rows
    s = svd(tcrossprod(M)) # svd(MM^t) = UD^2U^t
    s$d = sqrt(s$d)
    s$v = sweep(crossprod(M,s$u),2,s$d,FUN='/') 
  }
  else {
    # about the same rows/cols
    s = svd(M)
  }
  return(s)
}

###############################

center.rows = function(X) {
  # subtract the row-means of a matrix
  return(t(scale(t(X),scale=FALSE)))
}

###############################

center.cols = function(X) {
  # subtract the column-means of a matrix
  return(scale(X,scale=FALSE))
}

###############################
# advanced topics
###############################

ridge.regression.qr = function(N,E,lambda) {
  # ridge regression based on QR decomposition of N^tN
  
  QR = qr(crossprod(N) + diag(nrow(N)*lambda,ncol(N),ncol(N)))
  Ahat = qr.coef(QR,crossprod(N,E)) 
  Chi2 = colSums((E-N%*%Ahat)^2)   
  fov = 1-Chi2/colSums(E^2)
  
  # compute Z-score
  R = qr.R(QR) # the R matrx of the QR decomposition
  Rinv = backsolve(R,diag(1,nrow(R),nrow(R))) # R inverse
  C = Rinv %*% t(qr.Q(QR)) # (NtN + I*lambda)^-1
  AhatSE = sqrt(diag(C) %x% t(Chi2/nrow(E)))
  Zscore = Ahat/AhatSE
  combined.Zscore = sqrt(rowMeans(Zscore^2))  
  
  fit = list(Ahat=Ahat,Zscore=Zscore,combined.Zscore=combined.Zscore,fov=fov)
    
  return(fit)
}

###############################

compute.lml = function(N,E,NtN,NtE,lambda,idx=NULL) {
  # compute the log marginal likelihood log[P(E|N,lambda)] where
  # P(E|N,lambda) = int_sigma 1/sigma int_A P(E|N,A,sigma,lambda) P(A|lambda,sigma)  
  # This can be used to compare different models i.,e. site-count
  # matrices N as done below in the variable selection
  
  if (is.null(idx))
    idx = rep(TRUE,ncol(N))
  if (sum(idx)==0) { # empty model
    f = ncol(E)*nrow(E)/2
    return((f-1)*log(2) + lgamma(f) - f*log(sum(E^2)))
  }
  s = svd(NtN[idx,idx,drop=FALSE])
  dia = s$d + nrow(N)*lambda
  M = sweep(s$v,2,1/dia,FUN='*') %*% t(s$v)
  Ahat = M %*% NtE[idx,,drop=FALSE]  
  Chi2 = sum((E - tcrossprod(N[,idx,drop=FALSE],t(Ahat)))^2) # a bit faster
  f = ncol(E)*(nrow(E)-nrow(Ahat))/2
  log.det = sum(log(dia)) # log determinant 
  log.marg.lik = ncol(E)*0.5*(nrow(Ahat)*log(2*pi) - log.det) + (f-1)*log(2) + lgamma(f) - f*log(Chi2) 
  return(log.marg.lik)
}

###############################

variable.selection = function(N,E,start.set=NA) {
  # do a stepwise variable selection (add/drop variables) using the
  # Bayes factor (ratio of marginal likelihoods)

  NtN = crossprod(N)
  NtE = crossprod(N,E)
  # initialize
  if (all(is.na(start.set))) {
    active = rep(FALSE,ncol(N)) # inital model is empty
    lambda = 0
  } else {
    active = rep(FALSE,ncol(N))
    active[start.set] = TRUE
    lambda = optimize.lambda(N[,active,drop=FALSE],E)$lambda.opt
  } 
  lml = compute.lml(N,E,NtN,NtE,lambda,active) 
  ok = TRUE
  while (ok) {
    # determine lml of all neighboring models
    alml = rep(-Inf,ncol(N))
    for (n in 1:ncol(N)) {
      active[n] = !active[n] # add/drop
      alml[n] = compute.lml(N,E,NtN,NtE,lambda,active)
      active[n] = !active[n] # flip back
    }
    lbf = alml-lml # log Bayes factors
    # pick best neighboring model
    midx = which.max(lbf)
    ok = lbf[midx]>0
    if (ok) {
      # update model
      active[midx] = !active[midx]
      lml = alml[midx]
      cat(sprintf("%s\t%s\tlog[BF] = %g\n",ifelse(active[midx],"add","drop"),colnames(N)[midx],lbf[midx]))
      # update lambda
      lambda = optimize.lambda(N[,active,drop=FALSE],E)$lambda.opt
    } else {
      cat(sprintf("\nIteration stoped at max(log[BF]) = %g\nwith %d features included\n",lbf[midx],sum(active)))
    }
  }
  active = colnames(N)[active]
  return(active)
}

###############################

fit.promoter.scale = function(N,E,nmax=20,fov.eps=0.01) {
  # Fit optimal scaling factor for a given promotor, i.e. a model of
  # the form: E_ps = exp(beta_p) (\sum_m N_pm A_ms) by iteratively
  # estimating A (via ridge regression) and b_p. There are two
  # regularization factors: lambda1 for |A| and lambda2 for |exp(beta)|

  n = 0
  promoter.scales = rep(1,nrow(E))
  fov = c()
  while (n<nmax) {
    n = n + 1

    # apply promoter.scales and center 
    Nl = scale(sweep(N,1,promoter.scales,FUN='*'),scale=FALSE)

    # update lambda 1 + re-estiamte activities
    lambda1 = optimize.lambda(Nl,E)$lambda.opt
    r = ridge.regression(Nl,E,lambda1)
    fov[n] = mean(r$fov)
    print(sprintf('lambda1 = %g, FOV = %g',lambda1,fov[n]),quote=F)
    if ((n>1) && (fov[n]-fov[-1])/fov[n] > fov.eps) {
      break
    }
    # estimate promotor scaling: find a lambda2 to fix |promoter.scales| = ngene
    Ehat = N %*% r$Ahat # Ehat without promoter scaling 
    s1 = apply(E*Ehat,1,sum)
    s2 = apply(Ehat^2,1,sum)
    lambda2.root = function(lambda2,s1,s2,nprom) {
      f = (s1/(lambda2 + s2))
      f[f<0] = 0
      return(sum(f^2)-nrow(E))
    }
    lambda2 = uniroot(lambda2.root,c(0,1e2),s1,s2,nrow(E))$root
    print(sprintf('lambda2 = %g',lambda2),quote=F)
    
    # 2. estimate promoter scales with optimal lambda2
    promoter.scales = s1/(lambda2 + s2)
    promoter.scales[promoter.scales<0] = 0
  }
  return(list(promoter.scales=promoter.scales,Ahat=r$Ahat)) 
}


