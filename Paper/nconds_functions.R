

# nconds_functions.R
##############################################################


##################################################
make_circ_coord = function(t,x,ttot) {
  dt=(t[2]-t[1])*.45
  a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
  h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
  list(angles=a,heights=h)
}

circular_phase24H_histogram<-function(x,name,ttot){
  color_hist = rgb(0.6,0,0.2)
  br=0:ttot
  h=hist(x, br=br,plot=F)
  co=make_circ_coord(br[-1],h$counts,ttot)
  radial.plot(co$heights,co$angle,br[-1]-br[2],
              clockwise=T,start=pi/2,main=paste("",name),
              rp.type='p',poly.col=color_hist, xlab = "",ylab = "", show.grid.labels=0)
}
#################################################################
comb = function(n,k){
  factorial(n)/(factorial(k)*factorial(n-k))
  
}
##################################################################
nbt = function(x){
  l=length(which(x) == TRUE)
  l
}
#####################################
simply_it = function(x){
  a = 0
  for(i in x){
    a= paste(a,paste(which(x == as.numeric(i)), collapse = "",sep = ""), collapse = "", sep = "")  
  }
  a
}
######################################



creat_matrix_list = function(t, conds, n.co, period){
  
  my_matrix = list()
  
  c <- cos(2*pi*t/period)
  s <- sin(2*pi*t/period)   
  
  MAT <- cbind(rep(1,length(t)),c[1:length(t)],s[1:length(t)])
  GMAT <- matrix(NA,ncol=3*n.co, nrow =length(t))
  rownames(GMAT) <- conds
  colnames(GMAT) <- c(paste(c('u','a','b'),rep(1:n.co,each =3), sep = "."))
  
  it <- 1
  for(i in unique(rownames(GMAT))){
    GMAT[rownames(GMAT)==i,grep(paste0('.',it,'$'),colnames(GMAT))] = MAT[rownames(GMAT)==i,]
    it=it+1
  }
  
  vn = rep(F,n.co)
  for(i in 1:n.co){
    g = rep(F,n.co)
    g[1:i] = TRUE
    p = unique(combinat::permn(g))
    v = matrix(unlist(p),ncol = n.co,byrow = TRUE)
    vn = rbind(vn,v)
    
  }
  
 
  vn = vn[,rep(1:n.co,each=3)]
  vn[,seq(1,3*n.co,3)] = TRUE
  vn = data.frame(vn,row.names= NULL)
  vn[,dim(vn)[2] + 1]=(apply(vn,1,nbt)-n.co)/2
  colnames(vn) = c(paste(c('u','a','b'),rep(1:n.co,each =3), sep = "."),'nb_cycl')
  
  model = 1
  for(g in 0:n.co){
    
    
    nb_cycl =g
    com = expand.grid(rep(list(1:nb_cycl),nb_cycl))
    simply = apply(com,1,simply_it)
    poss =match(unique(simply),simply)
    com_l = com[poss,]
    pos = which(vn$nb_cycl == g)
    
    for(k in pos){
      if(g > 1){
        for(v in 1:nrow(com_l)){
          gmat = GMAT[,unlist(vn[k,-dim(vn)[2]])]
          ve = as.numeric(com_l[v,])
          id =1
          sa = ve
          while(length(ve) !=0){
            
            poc = which(sa == ve[1])
            po = which(ve ==ve[1])
            if(length(poc) !=1){
              poch =c(2*poc-1,2*poc)
              poch =poch[order(poch)]
              he = grep("[ab]",colnames(gmat))
              he = he[poch]
              pp=0
              for(z in 1:((length(he)-2)/2)){
                repl1 = which(gmat[,he[2*z+1]]!='NA')
                repl2 = which(gmat[,he[2*z+2]]!='NA')
                gmat[repl1,he[1]] = gmat[repl1,he[2*z+1]]
                gmat[repl2,he[2]] = gmat[repl2,he[2*z+2]]
                colnames(gmat)[he[1]]= paste(colnames(gmat)[he[1]],colnames(gmat)[he[2*z+1]],sep=',')
                colnames(gmat)[he[2]]= paste(colnames(gmat)[he[2]],colnames(gmat)[he[2*z+2]],sep=',')
                gmat[repl1,he[2*z+1]] =NA
                gmat[repl2,he[2*z+2]]=NA
                pp = pp+2
              }
              id = id+1
              ve = ve[-po]
            }else{
              ve = ve[-1]
            }
            
          }
          gmat[is.na(gmat)] =0
          del=which(apply(gmat,2,function(x) length(which(x == 0))) == length(t))
          if(length(del)!=0){
            gmat = gmat[,-del]
          }
          my_matrix[[model]] = gmat
          model = model + 1
        }
      }else{
        gmat = GMAT[,unlist(vn[k,-dim(vn)[2]])]
        gmat[is.na(gmat)] =0
        del =which(apply(gmat,2,function(x) length(which(x == 0))) == length(t))
        if(length(del)!=0){
          gmat = gmat[,-del]
        }
        my_matrix[[model]] = gmat
        model = model +1
      }
    }
    
  }
  
  
  
  return(my_matrix)
}
#################################3
compute_RSS = function(x, matX,kk){
  matX=matX[kk,]
  xx = solve(t(matX)%*%matX)
  y = xx %*% t(matX) %*% as.numeric(x)
  y = as.matrix(y)
  rownames(y)  = colnames(matX)
  RSS = t(x) %*% x -t(x) %*% matX %*% xx %*% t(matX) %*% x
  list(param=y,RSS = RSS)
}
#############################################
compute_AIC = function(A,n){
  
  p = length(A$param)
  AIC= n * log(A$RSS/n) + 2*p

  list(AIC = AIC, param = A$param)
  
} 
###############################################
extract_minAIC = function(x){
  ex = sapply(x, "[[", 1)
  pos = which.min(ex)
  
  A_w = 1/sum(exp((ex[pos]-ex)/2))
  c(model=pos,x[[pos]],AICW=A_w)
}
###############################################
compute_param = function(bf,n.co,period,conds){
  
  param = c(paste(c('u','a','b'),rep(1:n.co,each =3), sep = "."))
  paramout = rep(0,n.co*4)
  for(i in 1:n.co){
    u=bf$param[grep(param[3*i -2], rownames(bf$param))]
    a=bf$param[grep(param[3*i -1], rownames(bf$param))]
    b=bf$param[grep(param[3*i], rownames(bf$param))]
    if(length(u) ==0) u=0
    if(length(a) ==0) a=0
    if(length(b) ==0) b=0
    phase=period/(2*pi)*atan2(b,a)
    if(phase<0) phase=phase+period
    if(phase>period) phase=phase-period
    phase=phase%%period
    
    amp =2*sqrt(a^2+b^2)
    
    relamp=0.5*amp/u
    paramout[(1:4 + 4*(i-1))] = c(u,amp,relamp,phase)
  }
  paramout=c(bf$model,bf$AICW,paramout)
  names(paramout) = c('model','AICW',paste(c('mean','amp','relamp','phase'),rep(unique(conds),each =4), sep = "_"))
  paramout
}
######################################################################################
do_all = function(x,t,n.co,period,my_mat,conds){
  x = as.numeric(x)
  kk = which(is.na(x)==FALSE)
  x = x[kk]
  t = t[kk]
  n = length(x)
  
  my_fit = lapply(my_mat,compute_RSS, x = x,kk)
  my_AIC =lapply(my_fit,compute_AIC,n=n)
  
  bestfit = extract_minAIC(my_AIC)
  
  OUT = compute_param(bestfit,n.co,period,conds)
  OUT
}

InsertFitToMat <- function(fit, dat,n.co){
  fit_mat = matrix(unlist(fit), nrow = dim(dat)[1], byrow = T)
  rownames(fit_mat) = names(fit)
  pos = match(rownames(dat), rownames(fit_mat))
  data = data.frame(dat,fit_mat[pos,])
  names(data)[(dim(data)[2] - 4*n.co - 1):dim(data)[2]] = names(unlist(fit[[1]]))
  return(data)
}

