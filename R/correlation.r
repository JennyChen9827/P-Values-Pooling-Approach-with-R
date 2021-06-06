zcorr <- function(x,y,alternative=c("two.sided","less","great")){
  # method3 Z-test by correlation among the n1 match pair
  data = data.frame(x,y)
  xy1 = data[complete.cases(data),]  #one sample test with n1
  xy2 = data[is.na(data$x)|is.na(data$y), ]
  x2 = na.omit(xy2$x)
  y2 = na.omit(xy2$y)
  x1 = xy1$x
  y1 = xy1$y
  d = x1-y1
  x = c(x1,x2)
  y = c(y1,y2)
  
  n1 = length(d)
  n2 = length(x2)
  n3 = length(y2)
  d.mean = mean(d)
  t1.mean = mean(x1) #T1_bar(paried tumor)
  n1.mean = mean(y1) #N1_bar(paried normal)
  t.mean = mean(x2) #T_bar(independented tumor)
  n.mean = mean(y2) #N_bar(independented normal)
  t.star.mean =mean(c(x1,x2)) #T.star_bar(combined tumor for n1+n2)
  n.star.mean = mean(c(y1,y2)) #N_star_bar(combined normal for n1+n3)
  nh = 2/(1/n2+1/n3)
  
  d.sd = sd(d)    # S(D)=sd(paried of n1)
  t1.sd=sd(x1)   #S(T1)=sd(paired tumor of n1)
  n1.sd=sd(y1)   #S(N1)=sd(paried normal of n1)
  t.sd = sd(x2)   # S(T)=sd(independed tumor of n2)
  n.sd = sd(y2)   # S(N)=sd(independed normal of n3)
  t.star.sd=sd(c(x1,x2))  #S(T.star)=sd(combined tumor of n1+n2)
  n.star.sd=sd(c(y1,y2)) #S(N.star)=sd(combined normal of n1+n3)
  t1.sd=sd(x1)   #S(T1)=sd(paired tumor of n1)
  n1.sd=sd(y1)   #S(N1)=sd(paried normal of n1)
  tn1.sd=cov(x1,y1) #S(TN1)=covariance of n1
  r = tn1.sd/(t1.sd*n1.sd)
  
  zcor=(t.star.mean-n.star.mean)/sqrt(t.star.sd^2/(n1+n2)+n.star.sd^2/(n1+n3)-2*n1*tn1.sd/((n1+n2)*(n1+n3)))
  
  if (alternative=="two.sided"){
     p.value=2*pnorm(abs(zcor), lower.tail=FALSE)
  }
  else if (alternative=="great"){
    p.value=pnorm(zcor, lower.tail=FALSE)
    
  }
  
  else{
    p.value=pnorm(zcor, lower.tail=TRUE)
  }
  
  cat("n1_paired match = ",n1,"\n")
  cat("n2_tumor sample = ",n2,"\n")
  cat("n3_matched normal =",n3,"\n")
  cat("p.value=",p.value)
  cat("\n")
}

