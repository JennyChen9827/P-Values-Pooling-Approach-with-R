mle_homo <- function(x,y,alternative=c("two.sided","less","great")){
#method5 MLE under homoscedasticity var(tumor)=var(normal)
  # approximated t score with df=n1
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
  t.star.mean =mean(x) #T.star_bar(combined tumor for n1+n2)
  n.star.mean = mean(y) #N_star_bar(combined normal for n1+n3)
  d.sd = sd(d)    # S(D)=sd(paried of n1)
  t1.sd=sd(x1)   #S(T1)=sd(paired tumor of n1)
  n1.sd=sd(y1)   #S(N1)=sd(paried normal of n1)
  t.star.sd=sd(x)  #S(T.star)=sd(combined tumor of n1+n2)
  n.star.sd=sd(y) #S(N.star)=sd(combined normal of n1+n3)
  nh = 2/(1/n2+1/n3)
  tn1.sd=cov(x1,y1)#S(TN1)=covariance of n1
  r = tn1.sd/(t1.sd*n1.sd)
  
    t.mean = mean(x2) #T_bar(independented tumor)
    n.mean = mean(y2) #N_bar(independented normal)
    t.sd = sd(x2)   # S(T)=sd(independed tumor of n2)
    n.sd = sd(y2)   # S(N)=sd(independed normal of n3)
    
    g.star = n1*(n1+n2+n3*r)*((n1+n2)*(n1+n3)-n2*n3*r^2)^(-1)
    f.star = n1*(n1+n3+n2*r)*((n1+n2)*(n1+n3)-n2*n3*r^2)^(-1)
    est.var = (t1.sd^2*(n1-1)+n1.sd^2*(n1-1)+(1+r^2)*(t.sd^2*(n2-1)+n.sd^2*(n3-1)))/(2*(n1-1)+(1+r^2)*(n2+n3-2))
    v1.star = est.var*(2*n1*(1-r)+(n2+n3)*(1-r^2))/((n1+n2)*(n1+n3)-n2*n3*r^2)
    ze = (f.star*(t1.mean-t.mean)-g.star*(n1.mean-n.mean)+t.mean-n.mean)/sqrt(v1.star)

    t.sd = 0   # S(T)=sd(independed tumor of n2)
    n.sd = 0   # S(N)=sd(independed normal of n3)
    t.mean = 0 #T_bar(independented tumor)
    n.mean = 0 #N_bar(independented normal)
    g.star = n1*(n1+n2+n3*r)*((n1+n2)*(n1+n3)-n2*n3*r^2)^(-1)
    f.star = n1*(n1+n3+n2*r)*((n1+n2)*(n1+n3)-n2*n3*r^2)^(-1)
    est.var = (t1.sd^2*(n1-1)+n1.sd^2*(n1-1)+(1+r^2)*(t.sd^2*(n2-1)+n.sd^2*(n3-1)))/(2*(n1-1)+(1+r^2)*(n2+n3-2))
    v1.star = est.var*(2*n1*(1-r)+(n2+n3)*(1-r^2))/((n1+n2)*(n1+n3)-n2*n3*r^2)
    ze = (f.star*(t1.mean-t.mean)-g.star*(n1.mean-n.mean)+t.mean-n.mean)/sqrt(v1.star)
  
  if (alternative=="two.sided"){
    p.value=2*pt(q=ze, df=n1,lower.tail=FALSE)
  }
  else if (alternative=="great"){
    p.value=pt(q=ze, df=n1,lower.tail=FALSE) 
  }
  
  else{
    p.value=pt(q=ze, df=n1,lower.tail=FALSE)
  }
  
  cat("n1_paired match = ",n1,"\n")
  cat("n2_tumor sample = ",n2,"\n")
  cat("n3_matched normal =",n3,"\n")
  cat("p.value=",p.value)
  cat("\n")
}
