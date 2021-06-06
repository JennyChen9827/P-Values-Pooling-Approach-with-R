mle_hetero <- function(x,y,alternative=c("two.sided","less","great")){
# method4 MLE under heteroscedasticity with t df=n1 var(tumor)=2*var(normal)
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
  t.mean = mean(x2) #T_bar(independented tumor)
  n.mean = mean(y2) #N_bar(independented normal)
  t.star.mean =mean(x) #T.star_bar(combined tumor for n1+n2)
  n.star.mean = mean(y) #N_star_bar(combined normal for n1+n3)
  nh = 2/(1/n2+1/n3)

  d.sd = sd(d)    # S(D)=sd(paried of n1)
  t1.sd=sd(x1)   #S(T1)=sd(paired tumor of n1)
  n1.sd=sd(y1)   #S(N1)=sd(paried normal of n1)
  t.sd = sd(x2)   # S(T)=sd(independed tumor of n2)
  n.sd = sd(y2)   # S(N)=sd(independed normal of n3)
  t.star.sd=sd(x)  #S(T.star)=sd(combined tumor of n1+n2)
  n.star.sd=sd(y) #S(N.star)=sd(combined normal of n1+n3)
  t1.sd=sd(x1)   #S(T1)=sd(paired tumor of n1)
  n1.sd=sd(y1)   #S(N1)=sd(paried normal of n1)
  tn1.sd=cov(x1,y1)#S(TN1)=covariance of n1
  r = tn1.sd/(t1.sd*n1.sd)
  
  f = n1*(n1+n3+n2*tn1.sd/(t1.sd^2))*(((n1+n2)*(n1+n3)-n2*n3*(r^2))^(-1))
  g = n1*(n1+n2+n3*tn1.sd/n1.sd^2)*((n1+n2)*(n1+n3)-n2*n3*r^2)^(-1)
  v1 = (f^2/n1+(1-f)^2/n2)*t1.sd^2*(n1-1)+(g^2/n1+(1-g)^2/n3)*n1.sd^2*(n1-1)-2*f*g*tn1.sd*(n1-1)/n1/(n1-1)
  zls = (f*(t1.mean-t.mean)-g*(n1.mean-n.mean)+t.mean-n.mean)/sqrt(v1)
  
  if (alternative=="two.sided"){
    p.value=2*pt(q=zls, df=n1,lower.tail=FALSE)
  }       
  else if (alternative=="great"){
    p.value=pt(q=zls, df=n1,lower.tail=FALSE) 
  }
  else{
    p.value=pt(q=zls, df=n1,lower.tail=TRUE) 
  }
  
  cat("n1_paired match = ",n1,"\n")
  cat("n2_tumor sample = ",n2,"\n")
  cat("n3_matched normal =",n3,"\n")
  cat("p.value=",p.value)
  cat("\n")
  }
  
