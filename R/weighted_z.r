weighted_z <- function(x,y,alternative=c("two.sided","less","great")){
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

  if (n1<3 & n2>3 & n3>3){
    p2=t.test(x2,y2,alternative = "great")$p.value
    zc=qnorm((1-p2),lower.tail = FALSE)
  }
  else if (n2<3 | n3<3){
    p1=t.test(d,alterntive="great")$p.value
    zc=qnorm((1-p1),lower.tail = FALSE)
  }
  else{
    w1=sqrt(2*n1)
    w2=sqrt(n2+n3)
    p1=t.test(d,alterntive="great")$p.value
    p2=t.test(x2,y2,alternative = "great")$p.value
    z1=qnorm((1-p1),lower.tail = FALSE)
    z2=qnorm((1-p2),lower.tail = FALSE)
    zc=(w1*z1+w2*z2)/(sqrt(w1^2+w2^2))
  }
  
 if (alternative=="two.sided"){
  pc=1-pnorm(zc,lower.tail = FALSE)
  p.value=ifelse(pc < 0.5, 2*pc, 2*(1-pc))
 }
 else if (alternative=="great"){
  cat("Great")
  p.value=1-pnorm(zc,lower.tail = FALSE)
 }
 else{
  cat("less")
  p.value=1-pnorm(zc,lower.tail = TRUE)
 }
  cat("n1_paired match = ",n1,"\n")
  cat("n2_tumor sample = ",n2,"\n")
  cat("n3_matched normal =",n3,"\n")
  cat("p.value=",p.value)
  cat("\n")
}

