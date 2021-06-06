#' main function
#' 
#' This is the function allows you to calculate p value 
#' @param x  a (non-empty) numeric vector of data values(tumor sample)
#' @param y  a (non-empty) numeric vector of data values(matched normal sample)
#' @param method a character string indicating which method is to be computed. One of the "weighted"(default),"mle_hetero","mle_hoto","zcorr",or "modified_t".
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' 
#' @return statistics p-value of several methods for handing partially matched sample (tumor and matched normal sample) test by using meta-analysis framework.
#' @examples 
#' 
#' set.seed(123)
#' n=20
#' x=rnorm(n) 
#' y=rnorm(n)
#' x[1:5]=NA
#' y[7:10]=NA
#' pmchen(x,y,method="weighted_z",alternative="two.sided")
#' pmchen(x,y,method="mle_hetero",alternative="two.sided")
#' pmchen(x,y,method="mle_homo",alternative="two.sided")
#' pmchen(x,y,method="zcorr",alternative="two.sided")
#' pmchen(x,y,method="modified_t",alternative="two.sided")
#' pmchen(x,y)
#' pmchen(x,y,method="modified_t",alternative="great")
#' @export
pmchen <- function(x,y,method=c("weighted_z","mle_hetero","mle_homo","modified_t","zcorr"),alternative="two.sided"){
  
  data = data.frame(x,y)
  xy1 = data[complete.cases(data),]  #one sample test with n1
  xy2 = data[is.na(data$x)|is.na(data$y), ]
  x2 = na.omit(xy2$x)
  y2 = na.omit(xy2$y)
  x1 = xy1$x
  y1 = xy1$y
  d = y1-x1
  
  n1 = length(x1)
  n2 = length(x2)
  n3 = length(y2)
  

  if (n1<3 & n2<3 & n3<3){
    cat("Can't perform samples test since not enough observations, please enter more observations\n")
    return(-1)
  }

  else if(n2 <=3 | n3<=3){
    cat("Perform one sample t-test of before treatment and after treatment\n\n")
      p.value = t.test(d)$p.value
      cat("n1_paired match = ",n1,"\n")
      cat("p.value=",p.value)
      cat("\n")
  }

  else if(n1 <=3 ) {
    cat("Perform two samples t-test of before treatment and after treatment\n\n")
    p.value = t.test(x2,y2)$p.value
    cat("n2_tumor sample = ",n2,"\n")
    cat("n3_matched normal =",n3,"\n")
    cat("p.value=",p.value)
    cat("\n")
  }
  
  else if(missing(method)) {
    cat("Perform weighted z-test\n")
    weighted_z(x,y,alternative)
  }
  else if(method == "weighted_z") {
    cat("Perform weighted z-test\n")
    weighted_z(x,y,alternative)
  }
  else if(method == "mle_hetero"){
    p = var.test(x,y)$p.value
    if(p > 0.05){
      cat("Your input is equal variance, we perform this test with wrong result \n")
      mle_hetero(x,y,alternative)
    }
    
    else{
      cat("Perform the MLE based test statistic under heteroscedasticity\n")
      mle_hetero(x,y,alternative)
    }
  }
  else if(method == "mle_homo"){
    p = var.test(x,y)$p.value
    if( p < 0.05){
      cat("Your input x and y is unequal variance, we perform this test with wrong result\n")
      mle_homo(x,y,alternative)
    }
    else{
      cat("Perform the MLE based test statistic under homoscedasticity\n")
      mle_homo(x,y,alternative)
    }
  }
  else if(method == "modified_t"){
    cat("Perform the modified t-statistic\n")
    modified_t(x,y,alternative)
  }
  else if(method == "zcorr"){
    cat("Perform corrected Z-test\n")
    zcorr(x,y,alternative)
  }
}


