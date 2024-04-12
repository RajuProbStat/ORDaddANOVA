################################################################################
#install.packages("Iso")
#install.packages("MASS")
#install.packages("NlcOptim")
library(Iso)
library(NlcOptim)
library(MASS)
################################################################################
#' Evaluate the observed test statistic value, critical value and P value for the test statistic W_1 from a given data
#'
#' @param a the number of levels of the row factor A
#' @param b the number of levels of the column factor B
#' @param size the level of significance
#' @param n matrix having a rows and b columns, and (i,j)-th element is the sample size in the (i,j)-th cell
#' @param mean matrix having a rows and b columns, and (i,j)-th element is the sample mean in the (i,j)-th cell
#' @param var matrix having a rows and b columns, and (i,j)-th element is the sample variance in the (i,j)-th cell
#'
#' @return observed, critical and P values along with the decision rule for the test W_1
#' @export
#'
#' @examples
#' W_1(3,2,0.05,rbind(c(45,50),c(45,50),c(45,50)),rbind(c(2.3709,2.6140),c(2.5236,2.7294),c(2.6444,2.7786)),rbind(c(0.3063,0.5912),c(0.2145,0.4684),c(0.4557,0.3797)))
W_1<-function(a,b,size,n,mean,var){
    func_MLE<-function(a,b,n,mean,var){
    ############################################################################
    ##MLEs of the parameters under the null parameter space
    mu_null_0<-mean(mean)
    var_null_0<-var
    var_null_1<-array(NA,dim=c(a,b))
    repeat
    {
      ##value of the log-likelihood function under the null space at (n-1)th iteration
      lik_null_0_s1<-0
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          lik_null_0_s1<- lik_null_0_s1-(n[i,j]/2)*log(2*pi*var_null_0[i,j])
        }
      }
      lik_null_0_s2<-0
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          lik_null_0_s2<-lik_null_0_s2-(n[i,j]/(2*var_null_0[i,j]))*(var[i,j]+(mean[i,j]-mu_null_0)^2)
        }
      }
      value_0_null<-lik_null_0_s1+lik_null_0_s2
      #print(value_0_null)
      w<-array(NA,dim=c(a,b))
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          w[i,j]<-(n[i,j]/var_null_0[i,j])
        }
      }
      sum1<-0
      sum2<-0
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          sum1<-sum1+w[i,j]*mean[i,j]
          sum2<-sum2+w[i,j]
        }
      }
      mu_null_1<-sum1/sum2                                   #updated value of mu
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          var_null_1[i,j]<-var[i,j]+(mean[i,j]-mu_null_1)^2  #updated values of variances
        }
      }
      ##########################################################################
      ##value of the log-likelihood function under the null space at n-th iteration
      lik_null_1_s1<-0
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          lik_null_1_s1<- lik_null_1_s1-(n[i,j]/2)*log(2*pi*var_null_1[i,j])
        }
      }
      lik_null_1_s2<-0
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          lik_null_1_s2<-lik_null_1_s2-(n[i,j]/(2*var_null_1[i,j]))*(var[i,j]+(mean[i,j]-mu_null_1)^2)
        }
      }
      value_1_null<-lik_null_1_s1+lik_null_1_s2
      if(abs(value_1_null-value_0_null)<0.000001)  #convergence criteria in terms of the log-likelihood value
      {
        break
      }
      mu_null_0<-mu_null_1
      var_null_0<-var_null_1
    }
    #print(c(mu_null_0,var_null_0))
    ############################################################################ MLEs of the parameters under the full parameter space
    mean_alpha<-rep(NA,a)
    mean_beta<-rep(NA,b)
    for(i in 1:a)
    {
      mean_alpha[i]<-mean(mean[i,])
    }
    for(j in 1:b)
    {
      mean_beta[j]<-mean(mean[,j])
    }
    ############################################################################ initialization of alpha, beta and sigma values
    alpha_dash_0_full<-mean_alpha
    beta_0_full<-mean_beta-mean(mean_beta)
    var_0_full<-var
    ############################################################################ maximization of beta when alpha and sigma values are fixed
    alpha_dash_1_full<-rep(NA,a)
    beta_1_full<-rep(NA,b)
    var_1_full<- array(NA,dim=c(a,b))
    repeat
    {
    ############################################################################ log-likelihood value at the (n-1)-th iteration
      lik_full_0_s1<-0
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          lik_full_0_s1<- lik_full_0_s1-(n[i,j]/2)*log(2*pi*var_0_full[i,j])
        }
      }
      lik_full_0_s2<-0
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          lik_full_0_s2<-lik_full_0_s2-(n[i,j]/(2*var_0_full[i,j]))*(var[i,j]+(mean[i,j]-alpha_dash_0_full[i]-beta_0_full[j])^2)
        }
      }
      value_0_full<-lik_full_0_s1+lik_full_0_s2
      #print(value_0_full)
      w1<-rep(NA,a)
      x1_dash<-rep(NA,a)
      for(i in 1:a)
      {
        sum1<-0
        sum2<-0
        for(j in 1:b)
        {
          sum1<-sum1+(n[i,j])/(var_0_full[i,j])
          sum2<-sum2+((n[i,j])/(var_0_full[i,j]))*mean[i,j]
        }
        w1[i]<-sum1
        sum3<-0
        for(j in 1:(b-1))
        {
          sum3<-sum3+((n[i,b]/var_0_full[i,b])-(n[i,j]/var_0_full[i,j]))*beta_0_full[j]
        }
        x1_dash[i]<-(1/sum1)*(sum2+sum3)
      }
      alpha_dash_1_full<-pava(x1_dash, w1, decreasing=FALSE, long.out=FALSE, stepfun=FALSE)  #updated value of alpha_dash
      ######################################################################### for the initialization of the beta values
       w2<-rep(NA,b)
       g<-rep(NA,b)
       for(j in 1:b)
       {
         sum1<-0
         sum2<-0
         for(i in 1:a)
         {
           sum1<-sum1+((n[i,j])/(var_0_full[i,j]))
           sum2<-sum2+((n[i,j])/(var_0_full[i,j]))*(mean[i,j]-alpha_dash_1_full[i])
         }
         w2[j]<-sum1
         g[j]<-(sum2/sum1)
       }
       beta_tilde_full<-pava(g, w2, decreasing=FALSE, long.out=FALSE, stepfun=FALSE)
       beta_0_initial_full<-beta_tilde_full-mean(beta_tilde_full)
       ########################################################################## maximization of beta's when alpha and sigmas are fixed
       objfun=function(y){
         beta<-rep(NA,b)
         for(j in 1:b)
         {
           beta[j]<-y[j]
         }
         sum_beta<-0
         for(j in 1:b)
         {
           sum1_beta<-0
           sum2_beta<-0
           for(i in 1:a)
           {
             sum1_beta<-sum1_beta+(n[i,j]/var_0_full[i,j])
             sum2_beta<-sum2_beta+(n[i,j]/var_0_full[i,j])*(mean[i,j]-alpha_dash_1_full[i])
           }
           sum_beta<-sum_beta+((sum2_beta/sum1_beta)-beta[j])^2*sum1_beta
         }
         return(sum_beta)
       }
       #constraint function
       confun=function(y){
         beta<-rep(NA,b)
         for(j in 1:b)
         {
           beta[j]<-y[j]
         }
         f=NULL
         f=rbind(f,sum(beta))
         return(list(ceq=f,c=NULL))
       }
       A_matrix<-array(NA,dim=c(b-1,b))
       for(i in 1:(b-1))
       {
         for(j in 1:b)
         {
           if(i==j)
           {
             A_matrix[i,j]<-1
           }
           else
           {
             if(j==i+1)
             {
               A_matrix[i,j]<--1
             }
             else
             {
               if(i==(b-1) & j==b)
               {
                 A_matrix[i,j]<--1
               }
               else
               {
                 A_matrix[i,j]<-0
               }
             }

           }
         }
       }
      B_matrix<-rep(0,b-1)
      beta_1_full<-solnl(beta_0_initial_full,objfun=objfun,confun=confun,A = A_matrix,B = B_matrix)$par #updated values of beta
      #print(beta_1_full)
      ########################################################################## maximization of sigma values when alpha and beta values are fixed
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          var_1_full[i,j]<-var[i,j]+(mean[i,j]-alpha_dash_1_full[i]-beta_1_full[j])^2                    #updated values of variances
        }
      }
      ########################################################################## log-likelihood value at the n-th iteration
      lik_full_1_s1<-0
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          lik_full_1_s1<- lik_full_1_s1-(n[i,j]/2)*log(2*pi*var_1_full[i,j])
        }
      }
      lik_full_1_s2<-0
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          lik_full_1_s2<-lik_full_1_s2-(n[i,j]/(2*var_1_full[i,j]))*(var[i,j]+(mean[i,j]-alpha_dash_1_full[i]-beta_1_full[j])^2)
        }
      }
      value_1_full<-lik_full_1_s1+lik_full_1_s2
      #print(value_1_full)
      if(abs(value_1_full-value_0_full)<0.000001)
      {
        break
      }
      alpha_dash_0_full<-alpha_dash_1_full
      beta_0_full<-beta_1_full
      var_0_full<-var_1_full
    }
    return(c(exp(value_1_null-value_1_full),var))
}
  LRT_value_var<-func_MLE(a,b,n,mean,var)
  B<-10000                              #number of bootstrap replicates
  set.seed(171197)
  LRT_boot<-rep(NA,B)
  for(l in 1:B)
  {
    mean_boot<-array(NA,dim=c(a,b))
    var_boot<-array(NA,dim=c(a,b))
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        data_boot<-rnorm(n[i,j],0,sqrt(LRT_value_var[i+j+1]))
        mean_boot[i,j]<-mean(data_boot)
        var_boot[i,j]<-var(data_boot)
      }
    }
    LRT_boot[l]<-func_MLE(a,b,n,mean_boot,var_boot)[1]#bootstrap values of the LRT statistic W_1
  }
  critical_value<-quantile(LRT_boot,size,names=FALSE) #critical_value of the LRT statistic W_1
  LRT_value<- LRT_value_var[1]                        #observed value of the LRT statistic W_1
  P_value<-mean(LRT_boot< LRT_value)                  #P_value for W_1
  cat("Test statistic value, critical_value and P_value for the LRT statistic W_1 are respectively \n ",LRT_value, critical_value, P_value)
  if(LRT_value<critical_value)
  {
    cat("\n W_1 rejects the null hypothesis\n")
  }
  else
  {
    cat("\n W_1 do not reject the null hypothesis\n")
  }
  cat("*********************************************\n")
}
################################################################################
#' Evaluate the observed test statistic values, critical values and P values for the test statistics W_2 and W_3 from a given data
#'
#' @param a the number of levels of the row factor A
#' @param b the number of levels of the column factor B
#' @param size the level of significance
#' @param n matrix having a rows and b columns, and (i,j)-th element is the sample size in the (i,j)-th cell
#' @param mean matrix having a rows and b columns, and (i,j)-th element is the sample mean in the (i,j)-th cell
#' @param var matrix having a rows and b columns, and (i,j)-th element is the sample variance in the (i,j)-th cell
#'
#' @return observed, critical and P values along with the decision rules for the tests W_2 and W_3
#' @export
#'
#' @examples
#' W_2_W_3(3,2,0.05,rbind(c(45,50),c(45,50),c(45,50)),rbind(c(2.3709,2.6140),c(2.5236,2.7294),c(2.6444,2.7786)),rbind(c(0.3063,0.5912),c(0.2145,0.4684),c(0.4557,0.3797)))
W_2_W_3<-function(a,b,size,n,mean,var)
{
  B<-10000                                #number of bootstrap replicates
  set.seed(171197)
  mean_alpha<-rep(NA,a)
  mean_beta<-rep(NA,b)
  mean_boot<-array(NA,dim=c(a,b))
  var_boot<-array(NA,dim=c(a,b))
  mean_alpha_boot<-rep(NA,a)
  mean_beta_boot<-rep(NA,b)
  W2_boot<-rep(NA,B)
  W3_boot<-rep(NA,B)
  for(i in 1:a)
  {
    mean_alpha[i]<-mean(mean[i,])
  }
  for(j in 1:b)
  {
    mean_beta[j]<-mean(mean[,j])
  }
  var1<-rep(0,a-1)
  var2<-rep(0,b-1)
  T1<-rep(NA,a-1)
  T2<-rep(NA,b-1)
  for(i in 1:a-1)
  {
    var1_0<-0
    for(j in 1:b)
    {
      var1[i]<-var1_0+(var[i+1,j]/n[i+1,j])+(var[i,j]/n[i,j])
      var1_0<-var1[i]
    }
    T1[i]<-(mean_alpha[i+1]-mean_alpha[i])/sqrt(var1[i]/b^2)
  }
  for(j in 1:b-1)
  {
    var2_0<-0
    for(i in 1:a)
    {
      var2[j]<-var2_0+(var[i,j+1]/n[i,j+1])+(var[i,j]/n[i,j])
      var2_0<-var2[j]
    }
    T2[j]<-(mean_beta[j+1]-mean_beta[j])/sqrt(var2[j]/a^2)
  }
  W2_ob<-max(c(max(c(T1)),max(c(T2))))    #observed value of W_2
  W3_ob<-max(c(min(c(T1)),min(c(T2))))    #observed value of W_3
  for(l in 1:B)
  {
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        data_boot<-rnorm(n[i,j],0,sqrt(var[i,j]))
        mean_boot[i,j]<-mean(data_boot)
        var_boot[i,j]<-var(data_boot)
      }
      mean_alpha_boot[i]<-mean(mean_boot[i,])
    }
    for(j in 1:b)
    {
      mean_beta_boot[j]<-mean(mean_boot[,j])
    }
    var1_boot<-rep(0,a-1)
    var2_boot<-rep(0,b-1)
    T1_boot<-rep(NA,a-1)
    T2_boot<-rep(NA,b-1)
    for(i in 1:a-1)
    {
      var1_0<-0
      for(j in 1:b)
      {
        var1_boot[i]<-var1_0+(var_boot[i+1,j]/n[i+1,j])+(var_boot[i,j]/n[i,j])
        var1_0<-var1_boot[i]
      }
      T1_boot[i]<-(mean_alpha_boot[i+1]-mean_alpha_boot[i])/sqrt(var1_boot[i]/b^2)
    }
    for(j in 1:b-1)
    {
      var2_0<-0
      for(i in 1:a)
      {
        var2_boot[j]<-var2_0+(var_boot[i,j+1]/n[i,j+1])+(var_boot[i,j]/n[i,j])
        var2_0<-var2_boot[j]
      }
      T2_boot[j]<-(mean_beta_boot[j+1]-mean_beta_boot[j])/sqrt(var2_boot[j]/a^2)
    }
    W2_boot[l]<-max(c(max(c(T1_boot)),max(c(T2_boot)))) #bootstrap values of W2
    W3_boot[l]<-max(c(min(c(T1_boot)),min(c(T2_boot)))) #bootstrap values of W3
  }
  critical_W2<-quantile(W2_boot,1-size,names=FALSE)     #critical value of W_2
  critical_W3<-quantile(W3_boot,1-size,names=FALSE)     #critical value of W_3
  P_value_W2<-mean(W2_boot>W2_ob)                       #P_value of W_2
  P_value_W3<-mean(W3_boot>W3_ob)                       #P_value of W_3
  cat("Test statistic value, critical_value and P_value for the test statistic W_2 are respectively \n",W2_ob, critical_W2, P_value_W2)
  if(W2_ob>critical_W2)
  {
    cat("\n W_2 rejects the null hypothesis")
  }
  else
  {
    cat("\n W_2 do not reject the null hypothesis")
  }
  cat("\n Test statistic value, critical_value and P_value for the test statistic W_3 are respectively \n",W3_ob, critical_W3, P_value_W3)
  if(W3_ob>critical_W3)
  {
    cat("\n W_3 reject the null hypothesis")
  }
  else
  {
    cat("\n W_3 do not reject the null hypothesis")
  }
  cat("\n*********************************************")
}
################################################################################














