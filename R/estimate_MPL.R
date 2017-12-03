estimate_MPL <- function(mu,dt_array,single_mut,double_mut)
{
  
  L <- ncol(single_mut)
  T <- length(dt_array)
  sum_cov <- array(0, c(L, L));
  sum_b <- array(0,c(1,L))
  
  for (ite_T in 2:T) {
    dt <-dt_array[ite_T] - dt_array[ite_T-1]
    
    prev_single_mut = single_mut[ite_T-1,]
    mut_mat = double_mut[ite_T-1,,]
    
    cov_mat <- mut_mat - prev_single_mut %o% prev_single_mut # covariance matrix
    sum_cov <- sum_cov + cov_mat*dt
    sum_b <- sum_b + dt*(1 - 2*prev_single_mut)  
  }
  
  reg_term <- 0.1
  
  if (L==1) {
    b_vec <- single_mut[T]  - single_mut[1] - mu*sum_b
    s_MPL <- b_vec/(sum_cov+reg_term) 
  } else{
    b_vec <- single_mut[T,]  - single_mut[1,] - mu*sum_b
    
    s_MPL <- solve(sum_cov+reg_term*diag(L),t(b_vec))
  }
  
  return(s_MPL)
  
}