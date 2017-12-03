# function which creates a binary representation of a decimal "number" with "noBits" bits
de2bi = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}

# 
WF_sim_traj <- function(s,mu,L,N,p_init,dt_array)
{
  
  K <- 2^L
  T_tot <- max(dt_array)
  
  # generate binary matrix
  U<-matrix(0,nrow=K,ncol=L)
  for (ind_count in 1:K){
    U[ind_count,] <- de2bi(ind_count-1,L)
  }
  
  g = U%*% s # genotype selection coefficients
  w = g + array(1,c(K,1))
  
  # mutation matrix
  Q<-matrix(0,nrow=K,ncol=K)
  for (indi in 1:K){
    for (indj in 1:K){
      hamming_dist = sum(abs(U[indi,]-U[indj,]));
      Q[indi,indj] = (mu^hamming_dist)*((1-mu)^(L-hamming_dist));
      
    }
  }
  
  Wmat <- Q%*%diag(c(w))
  
  # initialize first time point
  p  <- p_init
  
  mut_mat <-t(U) %*% diag(p) %*% U 
  
  double_mut<- array(0, c(T_tot,L, L))
  single_mut <- array(0,c(T_tot,L))
  
  double_mut[1,1:L,1:L] <- mut_mat
  single_mut[1,1:L] <- diag(mut_mat) 
  
  for (ite_T in 2:T_tot){
    
    # calculate mean vector
    f = drop(Wmat %*% p)
    f = f/sum(f)
    
    # calculate frequencies in the next generaiton
    p<-drop(rmultinom(1,N,f))/N
    
    mut_mat <-t(U) %*% diag(p) %*% U 
    double_mut[ite_T,1:L,1:L] <-mut_mat
    single_mut[ite_T,1:L] <- diag(mut_mat)
    
  }
  
  double_mut <- double_mut[dt_array,1:L,1:L]
  single_mut <- single_mut[dt_array,1:L]
  
  return(list(single_mut,double_mut))
  
}