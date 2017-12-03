rm(list = setdiff(ls(), lsf.str())) # clear everything except functions

# load packages
require(ggplot2)
library(reshape)

source("WF_sim_traj.R")
source("estimate_MPL.R")

# function which creates a binary representation of a decimal number
de2bi = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}

########################################################################3
# set initial parameters

N <- 1000 # population size
L <- 7 # number of residues 
dt_array <- seq(1,10000,by=10) # generations sampled
mu <- 1/N # mutation probability

# selection coefficients (0, -, +)
s <- array(0,c(L,1))
# s[,1] <- c ( array(0,c(floor(L/3),1)) , -abs(rnorm(floor(L/3)))/100 , abs(rnorm(L-2*floor(L/3))/100) )
s[,1] <- c ( array(0,c(floor(L/3),1)) , array(-10/N,c(floor(L/3),1)), array(10/N,c(L-2*floor(L/3),1)))

no_runs <- 20 # number of WF trajectories to generate, for testing purposes

########################################################################3
# Generate WF runs and estimate the selection coefficients

K <- 2^L

p_init  <- drop(array(1/K,c(1,K)))
p_init <- p_init/sum(p_init)

ind_pos = which(s>0);
ind_neg = which(s<0);


s_MPL_array<- array(0,c(no_runs,L)) # MPL estimates


for (ind_run in 1:no_runs){
  
  start_time <- Sys.time()
  
  ########################################################################
  # Generate trajectories and estimate selection coefficients
  
  mut_list <- WF_sim_traj(s,mu,L,N,p_init,dt_array)
  
  single_mut<- mut_list[[1]]
  double_mut<- mut_list[[2]]
  

  s_MPL <- estimate_MPL(mu,dt_array,single_mut,double_mut)

  s_MPL_array[ind_run,1:L] <- drop(s_MPL)
  
  ########################################################################
  # Evaluate performance
  
  # Calculate NRMSE (normalized root mean square error)
  
  
  # Calculate AUROC
  
  end_time <- Sys.time()
  
  time_run <- end_time-start_time
  

  cat('Run number = ',ind_run,'/', no_runs , ' run time = ',time_run,' sec. \n')
  
  
}

id_pos_neg <- array("Neutral",c(L,1))
id_pos_neg[ind_pos] <- array("Positive",c(length(ind_pos),1))
id_pos_neg[ind_neg] <- array("Negative",c(length(ind_neg),1))

s_true_df = data.frame(s,id_pos_neg)
colnames(s_true_df) <- c("s","variable")

s_MPL_df = data.frame(s_MPL_array)
colnames(s_MPL_df) <- id_pos_neg
s_MPL_melt_df <-melt(s_MPL_df)

p_plot<- ggplot(s_MPL_melt_df, aes(x=value, fill=variable)) +
  geom_histogram(binwidth=abs(min(s))/20,alpha=.5, position="identity") +
  geom_vline(data=s_true_df, aes(xintercept=s),
             linetype="dashed", size=1) +  
  xlab("Selection coefficient estimates") +
  ylab("Probability of selection coefficient estimates")

x11()
print(p_plot)