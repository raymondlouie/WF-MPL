function s_MPL = estimate_MPL(mu,dt_array,single_mut,double_mut)
% estimate_MPL(mu,dt_array,single_mut,double_mut)
% 
% Estimate the selection coefficients using the MPL algorithm

% Inputs:
%       mu: mutation rate
%       dt_array: generation numbers corresponding to the observed frequencies
%       single_mut: single mutant frequencies. T x L matrix, where T is the 
%                   number of generations, and L is the number of residues
%       double_mut: double mutant frequencies. T x L x L matrix, where T is the 
%                   number of generations, and L is the number of residues
%
% Outputs:
%       s_MPL: selection coefficient estimates
%        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = length(dt_array);
L = size(single_mut,2);

sum_cov = zeros(L,L);
sum_b = zeros(L,1);

for ind_T = 2:T
    dt = dt_array(ind_T) - dt_array(ind_T-1);
    mut_mat = squeeze(double_mut(ind_T-1,:,:));
    
    prev_single_mut = single_mut(ind_T-1,:)';
    cov_mat = mut_mat - prev_single_mut*prev_single_mut';
    
    sum_cov = sum_cov+ dt*cov_mat;
    
    sum_b = sum_b + dt*(1 - 2*prev_single_mut);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate selection coefficients

reg_term=0.1;
b_vec = single_mut(end,:)' - single_mut(1,:)' - mu*sum_b;
s_MPL = (sum_cov+ reg_term*eye(L))\b_vec;