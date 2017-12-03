function [single_mut double_mut] = WF_sim_traj(s,mu,L,N,p_init,dt_array)
% WF_sim_traj(s,mu,L,N,p_init,dt_array)
% 
% Generate the single and double mutant frequencies using a Wright-Fisher model

% Inputs:
%       s: selection coefficients (L x 1)
%       mu: mutation rate
%       L: number of residues
%       N: population size
%       p_init: initial genotype frequencies. T x 2^L matrix, where T is the 
%               number of observed generations
%       dt_array: generation numbers corresponding to the observed frequencies
%
% Outputs:
%       single_mut: single mutant frequencies
%       double_mut: single mutant frequencies
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


K = 2^L;
T_tot = max(dt_array);

U = sparse(zeros(K,L));
for ind_K=1:K
    U(ind_K,:) = fliplr(de2bi(ind_K-1,L));
end

g = U*s; % selection coefficients
w = g+ones(K,1); % genotype fitness


for indi=1:K
    for indj=1:K
        hamming_dist = sum(abs(U(indi,:)-U(indj,:)));
        Q(indi,indj) = (mu^hamming_dist)*((1-mu)^(L-hamming_dist));
    end
end

W_mat = Q*diag(w);

% initialize first time point
p=p_init;
mut_mat  = U'*diag(p)*U;

double_mut = zeros(T_tot,L,L);
single_mut = zeros(T_tot,L);


double_mut(1,:,:) = mut_mat;
single_mut(1,:)=diag(mut_mat);

for ind_T=2:T_tot
    
    % calculate mean vector
    f = W_mat*p;
    f = f/sum(f);
    
    % calculate frequencies in the next generation
    p = mnrnd(N,f)'/N;
    
    mut_mat  = U'*diag(p)*U;
    double_mut(ind_T,:,:) = mut_mat;
    single_mut(ind_T,:)=diag(mut_mat);
end

double_mut = double_mut(dt_array,:,:);
single_mut = single_mut(dt_array,:);