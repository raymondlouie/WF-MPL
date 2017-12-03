close all;
clear all;

% This script generates the Wright-Fisher mutant trajectories, and estimates
% the selection coefficients from these trajectories using the MPL method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial parameters

N = 1000; % population size
L=7; % number of residues
dt_array =1:10:10000; % generations sampled
mu=1/N; % mutation probability

% selection coefficients (0, -, +)
% s = [zeros(floor(L/3),1) ; -abs(randn(floor(L/3),1))/100 ; abs(randn(L - 2*floor(L/3),1))/100];
s = [zeros(floor(L/3),1) ; -(10/N)*ones(floor(L/3),1) ; (10/N)*ones(L - 2*floor(L/3),1)];

no_runs=20; % number of WF trajectories to generate, for testing purposes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate WF runs and estimate the selection coefficients

K = 2^L; % number of genotypes

% Set initial genotype frequencies
p_init = ones(K,1);
p_init = p_init/sum(p_init);

% Variables used to calculate AUROC
ind_pos = find(s>0);
labels_pos = zeros(1,L);
labels_pos(ind_pos) = 1;

ind_neg = find(s<0);
labels_neg = zeros(1,L);
labels_neg(ind_neg) = 1;

% Initialization
s_MPL_array = zeros(no_runs,L);
nrmse_s_MPL = zeros(1,no_runs);
auc_s_est_pos = zeros(1,no_runs);
auc_s_est_neg = zeros(1,no_runs);

for ind_run=1:no_runs
    time_run= tic();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate trajectories and estimate selection coefficients
    
    [single_mut double_mut] = WF_sim_traj(s,mu,L,N,p_init,dt_array); % generate WF trajectories
    s_MPL_array(ind_run,:) = estimate_MPL(mu,dt_array,single_mut,double_mut); % estimate selection coefficients
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate performance
    
    % Calculate NRMSE (normalized root mean square error)
    nrmse_s_MPL(ind_run) = sqrt(sum((s_MPL_array(ind_run,:)' -s).^2)/sum(s.^2));
    
    % Calculate AUROC
    
    [~,~,~,auc_s_est_pos(ind_run)] = perfcurve(labels_pos,s_MPL_array(ind_run,:),1);
    [~,~,~,auc_s_est_neg(ind_run)] = perfcurve(labels_neg,-s_MPL_array(ind_run,:),1);
    
    time_run = toc(time_run);
    
    fprintf('Run number = %.0f/%.0f, run time = %f sec. \n',ind_run,no_runs,time_run);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Results

% Plot histogram

s_MPL_pos=s_MPL_array(:,find(s>0));
s_MPL_neg=s_MPL_array(:,find(s<0));
s_MPL_neu=s_MPL_array(:,find(s==0));

figure
h1 = histogram(s_MPL_pos(:));hold on
h2 = histogram(s_MPL_neg(:));
h3 = histogram(s_MPL_neu(:));
legend('Positive', 'Deleterious', 'Neutral','Location','Best');
ylabel('Probability of selection coefficient estimates');
xlabel('Selection coefficient estimates');
h1.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';
h1.BinWidth = h1.BinWidth /5;
h2.BinWidth = h1.BinWidth;
h3.BinWidth = h1.BinWidth;
for ind_L=1:L
    line([s(ind_L), s(ind_L)], ylim, 'LineWidth', 1, 'Color', 'k','LineStyle','--');
end

% Plot NRMSE
figure
boxplot(nrmse_s_MPL,'MPL');hold on;
ylabel(['NRMSE of the selection coefficients'])

% Plot AUROC
figure
boxplot([auc_s_est_pos' auc_s_est_neg'],{'Pos','Neg'});hold on;
ylabel(['AUROC'])