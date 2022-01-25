%%% This code generates examples of different stats (first row of 
%%% Figure 5) and simulates their ARL vs. EDD (Figure 6(a)) for the size-8
%%% network. All model setup used in Section 5 is included in 
%%% 'configs/', so you can also generate data and plots for size-100
%%% network by changing certain file names in the code.
%% set up
load('configs/w80.mat')

%% generate hawkes process with change-point
T = 1000;
kappa = 500;
[t,u] = simulate_hawkes(A0,A1,mu0,mu1,kappa,T);

%% plot stats examples
% cusum
[~,~,statslog1] = cusum(t, u, phi, Phi, B, T, mu0, mu1, A0, A1p, Inf, gamma);
figure(1)
plot(statslog1(1,:),statslog1(2,:),'k-','Linewidth',1)
hold on
plot([kappa kappa],ylim,'k--','Linewidth',1)
xlabel('t')
ylabel('CUSUM')

% GLR
[~,~,statslog2] = GLR(t,u,phi,Phi,B,T,mu0,A0,w,Inf,gamma);
figure(2)
plot(statslog2(1,:),statslog2(2,:),'k-','Linewidth',1)
hold on
plot([kappa kappa],ylim,'k--','Linewidth',1)
xlabel('t')
ylabel('GLR')

% score stat
[~,~,statslog3] = score(t,u,phi,Phi,B,T,mu0,A0,w,Inf,gamma,Ihat);
figure(3)
plot(statslog3(1,:),statslog3(2,:),'k-','Linewidth',1)
hold on
plot([kappa kappa],ylim,'k--','Linewidth',1)
xlabel('t')
ylabel('Score Stat')

% shewhart chart
w = 120; % set window length
[statslog4] = shewhart(t,T,w);
figure(4)
plot(statslog4(1,:),statslog4(2,:),'k-','Linewidth',1)
hold on
plot([kappa kappa],ylim,'k--','Linewidth',1)
xlabel('t')
ylabel('Shewhart')


%% Prepare to simulate ARL vs. EDD.
%%% Clear the pre-saved data if you wish to simulate them yourself.
%%% Otherwise skip the next section.

% rmdir data s
% mkdir data
% mkdir data/cusum
% mkdir data/GLR
% mkdir data/score
% mkdir data/shewhart
%  

%% simulate and save ARL and EDD data for [stat_name],
%%% where [stat_name] is among {cusum,GLR,score,shewhart}.

%%% simulating for CUSUM takes around 40 seconds.

stat_name = 'cusum';
simulate_ARL('w80.mat',stat_name,20);% simulate 20 sample paths for ARL.
simulate_EDD('w80.mat',stat_name,100);% simulate 100 sample paths for EDD.

%%% simulating for GLR takes around 14 mins.

stat_name = 'GLR';
simulate_ARL('w80.mat',stat_name,20);% simulate 20 sample paths for ARL.
simulate_EDD('w80.mat',stat_name,100);% simulate 100 sample paths for EDD.

%%% simulating for score stat takes around 70 seconds.

stat_name = 'score';
simulate_ARL('w80.mat',stat_name,20);% simulate 20 sample paths for ARL.
simulate_EDD('w80.mat',stat_name,100);% simulate 100 sample paths for EDD.

%%% simulating for shewhart chart takes around 10 seconds.

stat_name = 'shewhart';
simulate_ARL('w120.mat',stat_name,20);% simulate 20 sample paths for ARL.
simulate_EDD('w120.mat',stat_name,100);% simulate 100 sample paths for EDD.

%%% simulate ARL and EDD for CUSUM under misspecification. Takes around 160
%%% seconds.

stat_name = 'CUSUM';
simulate_ARL('misspec_dbl_scale.mat',stat_name,20);
simulate_EDD('misspec_dbl_scale.mat',stat_name,100);

simulate_ARL('misspec_dbl_topo.mat',stat_name,20);
simulate_EDD('misspec_dbl_topo.mat',stat_name,100);

simulate_ARL('misspec_half_scale.mat',stat_name,20);
simulate_EDD('misspec_half_scale.mat',stat_name,100);

simulate_ARL('misspec_half_topo.mat',stat_name,20);
simulate_EDD('misspec_half_topo.mat',stat_name,100);

%% estimate and save ARL and EDD for a series of detecting threshold b from saved data
get_ARL_EDD('data/cusum/w80.mat',3:0.5:9);
get_ARL_EDD('data/GLR/w80.mat',20:35);
get_ARL_EDD('data/score/w80.mat',90:5:170);
get_ARL_EDD('data/shewhart/w120.mat',800:8:920);
get_ARL_EDD('data/cusum/misspec_dbl_scale.mat',3:0.5:9);
get_ARL_EDD('data/cusum/misspec_half_scale.mat',3:0.5:9);
get_ARL_EDD('data/cusum/misspec_dbl_topo.mat',3:0.5:9);
get_ARL_EDD('data/cusum/misspec_half_topo.mat',3:0.5:9);

%% plot ARL vs EDD
figure(5)
hold on
load('data/cusum/w80.mat','ARL','EDD')
plot(log(ARL),EDD,'LineWidth',1);
load('data/GLR/w80.mat','ARL','EDD')
plot(log(ARL),EDD,'x-','LineWidth',1,'MarkerIndices',1:40);
load('data/score/w80.mat','ARL','EDD')
plot(log(ARL),EDD,'d-','LineWidth',1,'MarkerIndices',1:40);
load('data/shewhart/w120.mat','ARL','EDD')
plot(log(ARL),EDD,'o-','LineWidth',1,'MarkerIndices',1:40);
load('data/cusum/misspec_dbl_scale.mat','ARL','EDD')
plot(log(ARL),EDD,':','LineWidth',2,'Color',[29 206 225]/255);
load('data/cusum/misspec_half_scale.mat','ARL','EDD')
plot(log(ARL),EDD,'--','LineWidth',1.5,'Color',[29 206 225]/255);
load('data/cusum/misspec_dbl_topo.mat','ARL','EDD')
plot(log(ARL),EDD,':','LineWidth',2,'Color',[0 114 189]/255);
load('data/cusum/misspec_half_topo.mat','ARL','EDD')
plot(log(ARL),EDD,'--','LineWidth',1.5,'Color',[0 114 189]/255);
xlim(log([500 50000]))
xticks(log([500 1000 2e3 5e3 1e4 2e4 5e4]))
xticklabels({500 1000 2e3 5e3 1e4 2e4 5e4})
grid on
xlabel('log(ARL)')
ylabel('EDD')
legend('CUSUM','GLR','Score','Shewhart','Magn1','Magn2','Topo1','Topo2')


%% report average running time (seconds) in computing different stats 
%%% in ARL simulation

data_files = {'w80.mat','w80.mat','w80.mat','w120.mat'};
stats = {'cusum','GLR','score','shewhart'};
for i=1:4
    load(['data/' stats{i} '/' data_files{i}],'simulation_time','ARL_stats','config_file');
    load(['configs/' config_file],'ARL_T');
    fprintf('%s takes %f seconds to compute over T = %d.\n',...
        stats{i}, simulation_time/length(ARL_stats), ARL_T);
end




