function [EDD_stats] = simulate_EDD(config_file,stat_name,N)
% simulate the EDD using N sample paths.
% Load from config_file, which contains the following:
% A0 the true and specified in CUSUM pre-change parameter
% A1 the true post-change parameter
% A1p the specified post-change parameter in CUSUM. May not equals to A1.
% mu0 the true and specified in CUSUM pre-change parameter
% mu1 the true post-change parameter
% mu1p the specified post-change parameter in CUSUM. May not equals to mu1.
% kappa the true change-point
% T the time horizon
% gamma the grid size
% Save the output 'stats_EDD' in save_file.
% EDD_stats: records each time the cusum hits a new highest for each
% sample path.
rng shuffle

switch stat_name(1:2)
    case {'cu','Cu','CU'}
        folder_name = 'cusum';
    case {'gl','GL'}
        folder_name = 'GLR';
    case {'sc','Sc'}  
        folder_name = 'score';
    case {'sh','Sh'}
        folder_name = 'shewhart';
    otherwise
        error('unexpected statistic\n');
end
if nargin < 3
    N = 400;
else if ~isnumeric(N)
        N = str2double(N);
    end
end
load(['configs/' config_file],'A0','A1','mu0','mu1','kappa','gamma',...
    'A1p','mu1p','EDD_T','B','w','Ihat');
phi = @fphi;
Phi = @cphi;
EDD_stats(N) = struct();
%% simulate
for i=1:N 
    r = rand()*gamma;
    [t,u] = simulate_hawkes(A0,A1,mu0,mu1,kappa+r,EDD_T);
    switch stat_name(1:2)
            case {'cu','Cu','CU'}
                [~,~,statslog] = cusum(t,u,phi,Phi,B,EDD_T,mu0,mu1p,A0,A1p,Inf,gamma);
            case {'gl','GL'}
                [~,~,statslog] = GLR(t,u,phi,Phi,B,EDD_T,mu0,A0,w,Inf,gamma);
            case {'sc','Sc'}  
                [~,~,statslog] = score(t,u,phi,Phi,B,EDD_T,mu0,A0,w,Inf,gamma,Ihat);
            case {'sh','Sh'}
                [statslog] = shewhart(t,EDD_T,w);
        end
    stat = [0;0];
    for k=1:size(statslog,2)
        if statslog(2,k)>stat(2,end)
            stat(:,end+1) = statslog(1:2,k);
        end
    end
    EDD_stats(i).summary = stat;
    EDD_stats(i).kappa = kappa + r;
    EDD_stats(i).T = EDD_T;
end

%% save file
save_file = config_file;
if exist(['data/' folder_name '/' save_file],'file')
    EDD_stats2 = EDD_stats;
    clear EDD_stats
    N2 = N;
    load(['data/' folder_name '/' save_file],'EDD_stats');
    if exist('EDD_stats','var')
        EDD_stats(end+1:end+N2) = EDD_stats2;
    else
        EDD_stats = EDD_stats2;
    end
    save(['data/' folder_name '/' save_file],'EDD_stats','-append');
else
    save(['data/' folder_name '/' save_file],'EDD_stats','config_file');
end
fprintf(['file saved at data/' folder_name '/' save_file '\n']);
end

