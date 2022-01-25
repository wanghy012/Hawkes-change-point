function [ARL_stats,simulation_time] = simulate_ARL(config_file, stat_name, N)
% simulate the ARL using N sample paths.
% Load from config_file, which contains the following:
% A0 the true and specified in CUSUM pre-change parameter
% A1 the true post-change parameter
% A1p the specified post-change parameter in CUSUM. May not equals to A1.
% mu0 the true and specified in CUSUM pre-change parameter
% mu1 the true post-change parameter
% mu1p the specified post-change parameter in CUSUM. May not equals to mu1.
% T the simulation length. Should approximately equals to or larger than ARL.
% Save the output 'stats_ARL' in save_file.
% ARL_stats: records each time the cusum hits a new highest for each
% sample path.
% simulation_time: average time needed to compute cusum stats
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
load(['configs/' config_file],'A0','A1','mu0','mu1','B',...
    'gamma','A1p','mu1p','Ihat','w','ARL_T');
phi = @fphi;
Phi = @cphi;
ARL_stats(N) = struct();
simulation_time = 0;
%% simulate
for i=1:N
        [t,u] = simulate_hawkes(A0,A1,mu0,mu1,ARL_T,ARL_T);
        switch stat_name(1:2)
            case {'cu','Cu','CU'}
                tic
                [~,~,statslog] = cusum(t,u,phi,Phi,B,ARL_T,mu0,mu1p,A0,A1p,Inf,gamma);
                simulation_time = simulation_time + toc;
            case {'gl','GL'}
                tic
                [~,~,statslog] = GLR(t,u,phi,Phi,B,ARL_T,mu0,A0,w,Inf,gamma);
                simulation_time = simulation_time + toc;
            case {'sc','Sc'}  
                tic
                [~,~,statslog] = score(t,u,phi,Phi,B,ARL_T,mu0,A0,w,Inf,gamma,Ihat);
                simulation_time = simulation_time + toc; 
            case {'sh','Sh'}
                tic
                [statslog] = shewhart(t,ARL_T,w);
                simulation_time = simulation_time + toc;
        end
        stat = [0;0];
        for k=1:size(statslog,2)
            if statslog(2,k)>stat(2,end)
                stat(:,end+1) = statslog(1:2,k);
            end
        end
        ARL_stats(i).summary = stat;
        ARL_stats(i).T = ARL_T;
end
%% save file
save_file = config_file;
if exist(['data/' folder_name '/' save_file],'file')
    ARL_stats2 = ARL_stats;
    clear ARL_stats
    simulation_time2 = simulation_time;
    clear simulation_time
    N2 = N;
    load(['data/' folder_name '/' save_file],'simulation_time','ARL_stats');
    if exist('ARL_stats','var')
        ARL_stats(end+1:end+N2) = ARL_stats2;
        simulation_time = simulation_time + simulation_time2;
    else
        ARL_stats = ARL_stats2;
        simulation_time = simulation_time2;
    end
    save(['data/' folder_name '/' save_file],'ARL_stats','simulation_time','-append');
else
    save(['data/' folder_name '/' save_file],'ARL_stats','simulation_time','config_file');
end
fprintf(['file saved at data/' folder_name '/' save_file '\n']);
end



