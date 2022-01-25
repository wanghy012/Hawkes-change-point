function [ARL,EDD] = get_ARL_EDD(stats_file,b)
% Load ARL_stats, EDD_stats from stats_file and compute ARL,EDD at threshold b
% b : the array of thershold b
% ARL: the simulated ARL for each b
load(stats_file,'ARL_stats','EDD_stats');
b(end+1) = Inf;
%% ARL
N = length(ARL_stats);
S = zeros(N,length(b)); 
for i=1:N
    j=1;
    for k=1:size(ARL_stats(i).summary,2)
        while ARL_stats(i).summary(2,k)>b(j)
            S(i,j) = ARL_stats(i).summary(1,k);
            j = j+1;
        end
    end
    S(i,j:end) = -1;
end
ARL = zeros(1,length(b)-1);
for k=1:length(b)-1
    ARL(k) = (sum(S(S(:,k)>0,k)) + sum(S(:,k)<0)*ARL_stats(1).T)/(N - sum(S(:,k)<0));
end
%% EDD
N = length(EDD_stats);
S = zeros(N,length(b));
for i=1:N 
    j=1;
    for k=1:size(EDD_stats(i).summary,2)
        while EDD_stats(i).summary(2,k)>b(j)
            if EDD_stats(i).summary(1,k) > EDD_stats(i).kappa
                S(i,j) = EDD_stats(i).summary(1,k) - EDD_stats(i).kappa;
            else
                S(i,j) = 0;
            end
            j = j+1;
        end
    end
    S(i,j:end) = -1;
end
EDD = zeros(1,length(b)-1);
for i=1:length(b)-1
    if min(S(:,i))>=-0.5
        EDD(i) = sum(S(S(:,i)>0,i))/sum(S(:,i)>0);
    else
        EDD(i) = NaN;
    end
end

%% save file
b = b(1:end-1);
save(stats_file,'ARL','EDD','b','-append');
fprintf(['file saved at ' stats_file '\n']);
end

