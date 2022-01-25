function [stop_time, tau_hat, statslog] = score(t,u,phi,Phi,B,T,mu0,A0,w,b,gamma,Ihat)
% t,u the hawkes process
% phi, Phi the influence kernel
% A0 mu0 the prechange parameters
% w the window size. should be approximately equal/larger than EDD
% b the threshold
% gamma the update rate
% M the topological structure of the network.
stop_time = -1;
tau_hat = -1;
K = ceil(w/gamma);
w = K*gamma;
D = length(mu0);
stop_time = -1;
tau_hat = -1;
n = length(t);
eta = [zeros(D,n);ones(1,n)];
mu0 = reshape(mu0,1,D);
t = reshape(t,1,n);
u = reshape(u,1,n);
U = u == (1:D)' ;
N = floor(T/gamma);
score = zeros(D+1,D,K);%score for each node [index2] during the [index3]-th grid
S = zeros(D+1,D);
statslog = zeros(2,N);
statslog(1,:) = (1:N)*gamma;
for i=1:D
    invI{i} = inv(Ihat{i});
end
i1 = 1; % the first event with time > (k-1)gamma-B
i2 = 1; % the first event with time > k*gamma-B
i3 = 1; % the first event with time > (k-1)gamma
i4 = 1; % the first event with time > k*gamma
t = [t Inf];
for k = 1:N
    %update i1 i2 i3 i4
    i1 = i2;
    i3 = i4;
    while t(i4)<k*gamma
        eta(1:D,i4) = sum(U(:,i1:i4-1).*phi(min(t(i4) - t(i1:i4-1),B)),2);
        i4 = i4+1;
    end
    while t(i2)<k*gamma - B
        i2 = i2+1;
    end
    %% update likelihood ratios
    lambda = (t(i1:i4-1)' > t(i3:i4-1)-B).*((i1:i4-1)'<(i3:i4-1))...
        .*phi(t(i3:i4-1) - t(i1:i4-1)');
    lambda0 = lambda.*(A0(u(i3:i4-1),u(i1:i4-1))'); %influence of event i1+[index1]-1 to event i3+[index2]-1, pre-change
    lambda0 = mu0(u(i3:i4-1)) + sum(lambda0); %lambda at event i3+[index2]-1, pre-change
    
    lag = [(k-1)*gamma k*gamma] - t(i1:i4-1)';
    lag = min(lag,B*ones(i4-i1,2));
    lag = max(lag,zeros(i4-i1,2));
    inte = Phi(lag(:,2)) - Phi(lag(:,1));
    linear_term = U(:,i1:i4-1)*inte;
    for i=1:D
        idx = find(u(i3:i4-1) == i);
        if k>K
            S(:,i) = S(:,i) - score(:,i,mod(k-1,K)+1);
        end
        if ~isempty(idx)
            score(:,i,mod(k-1,K)+1) = sum(eta(:,i3+idx-1)./lambda0(idx),2) - [linear_term;gamma];
        else
            score(:,i,mod(k-1,K)+1) = - [linear_term;gamma];
        end
        S(:,i) = S(:,i) + score(:,i,mod(k-1,K)+1);
        if k>=K
            statslog(2,k) = statslog(2,k) + S(:,i)'*invI{i}*S(:,i)/w;
        end
    end
    if statslog(2,k)>b
        stop_time = k*gamma;
        tau_hat = k*gamma - w;
        return;
    end
end

end

