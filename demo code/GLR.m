function [stop_time, tau_hat, statslog,n_EM_iter] = GLR(t,u,phi,Phi,B,T,mu0,A0,w,b,gamma,M)
% t,u the hawkes process
% phi, Phi the influence kernel
% B truncation level
% A0 mu0 the prechange parameters
% w the window size. should be approximately equal/larger than EDD
% b the threshold
% gamma the update rate
% M the topological structure of the network.

%assert(T>max(t),'Event time exceeds T');
%assert(prod(t(1:end-1)<t(2:end)),'Event time not sorted.')
Mu0 = sum(mu0);
K = ceil(w/gamma);
w = K*gamma;
N = floor(T/gamma);
D = length(mu0);
n_EM_iter = 0;
if nargin<12
    M = ones(D,D);
end
stop_time = -1;
tau_hat = -1;
n = length(t);
t = reshape(t,1,n);
u = reshape(u,1,n);
mu0 = reshape(mu0,1,D);
logl0 = zeros(1,N);
l0 = 0;
A1 = A0;
mu1 = mu0;
statslog = zeros(2,N);

i1 = 1; % the first event with time > (k-1)gamma-B
i2 = 1; % the first event with time > k*gamma-B
i3 = 1; % the first event with time > (k-1)gamma
i4 = 1; % the first event with time > k*gamma
i5 = 1; % the first event with time > k*gamma - w
t(end+1) = Inf;
for k = 1:N
    %update i1 i2 i3 i4
    i1 = i2;
    i3 = i4;
    while t(i4)<k*gamma
        i4 = i4+1;
    end
    while t(i2)<k*gamma - B
        i2 = i2+1;
    end
    while t(i5) < k*gamma - w
        i5 = i5+1;
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
    SA0 = sum(A0);
    inte0 = sum(SA0(u(i1:i4-1))'.*inte) + Mu0*gamma; %integral of influence from (k-1)gamma till event time t(i3+[index2]-1),prechange
    logl0(k) = sum(log(lambda0)) - inte0;
    if k>=K
        t2 = t(i5:i4-1) - k*gamma + w;
        u2 = u(i5:i4-1);
        [A1,mu1,logl,n_iter] = MLE_EM(t2,u2,phi,Phi,B,w,A1,mu1,M);
        n_EM_iter = n_EM_iter + n_iter;
    end
    
    %stats updated. record statslog
    l0 = l0 + logl0(k);
    if k>K
        l0 = l0 - logl0(k-K);
    end
    if k>=K
        statslog(:,k) = [k*gamma;sum(logl) - l0];
    end
    %check for stopping rule
    if statslog(2,k) > b
        stop_time = k*gamma;
        tau_hat = max(k*gamma - w,0);
        return;
    end
end
end

