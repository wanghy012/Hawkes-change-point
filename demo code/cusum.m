function [stopping,tau_hat,statslog] = cusum(t, u, phi, Phi, B, T, mu0, mu1, A0, A1, b, gamma)
% t,u the input event time and label. t should be sorted.
% phi the influence kernel
% Phi the cumulative phi
% B the truncation level
% T the total time
% mu0, A0 the prechange parameter
% mu1, A1 the postchange parameter
% b the detection threshold
% gamma the update rate


% tau the detected change-point. Returns -1 if there is no change.
% S the stopping time. Returns -1 if there is no change.

% update the CUSUM stat when event happens and at multiple of gamma

%assert(T>max(t),'Event time exceeds T');
%assert(prod(t(1:end-1)<t(2:end)),'Event time not sorted.')
if nargin <10
    gamma = 2/sum(mu0);
end

n = length(t);
t = reshape(t,1,n);
u = reshape(u,1,n);
t(end+1) = Inf;
Mu0 = sum(mu0);
Mu1 = sum(mu1);
D = length(mu0);
mu0 = reshape(mu0,1,D);
mu1 = reshape(mu1,1,D);
stopping = -1;
tau_hat = -1;
N = floor(T/gamma);

% stats: first row for index of change-point event , second row for log-likelihood ratio.
stats = [0;0];
statslog = zeros(3,N);
i1 = 1; % the first event with time > (k-1)gamma-B
i2 = 1; % the first event with time > k*gamma-B
i3 = 1; % the first event with time > (k-1)gamma
i4 = 1; % the first event with time > k*gamma
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
    %% compare stats for change-point tau<(k-1)gamma-B, preserve largest one
    l = sum(stats(1,:)<i1);
    [~,idx] = max(stats(2,1:l));
    stats = [stats(:,idx),stats(:,l+1:end)];
    %% add potential change-points
    stats(:,end+(1:i4-i3)) = [i3:i4-1;zeros(1,i4-i3)];
    %% update likelihood ratios
    lambda = (t(i1:i4-1)' > t(i3:i4-1)-B).*((i1:i4-1)'<(i3:i4-1))...
        .*phi(t(i3:i4-1) - t(i1:i4-1)');
    lambda0 = lambda.*(A0(u(i3:i4-1),u(i1:i4-1))'); %influence of event i1+[index1]-1 to event i3+[index2]-1, pre-change
    lambda0 = mu0(u(i3:i4-1)) + sum(lambda0); %lambda at event i3+[index2]-1, pre-change
    log_lambda0 = cumsum(log(lambda0)); %cumulative log-lambda till event i3+[index2]-1, prechange
    lambda1 = lambda.*(A1(u(i3:i4-1),u(i1:i4-1))'); %influence of event i1+[index1]-1 to event i3+[index2]-1, post-change
    for i=1:i4-i1-1
        lambda1(i4-i1-i,:) = lambda1(i4-i1-i,:) + lambda1(i4-i1-i+1,:);
    end
    lambda1(end+1,:) = zeros(1,i4-i3);
    lambda1 = lambda1 + mu1(u(i3:i4-1)); %lambda at event i3+[index2]-1 with change-point at event i1+[index1]-2
    log_lambda1 = log(lambda1);
    for i=1:i4-i3-1
        log_lambda1(:,end-i) = log_lambda1(:,end-i) + log_lambda1(:,end-i+1);
    end %cumulative log-lambda from event i3+[index2]-1 with change-point at event i1+[index1]-2
    lag = [(k-1)*gamma t(i3:i4-1) k*gamma]-t(i1:i4-1)';
    lag = min(lag,B*ones(i4-i1,i4-i3+2));
    lag = max(lag,zeros(i4-i1,i4-i3+2));
    inte = Phi(lag(:,2:end)) - Phi(lag(:,1));
    SA0 = sum(A0);
    SA1 = sum(A1);
    inte0 = sum(SA0(u(i1:i4-1))'.*inte) + Mu0*[(t(i3:i4-1)-(k-1)*gamma),gamma]; %integral of influence from (k-1)gamma till event time t(i3+[index2]-1),prechange
    inte1 = SA1(u(i1:i4-1))'.*inte; %integral of influence by event i1+[index1]-1 from (k-1)gamma till event time t(i3+[index2]-1),post-change
    for i=1:i4-i1-1
        inte1(i4-i1-i,:) = inte1(i4-i1-i,:) + inte1(i4-i1-i+1,:);
    end 
    inte1(end+1,:) = zeros(1,i4-i3+1);
    inte1 = inte1 + Mu1*[(t(i3:i4-1)-(k-1)*gamma),gamma]; %integral of lambda with change-point at event i1+[index1]-2 from (k-1)gamma till event time t(i3+[index2]-1),post-change
    %update prechange loglikelihoods
    if i4>i3
        stats(2,:) = stats(2,:) - log_lambda0(end);
    end
    stats(2,:) = stats(2,:) + inte0(end);
    %update post-change loglikelihoods
    %log-lambda part
    if i4>i3
        stats(2,i3-i1+2:end) = stats(2,i3-i1+2:end) + log_lambda0;
        stats(2,1:i3-i1+1) = stats(2,1:i3-i1+1) + log_lambda1(1:i3-i1+1,1)';
        stats(2,i3-i1+2:end-1) = stats(2,i3-i1+2:end-1) + diag(log_lambda1(i3-i1+2:end-1,2:end))';
    end
    %integral part
    stats(2,:) = stats(2,:) -  inte1(:,end)';
    stats(2,i3-i1+2:end) = stats(2,i3-i1+2:end) - inte0(1:end-1) + diag(inte1(i3-i1+2:end,1:end-1))';
    %stats updated. record statslog
    [~,idx] = max(stats(2,:));
    statslog(:,k) = [k*gamma;stats(2,idx);stats(1,idx)];
    %check for stopping rule
    if stats(2,idx) > b
        stopping = k*gamma;
        tau_hat = stats(1,idx);
        return;
    end
end
end
           
    
    
    
    
   



