function [A,mu,logl,n_iter] = MLE_EM(t,u,phi,Phi,B,T,A,mu,M)
% use A,mu as initialization
% compute the MLE estimator using data t,u,
% phi the influence kernel, Phi the cumulative
% B the truncation level
% total time T
% network topology M

%output: the MLE estimator A,mu, fitted log-likelihood logl.
%assert(T>max(t),'Event time exceeds T');
%assert(prod(t(1:end-1)<t(2:end)),'Event time not sorted.')
D = length(mu);
n_iter = 0;
A(A<1e-4) = 1e-4;
mu(mu<1e-2) = 1e-2;
if nargin<8
    M = ones(D,D);
end
n = length(t);
t = reshape(t,1,n);
u = reshape(u,1,n);
U = u == (1:D)' ;
valid = sum(U,2) > 0;
for i=1:D
    linear_term(i,1) = sum(Phi(max(T-t(u==i),B)));
end
linear_term(D+1,1) = T;
eta = zeros(D,n);
i1 = 1; % the first event with event time > t(i)-B
for i=1:n
    while t(i1)<t(i)-B
        i1 = i1 + 1;
    end
    eta(:,i) = sum(U(:,i1:i-1).*phi(min(t(i) - t(i1:i-1),B)),2);
end
eta(D+1,:) = ones(1,n);
logl = zeros(D,1);
for i=1:D
    if ~valid(i) % if no event happens on node i
        logl(i) = 0;
        continue;
    end
    x2 = [A(i,:)';mu(i)];
    x2 = [M(i,:)';1].*x2.*[valid;1];
    eta_i = eta(:,u==i);
    logl2 = - x2'*linear_term + sum(log(x2'*eta_i));
    logl(i) = logl2-1;
    % EM for each node
    while abs(logl(i)-logl2)>1e-3 % precision for EM
        n_iter = n_iter + 1;
        x = x2;
        logl(i) = logl2;
        %E-step
        p = eta_i.*x;
        p = p./sum(p);
        %M-step
        P = sum(p,2);
        x2 = P./linear_term;
        x2(~valid) = 0;
        logl2 = - sum(P) + sum(log(x2'*eta_i));
    end

    A(i,:) = x2(1:end-1)';
    mu(i) = x2(end);
    logl(i) = logl2;
end
end

