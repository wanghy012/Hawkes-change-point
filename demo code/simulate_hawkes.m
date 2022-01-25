function [t,u] = simulate_hawkes(A0,A1,mu0,mu1,tau,T,beta)
%% simulate Hawkes process with change-point.
D = size(A0,1);
if nargin<7
    beta=1;
end
%scale 
tau = tau*beta;
T = T*beta;
mu0 = mu0/beta;
mu1 = mu1/beta;
Mu0 = sum(mu0);
Mu1 = sum(mu1);
%simulate
dst = zeros(D,1);
t1 = 0;
u = [];

%%before change
while 1
    gap = get_next(Mu0,sum(A0*dst));
    if t1(end) + gap > tau
        break;
    end
    t1(end+1) = t1(end) + gap;
    dst = dst * exp(-gap);
    r = rand()*(Mu0+sum(A0*dst));
    v=1;
    s=0;
    for j=1:D-1
        s = s+mu0(j)+A0(j,:)*dst;
        if r>s
            v = v+1;
        else
            break;
        end
    end
    u(end+1) = v;
    dst(v) = dst(v) + 1;
end
dst = zeros(D,1);
t2 = tau;
%%post change
while 1
    gap = get_next(Mu1,sum(A1*dst));
    if t2(end) + gap > T
        break;
    end
    t2(end+1) = t2(end) + gap;
    dst = dst * exp(-gap);
    r = rand()*(Mu1+sum(A1*dst));
    v=1;
    s=0;
    for j=1:D-1
        s = s+mu1(j)+A1(j,:)*dst;
        if r>s
            v = v+1;
        else
            break;
        end
    end
    u(end+1) = v;
    dst(v) = dst(v) + 1;
end
t = [t1(2:end),t2(2:end)];
t = t/beta;
end

function [gap] = get_next(Mu,dst)
% r = -log(rand());
% r = min(r,11);
% a = 0;
% b = 11/Mu;
% while b-a > 1e-04
%     t = Mu*(a+b)/2 + dst - dst*exp(-(a+b)/2);
%     if t>r
%         b=(a+b)/2;
%     else
%         a = (a+b)/2;
%     end
% end
% gap = a;
gap = 0;
r = -log(rand())/(Mu+dst);
gap = gap +r;
while rand()>(Mu+dst*exp(-r))/(Mu+dst)
    dst = dst * exp(-r);
    r = -log(rand())/(Mu+dst);
    gap = gap +r;
end
end

