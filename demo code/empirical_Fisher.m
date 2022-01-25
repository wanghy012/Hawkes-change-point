function Ihat = empirical_Fisher(t,u,beta,T,A0,mu0)
% get the empirical Fisher Information on unit time.
N = length(t);
D = length(mu0);
eta = zeros(D,1);
Ihat = cell(D,1);
for i=1:D
    Ihat{i} = zeros(D+1,D+1);
end
for i=1:N-1
    %eta(:,i+1) = eta(:,i)*exp(beta*(t(i)-t(i+1)));
    %eta(u(i),i+1) = eta(u(i),i+1) + beta*exp(beta*(t(i)-t(i+1)));
    eta(u(i)) = eta(u(i))+beta;
    eta = eta*exp(beta*(t(i) - t(i+1)));
    lambda = A0(u(i+1),:)*eta + mu0(u(i+1));
    Ihat{u(i+1)} = Ihat{u(i+1)} + [eta;1]*[eta',1]/lambda/lambda;
end
for i=1:D
    Ihat{i} = Ihat{i}/T;
end
end

