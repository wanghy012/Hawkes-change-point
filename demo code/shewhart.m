function [statslog] = shewhart(t,T,w)
% t the event times
% T the time horizon
% w the window size. should be approximately equal/larger than EDD
% b the lower threshold and upper threshold
stop_time = -1;
tau_hat = -1;
i1 = 1;
count =0;
i2 = 1;
l = length(t);
statslog = zeros(2,2*l);
c = 1;
t = [t Inf];
while t(i2)<w
    i2 = i2+1;
end
statslog(:,1) = [w;i2-i1];
c = c+1;
while i2 <= l
    while t(i1)<t(i2)-w
        i1 = i1+1;
        statslog(:,c) = [t(i1-1)+w;i2-i1];
        c = c+1;
    end
    i2 = i2+1;
    statslog(:,c) = [t(i2-1);i2-i1];
    c = c+1;
end
while t(i1)<T-w
    i1 = i1+1;
    statslog(:,c) = [t(i1-1)+w;i2-i1];
    c = c+1;
end
statslog(:,c:end) = [];

end

