function result = change_series(alpha1,alpha2,rho1,rho2,sigma,T,TB)

y = zeros(T,1);
y(1) = normrnd(0,sigma);
for t = 2 : T;
    if t < TB
        y(t) = alpha1 + rho1*y(t-1) + normrnd(0,sigma);
    else
        y(t) = alpha2 + rho2*y(t-1) + normrnd(0,sigma);
    end
end

result = y;
end