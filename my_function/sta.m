function result = sta(alpha1,alpha2,beta1,beta2,sigma,T,TB)
t = 1 : T;
result = (alpha1+beta1*t).*(t<TB) + (alpha2+beta2*t).*(t>=TB) +...
    normrnd(0,sigma,1,1000);
end