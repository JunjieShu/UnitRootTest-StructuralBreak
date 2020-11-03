%% temp!!

function result = y_rnd(T)
y = zeros(T,1);
for t = 2 : T
    y(t) = y(t-1) + normrnd(0,1);
end
result = y;
end