%% 生成各个模型下的分位数表

% CC 模型用函数生成的分位数表效果并不好，用以前的

model = 'BB';
data = csvread([model,'.csv']);

alpha = [0.01 0.025 0.05 0.1 0.5 0.9 0.95 0.975 0.99];
cValue_table = zeros(size(data,2),length(alpha));
sorted_cc_cValue = sort(data);
simu_num = size(data,1);    % 模拟次数
for i = 1 : size(data,2)
    for j = 1 : length(alpha)
        cValue_table(i,j) = sorted_cc_cValue(simu_num*alpha(j),i);
    end
end

xlswrite(cValue_table, 'BB_cVale_table')