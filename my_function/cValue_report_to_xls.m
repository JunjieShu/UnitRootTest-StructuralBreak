function cValue_report_to_xls(model_type, data)
%cValue_report_to_xls 不会输出结果，只能在当前文件夹生成所选模型各个参数的分位数表(csv)
% 分位数所生成的分位数为：0.01 0.025 0.05 0.1 0.5 0.9 0.95 0.975 0.99下的分位数


alpha = [0.01 0.025 0.05 0.1 0.5 0.9 0.95 0.975 0.99];
cValue_table = zeros(size(data,2),length(alpha));
sorted_cValue = sort(data);
simu_num = size(data,1);    % 模拟次数
for i = 1 : size(data,2)
    for j = 1 : length(alpha)
        cValue_table(i,j) = sorted_cValue(round(simu_num*alpha(j)),i);
    end
end

cValue_table = [alpha; cValue_table];    % 生成表头
diyiliekongzhu = zeros(size(cValue_table,1),1);
cValue_table = [diyiliekongzhu,cValue_table];


xlswrite([model_type,'_cVale_table.xls'],cValue_table)


end