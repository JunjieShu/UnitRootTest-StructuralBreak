function result = shizheng_function_2(data, model_type, T, lag_max)
%SHIZHENG_FUNCTION_2 这个函数用于计算时间序列在2个间断点的情况下的各个参数值，要与对应
% 的临界值表相对照来判断该序列在这个模型下的显著性。
% 注意'data'这个参数。由“时间”、“数据”、“时间”、“数据”...按列交替

breaks_num = 2;

DUt_idc = zeros(breaks_num, 1);     % A突变,DUt的指标集合
DTt_idc = zeros(breaks_num, 1);     % B突变,DTt的指标集合
for simu_num_count = 1 : breaks_num
    switch model_type(simu_num_count)
        case 'A'
            DUt_idc(simu_num_count) = 1;
        case 'B'
            DTt_idc(simu_num_count) = 1;
        otherwise
            DUt_idc(simu_num_count) = 1;
            DTt_idc(simu_num_count) = 1;
    end
end

para_num = 3 + sum(DUt_idc) + sum(DTt_idc);  %%% 参数个数，不含t
time_spot_and_t_stat_and_p = zeros(size(data,2)/2, 2 + para_num + 1); %%% 记录结果

for y_num = 1 : size(data,2)/2
    t_stat_table = zeros(T,T,para_num + 1);    
    % 用于记录不同Tb1和Tb2下的t统计量
    % 第3维顺序分别为rho,mu,beta,DUt,DTt,滞后阶数p
    
    y = data(:, 2*(y_num-1)+2);
    
    year = data(:, 2*(y_num-1)+1);
    delta_y = diff(y);
    
    for Tb1 = floor(0.15*T) : ceil(0.85*T - 2)
        for Tb2 = Tb1 + 2 : ceil(0.85*T)

            if Tb1 == 0
                Tb1 =  Tb1 + 1;
            end
            if Tb2 == T | Tb2==T-1
                break
            end
        
            % 确定滞后阶数
            for inv_p = 1 : lag_max + 1
                p = -inv_p + lag_max + 1;       % 滞后阶数
                y_L1 = y(1:end-1);              % 一阶滞后
                mu = ones(size(delta_y,1),1);
                t = [2 : T]';
                DU1t = (t > Tb1)*1;
                DT1t = (t > Tb1).*(t-Tb1);
                DU2t = (t > Tb2)*1;
                DT2t = (t > Tb2).*(t-Tb2);
                   
                %%% 生成虚拟变量DUt和DTt
                DU1t = (t > Tb1)*1;
                DU2t = (t > Tb2)*1;
                DT1t = (t > Tb1).*(t-Tb1);
                DT2t = (t > Tb2).*(t-Tb2);
                
                DUt = [DU1t, DU2t];
                DTt = [DT1t, DT2t];
                
                temp1 = breaks_num;
                while(temp1>=1)
                    if DUt_idc(temp1) == 0
                        DUt(:,temp1)=[];
                    end
                    if DTt_idc(temp1) == 0
                        DTt(:,temp1)=[];
                    end
                    temp1 = temp1 - 1;
                end
                               
                if p == 0
                    lag_zhihouxiang = [];
                else
                    lag_zhihouxiang = zeros(size(delta_y,1), p);
                    for j = 1 : p
                        lag_zhihouxiang(:,j) = [zeros(j,1); delta_y(1:end-j)];
                    end
                end
                
                X = [y_L1, mu, t, DUt, DTt, lag_zhihouxiang];
                [b,~,r] = regress(delta_y,X);
                
                % 计算最高阶差分滞后项的t值
                if p ~= 0
                    sigma = sqrt(sum(r.^2)/(length(delta_y)-size(X,2)));
                    se = sqrt(sigma^2 * diag( (X'*X)^(-1) ) );
                    t_stat = b./se;
                end
                % 判断t统计量是否显著，真为显著
                if p == 0 || abs(t_stat(end)) >= 1.6
                    break
                end
            end
            
            % 计算当前模拟下各个参数的的t值(多元)
            if p == 0
                sigma = sqrt(sum(r.^2)/(length(delta_y)-size(X,2)));
                se = sqrt(sigma^2 * diag( (X'*X)^(-1) ) );
                t_stat = b./se;
            end
            
            if isreal(t_stat) ~= 1
                t_stat_table(Tb1, Tb2, :) = nan(para_num + 1,1);
                continue
            end
            
            t_stat_table(Tb1, Tb2, :) = [t_stat(1:para_num);p];
        end
    end
    
    % 找最小的t_\rho，用于判断是不是平稳的
    [a,b]=find(t_stat_table==min(t_stat_table(find(t_stat_table(:,:,1)~=0))));
    time_spot_and_t_stat_and_p(y_num,:) = [year(a), year(b), reshape(t_stat_table(a,b,:),1,para_num+1)];
end

result = time_spot_and_t_stat_and_p;
end
