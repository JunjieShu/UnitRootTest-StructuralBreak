




function result = Parron_Model(breaks_num, model_type, T, simu_num, lag_max)
if length(model_type) ~= breaks_num
    error('The breaks number do not match the model type.')
end
if breaks_num > T-2
    error('The number of sample is too small.')
end

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

%%% 需要记录的参数个数
para_num = 3 + sum(DUt_idc) + sum(DTt_idc);

%%% 
simu_resule = zeros(simu_num, para_num);

%%% 开始进行模拟
for simu_num_count = 1 : simu_num
    cValue_table = zeros(T, T, para_num);
    
    t = [2:T]';  % 模型中的变量t
    time_mat = repmat(t',[breaks_num,1]);
    time_mat_point = [1 : breaks_num];
    
% 以下的一个大循环先指定了time_mat_point（指针）的值，然后根据这个指针分别指定对应
% 的Tb1,Tb2 ··· Tbbreaks_num的值，然后在对应地计算出其t统计量的值。
while(1==1)
    for inv_point_point = 1 : breaks_num
        point_point = -inv_point_point + breaks_num + 1;
        time_mat_point(point_point) = time_mat_point(point_point) + 1;
        
        if time_mat_point(point_point) >= T-1 && point_point > 1
            time_mat_point(point_point-1) = time_mat_point(point_point-1) + 1;
            time_mat_point(point_point) = time_mat_point(point_point-1) + 1;
        end       
    end
    
    if time_mat(1, time_mat_point(1)) > (T-1)-breaks_num
        break
    end
    
end
    
    

    
    
    
    for Tb1 = 2 : T - 2
        for Tb2 = Tb1 + 1 : T - 1
            y = y_rnd(T);
            delta_y = diff(y);
            % 以下循环用于寻找恰当的差分滞后阶数
            for inv_p = 1 : lag_max + 1
                p = -inv_p + lag_max + 1;       % 滞后阶数
                y_L1 = y(1:end-1);              % 一阶滞后
                mu = ones(size(delta_y,1),1);
                t = [2 : T]';
                
                %%% 生成虚拟变量DUt和DTt
                %{
                DU1t = (t > Tb1)*1;
                DT1t = (t > Tb1).*(t-Tb1);
                DU2t = (t > Tb2)*1;
                DT2t = (t > Tb2).*(t-Tb2);
                %}
                
                DUt = zeros(t, breaks_num);
                DTt = zeros(t, breaks_num);
                for temp = 1 : breaks_num
                    %DUt(:, breaks_num) = 
                end
                clear temp
                
                
                
                
                if p == 0
                    lag_zhihouxiang = [];
                else
                    lag_zhihouxiang = zeros(size(delta_y,1), p);
                    for j = 1 : p
                        lag_zhihouxiang(:,j) = [zeros(j,1); delta_y(1:end-j)];
                    end
                end
                
                X = [y_L1, mu, t, DU1t, DT1t, DU2t, DT2t, lag_zhihouxiang];
                [b,~,r] = regress(delta_y,X);
                
                % 计算最高阶查分滞后项的t值
                if p ~= 0
                sigma = sqrt(sum(r.^2)/(length(delta_y)-size(X,2)));
            se = sqrt(sigma^2 * diag( (X'*X)^(-1) ) );
            t_stat = b(1:7)./se(1:7);
                end
                
                % 判断t统计量是否显著，真为显著
                if  p ==0 || abs(t_stat(end)) >= 1.6 
                    break
                end
            end
            
            % 计算当前模拟下各个参数的的t值(多元)
            if p == 0
            sigma = sqrt(sum(r.^2)/(length(delta_y)-size(X,2)));
            se = sqrt(sigma^2 * diag( (X'*X)^(-1) ) );
            t_stat = b(1:7)./se(1:7);
            end
            
            if isreal(t_stat) ~= 1
                cValue_table(Tb1, Tb2, 1) = inf;
                break
            end
            
            cValue_table(Tb1, Tb2, :) = t_stat;
        end
    end
    
    
end













end


