
%%% 开始进行模拟
tic
cValue_table = zeros(T, T, T, para_num);
simu_result = zeros(1, para_num)

w = waitbar(0, 'Please wait...');
for Tb1 = 2 : T - breaks_num
    for Tb2 = Tb1 + 2 : T - breaks_num + 1
        for Tb3 = Tb2 + 2 : T - breaks_num + 2
            
            y = y_rnd(T);
            delta_y = diff(y);
            
            % 以下循环用于寻找恰当的差分滞后阶数
            for inv_p = 1 : lag_max + 1
                p = -inv_p + lag_max + 1;       % 滞后阶数
                y_L1 = y(1:end-1);              % 一阶滞后
                mu = ones(size(delta_y,1),1);
                t = [2 : T]';
                
                %%% 生成虚拟变量DUt和DTt
                DU1t = (t > Tb1)*1;
                DU2t = (t > Tb2)*1;
                DU3t = (t > Tb3)*1;
                DT1t = (t > Tb1).*(t-Tb1);
                DT2t = (t > Tb2).*(t-Tb2);
                DT3t = (t > Tb3).*(t-Tb3);
                
                DUt = [DU1t, DU2t, DU3t];
                DTt = [DT1t, DT2t, DT3t];
                
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
                
                % 计算最高阶查分滞后项的t值
                if p ~= 0
                    sigma = sqrt(sum(r.^2)/(length(delta_y)-size(X,2)));
                    se = sqrt(sigma^2 * diag( (X'*X)^(-1) ) );
                    t_stat = b./se;
                end
                
                % 判断t统计量是否显著，真为显著
                if  p == 0 || abs(t_stat(end)) >= 1.6
                    break
                end
            end
            
            % 计算当前模拟下各个需要的参数的t值(多元)
            if p == 0
                sigma = sqrt(sum(r.^2)/(length(delta_y)-size(X,2)));
                se = sqrt(sigma^2 * diag( (X'*X)^(-1) ) );
                t_stat = b(1:para_num)./se(1:para_num);
            end
            
            if isreal(t_stat) ~= 1
                cValue_table(Tb1, Tb2, 1) = nan;
                continue
            end
            
            cValue_table(Tb1, Tb2, Tb3, :) = t_stat(1:para_num);
        end
    end
    
    
    waitbar(Tb1 /  (T - breaks_num-1), w)
end

[min_val, position] =min(cValue_table(find(cValue_table(:,:,:,1)~=0)));
[i,j,k,l] = ind2sub(size(cValue_table),position);
simu_result(simu_num_count,:) = cValue_table(i, j, k, : );






close(w)
toc