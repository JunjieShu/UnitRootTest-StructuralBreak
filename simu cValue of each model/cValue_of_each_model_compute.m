clock

clear all
addpath('../my_function')



result = Parron_Model_2( 'BB', 60, 2000, 5);
xlswrite('BB',result)
clear all
clock

result = Parron_Model_2( 'AA', 60, 2000, 5);
xlswrite('AA',result)
clear all
clock



