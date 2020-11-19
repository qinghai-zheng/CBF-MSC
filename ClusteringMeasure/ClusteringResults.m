function [NMI,ACC,F,AVG,P,RI] = ClusteringResults(U,V,gt,numC)
%ClusteringResults

% Syntax: [NMI,ACC,F,AVG,P,RI] = ClusteringResults(input)
%
% Copyright: zhengqinghai@stu.xjtu.edu.cn
% 2019/03/07

U_tmp = zeros(size(U{1}));
v = size(U,2);
for i = 1:v
    U_tmp = U_tmp+U{i};
end
U_tmp = U_tmp/v;
Z_tmp = U_tmp*V;
[NMI,ACC,F,AVG,P,RI]=clustering(abs(Z_tmp)+abs(Z_tmp'), numC, gt);