function [U,V] = CBFMSC(X,opts)
%CBFMSC - Constrained Bilinear Factorization Multi-view Subspace Clustering
%
% Syntax: [U,V]  = CBFMSC(X,opts)
%
% Copyright: zhengqinghai@stu.xjtu.edu.cn
% 2019/06/13

%% Initialization variables
v = size(X,2); % number of views
n = size(X{1},2); % number of samples

lambda_1 = opts.lambda_1;
dim_V = opts.dim_V;
rho = opts.rho;

mu = 1e-4;
max_mu = 1e6;

E = cell(1,v);
U = cell(1,v);
Z = cell(1,v);
V = zeros(dim_V,n);
L = zeros(dim_V,n);
Y_1 = cell(1,v);
Y_2 = cell(1,v);
Y_3 = zeros(dim_V,n);

dim_all_view = zeros(v,1);
for i = 1:v
    dim_all_view(i) = size(X{i},1);
    E{i} = zeros(dim_all_view(i),n);
    U{i} = zeros(n,dim_V);
    Z{i} = zeros(n,n);
    Y_1{i} = zeros(dim_all_view(i),n);
    Y_2{i} = zeros(n,n);
    Z{i} = rand(n,n);
end

%% Optimization
iter_max = 100;
iter_curr = 0;
conv_flag = 0;
conv_threshold = 1e-6;

while conv_flag==0 && iter_curr<iter_max

    % updating V:
    V_a = mu*eye(dim_V);
    V_b = mu*L-Y_3;
    for i = 1:v
        V_a = V_a + mu*U{i}'*U{i};
        V_b = V_b + U{i}'*Y_2{i} + mu*U{i}'*Z{i};
    end
    V = V_a\V_b;

    % updating L:
    L_a = V + Y_3/mu;
    L = softth(L_a,lambda_1/mu);

    % updating variables that need to be updated for each views
    for i = 1:v
        % updating E
        E_a = X{i} - X{i}*Z{i} + Y_1{i}/mu;
        E{i} = solve_l1l2(E_a,1/mu);

        % updating U:
        U_a = Z{i} + Y_2{i}/mu;
        U_b = V*U_a';
        [svd_U,~,svd_V] = svd(U_b,'econ');
        U{i} = svd_V*svd_U';

        % updating Z:
        Z_a = mu*X{i}'*X{i} + mu*eye(n);
        Z_b = X{i}'*Y_1{i} + mu*(X{i}'*X{i}-X{i}'*E{i}) -Y_2{i} + mu*(U{i}*V);
        Z{i} = Z_a\Z_b;

        % updating Y_1 and Y_2:
        Y_1{i} = Y_1{i} + mu*(X{i} - X{i}*Z{i} - E{i});
        Y_2{i} = Y_2{i} + mu*(Z{i} - U{i}*V);
    end

    % updating Y_3:
    Y_3 = Y_3 + mu*(V - L);

    % updating mu:
    mu = min(rho*mu,max_mu);

    % check the convergence conditions
    conv_count_tmp = 0;
    for i = 1:v
        if norm(X{i}-X{i}*Z{i}-E{i},inf) < conv_threshold
            conv_count_tmp = conv_count_tmp + 1;
        end
        if norm(Z{i}-U{i}*V,inf) < conv_threshold
            conv_count_tmp = conv_count_tmp + 1;
        end
        if norm(V-L,inf) < conv_threshold
            conv_count_tmp = conv_count_tmp + 1;
        end
    end

    if conv_count_tmp == 3*v
        conv_flag = 1;
        fprintf('Convergence conditions achieved\n');
    end

    iter_curr = iter_curr + 1;
end