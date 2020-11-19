clear;clc;
addpath('./ClusteringMeasure');
load('./data/bbcsport_2view.mat');
fprintf('CBFMSC on BBCSport dataset\n');
gt = truth;
numC = size(unique(gt),1);
X = cell(1,2);
X{1} = data{1}';
X{2} = data{2}';

for i = 1:2
    X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1);
end

NMI_all = [];
ACC_all = [];
F_all = [];
AVG_all = [];
P_all = [];
RI_all = [];

opts.rho = 1.9;
opts.lambda_1 = 10;
opts.dim_V = 18;

fprintf('lambda_1 = %f, dim_V = %d\n', opts.lambda_1,opts.dim_V);
for i = 1:30
    fprintf('Clustering in %d-th iteration\n',i);
    [U,V]  = CBFMSC(X,opts);
    [NMI,ACC,F,AVG,P,RI] = ClusteringResults(U,V,gt,numC);
    fprintf('\tNMI: %f, ACC: %f, F: %f, AVG: %f, P: %f, RI: %f\n',NMI,ACC,F,AVG,P,RI);
    NMI_all = [NMI_all, NMI];
    ACC_all = [ACC_all, ACC];
    F_all = [F_all, F];
    AVG_all = [AVG_all, AVG];
    P_all = [P_all, P];
    RI_all = [RI_all, RI];
end
fprintf('---------------Average Results--------------\n');
fprintf('NMI: %f(%f), ACC: %f(%f), F: %f(%f), AVG: %f(%f), P: %f(%f), RI: %f(%f)\n', mean(NMI_all),std(NMI_all),mean(ACC_all),std(ACC_all),mean(F_all),std(F_all),mean(AVG_all),std(AVG_all),mean(P_all),std(P_all),mean(RI_all),std(RI_all));
