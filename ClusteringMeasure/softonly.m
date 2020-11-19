function E = softonly(F,lambda)
% 1 norm and F norm
%update J
%temp = F;
%[U S V] = svd(temp, 'econ');
       
diagS = diag(F);
svp = length(find(diagS > lambda));
diagS = max(0,diagS - lambda);
        
if svp < 0.5 %svp = 0
   svp = 1;
end
E = diag(diagS(1:svp)); 