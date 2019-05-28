function [x_real,x_imag] = matrix_completion(X,SparseAreal_index,SparseAimag_index,Sparseval_real,Sparseval_imag,Bv_real,Bv_imag,alternativeweight,nuclearnormweight)
% addpath(genpath(pwd));
%Global constants and defaults
admmpara = 0.0;
completionweight = 5000;
% nuclearnormweight = 0.0;
sz=size(X);
n = prod(sz);
MAX_ITER = 1;

%这里将隐式曲面重构左端矩阵做了行列排列变换
permutation_index = zeros(1,n);
sparseA_real = sparse(SparseAreal_index(:,1),SparseAreal_index(:,2),Sparseval_real,n,n);
sparseA_imag = sparse(SparseAimag_index(:,1),SparseAimag_index(:,2),Sparseval_imag,n,n);
for i=1:1:sz(1)
    for j=1:1:sz(2)
        permutation_index((j-1)*sz(1)+i) = (i-1)*sz(2)+j;
    end
end

%重构左端稀疏矩阵和右端项
sparseA_real = sparseA_real(permutation_index,permutation_index);
sparseA_imag = sparseA_imag(permutation_index,permutation_index);
sparseA = sparseA_real - sqrt(-1)*sparseA_imag;

%观测元素在原矩阵中的位置
Ind = zeros(2*(sz(2)-2)+2*sz(1),1);
for i=1:sz(1)
   Ind(i) = i;
end
for j=2:sz(2)-1
    Ind(sz(1)+2*(j-2)+1) = sz(1)*(j-1)+1;
    Ind(sz(1)+2*(j-1)) = sz(1)*j;
end
for i=1:sz(1)
    Ind(sz(1)+2*(sz(2)-2)+i) = sz(1)*(sz(2)-1)+i;
end

Bv = zeros(n,1);
Bv(Ind) = Bv_real + sqrt(-1)*Bv_imag;
Omega = sparse(Ind,Ind,1,n,n);
sparse_AA = 2*alternativeweight*sparseA + 2*completionweight*Omega + admmpara*speye(n);
sparse_b = 2*completionweight*Bv;

%ADMM
%initializing
Z=zeros(sz);
U=zeros(sz);
for k=1:1:MAX_ITER
    
    %step1: 求x_k+1
    X_temp = admmpara*(Z-U);
    x = sparse_AA\(sparse_b + X_temp(:));
    X = reshape(x,sz);
    
    %step2: 求Z_k+1
%     [Z,~] = softth(X+U,nuclearnormweight/admmpara); %Y其实就是\alpha
    
    %step3: 求U_k+1
    V = X - Z;
    U = U + V;
end

x_real = real(X(:));
x_imag = imag(X(:));
end