function [x_real,x_imag] = quadratic_constrained(X,SparseA_index,Sparseval,Sparseb)
sz = size(X);
sparseA = sparse(SparseA_index(:,1),SparseA_index(:,2),Sparseval,sz(1),sz(2));
opts = optimset('Algorithm','interior-point-convex','Display','off');
x=quadprog(sparseA+sparseA',Sparseb,[],[],[],[],-0.6*ones(sz(1),1),0.6*ones(sz(1),1),[],opts);
% x=quadprog(sparseA+sparseA',Sparseb,[],[],[],[],[],[],[],opts);
% Sparseb = 0.5*Sparseb;
% x = -sparseA\Sparseb;
x_real = x(1:sz(1)/2);
x_imag = x(sz(1)/2+1:end);
end
