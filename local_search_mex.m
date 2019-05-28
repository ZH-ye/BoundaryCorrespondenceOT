function [fp,E,xx] = local_search_mex(data,fp,cost_type,lambda,selectable)
% call the mex function to do the local search
N = size(data,2);
if nargin<4
    lambda= [0.1,0.1];
end
if nargin<5
    selectable=1:N;
end
fp = sort(fp,1);
curv_length = to_curvature_length(data);
curvature = curv_length(1,:);
edge_length = curv_length(2,:);
position = cumsum(edge_length);
position = position/position(end);

if sum(curvature)<1e-5
    % too noisy!
    fprintf('the input is not a boundary of a domain!\n');
    E=[];
else
t1 =tic;

min_iter=20;
max_iter=1000;

show_log = false;
[xx,E] = mexSinkhorn(curvature,position,position(fp'),lambda,cost_type,min_iter,max_iter,show_log);

toc(t1);
[E,I]= sort(E);
xx=xx(I,:);
xx = sort(mod(xx,1.0),2); % chord length parameters

%%%%%%%% post process %%%%%%%%%%%
if isempty(selectable)
    [fp,~] = project_to(data,xx);
else
[fp,~] = project_to(data,xx,selectable);
end
%[fp,~] = project_to(data,xx);
% removed invalid ones
valid = all(diff(fp,1)>0,1);
fp = fp(:,valid);
E = E(valid);
% unique column
[fp,ia,~] = unique(fp','stable','rows');
fp = fp';
E = E(ia)
% fp = 1+mod(fp+199-1,N);
fp = sort(fp,1);
end
end
