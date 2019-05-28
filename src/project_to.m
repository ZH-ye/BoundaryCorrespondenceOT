function [fp,selectable]=project_to(data,x,selectable)
% project the chord length parameters to polygon vertices
if nargin <3
    selectable = get_feature(data);
end
curv_length = to_curvature_length(data);
curvature = curv_length(1,:);
I = curvature(selectable)>-1e-10;
selectable=selectable(I);


%curvature_5 = curvature_neighbor(curvature,2);
edge_length = curv_length(2,:);
%position = cumsum(edge_length);
position = cumsum(edge_length);
position = position/position(end);
ps = reshape(position(selectable),1,[]);
xl = reshape(mod(x,1),[],1);
d = abs(ps-xl);
d = min(d,1-d);
[~,I] = min(d,[],2);
x_result = reshape(ps(I),size(x));
fp = selectable(I);
fp = reshape(fp,size(x));
fp = fp';
fp = sort(fp,1);

end

function features = get_feature(data,tol)
if nargin <2
    tol=1e-4;
end
oldp = path();
addpath('../');
C = data';
C = [C;data(:,1)'];
[C_out,i_rem,CI]=DecimatePoly(C,[tol,1],false);
%[C_out,i_rem,CI]=DecimatePoly(C,[Np/1025,2],false);
%p = C_out';
N = size(data,2);
f = (1:N);
i_rem = i_rem(1:N);
f = f(i_rem == false);
%features = f;
features = farest_admissible_alt(data,f,true);
path(oldp);
end