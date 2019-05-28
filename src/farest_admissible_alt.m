function chosen = farest_admissible_alt(data,fp,add_prominent,Np)
% simplify the polygon, so that about 20 points are left
if nargin <3
    add_prominent = true;
end
if nargin <4
    Np = 4;
end
N = size(data,2);
curv_length = to_curvature_length(data);
curvature = curv_length(1,:);
curvature_5 = curvature_neighbor(curvature,2);
edge_length = curv_length(2,:);
position = cumsum(edge_length);

admissible = curvature>-1e-10;
levels = [1e-10,pi/2*[0.01,0.1,0.5,0.7,1]];
for level = levels
    prominent = curvature>level;
    if nnz(prominent)<0.01*N
        break;
    end
end
fp = fp(fp>0);
chosen = fp(admissible(fp));

if size(chosen,1)>size(chosen,2)
    chosen = chosen';
end
if add_prominent
% pp = 1:N;
% pp = pp(curvature_5>pi/2*0.5);
% chosen = [chosen,pp];
[~,II] = sort(-curvature);
chosen = [chosen,II(1:Np)];
end
if isempty(chosen)
[~,I] = max(curvature);
chosen=[I];
end
chosen = unique(sort(chosen));
distences = get_dist(chosen,position);
for ii = 1:(12-numel(chosen))
    %[chosen,distences] = farest_admissible_helper(chosen,distences,position,prominent);
end
if Np>0
    for ii = 1:Np
        [chosen,distences] = farest_admissible_helper(chosen,distences,position,admissible);
    end
else
    for ii = 1:min(4,16-numel(chosen))
        [chosen,distences] = farest_admissible_helper(chosen,distences,position,admissible);
        
    end
end
chosen = unique(sort(chosen));
end

function dist = get_dist(chosen, position)
N=numel(position);
total_length = position(end);
position = [0,position(1:end-1)];
dist = 10000000* ones(size(position));
n = numel(chosen);
for ii = 1:n
    p_new = chosen(ii);
        dist_one_side = abs(position-position(p_new));
    
    dist_new = min(dist_one_side,...
        total_length-dist_one_side);
    dist = min(dist,dist_new);
end
end

function [chosen,dist] = farest_admissible_helper(chosen,dist,position,admissible)
    % one point is newly taken, cauculate the distences
    N=numel(admissible);
    p_new = chosen(end);
    total_length = position(end);
    position = [0,position(1:end-1)];
    dist_one_side = abs(position-position(p_new));
    
    dist_new = min(dist_one_side,...
        total_length-dist_one_side);
    dist = min(dist,dist_new);
    % find the index of farest admissible point
    dist_index = [dist;(1:N)];
    adm = dist_index(:,admissible);
    [~,I] = max(adm(1,:));
    idx = adm(2,I);
    chosen = [chosen,idx];
end

function c = curvature_neighbor(c,n)
N = numel(c);
idx = 1:N;
cc=zeros(size(c));
for ii = -n:n
   cc = cc+ cc(1+mod(idx-ii-1,N));
end
end
