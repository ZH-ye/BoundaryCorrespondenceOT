function [fp,features,E]=feature_and_filter(data,get_features,energy,batch,topN,getpermu)
% get feature points, find preferable permutations and sort by energy
if nargin < 4
    batch = 1;
end
if nargin <5
    topN = 10;
end
N = size(data,2);
curv_length = to_curvature_length(data);
curvature = curv_length(1,:);
%curvature_5 = curvature_neighbor(curvature,2);
edge_length = curv_length(2,:);
position = cumsum(edge_length);
position = position/position(end);
features = get_features(data);
if nargin < 6
%permu = nchoosek(features,4);
    getpermu = @(features,N)crossing_pairs(features,N);   
end

%permu = getpermu(get_features(data),N);
permu_pos = getpermu(position(features),1);
permu=zeros(size(permu_pos));
for ii =features
    permu(permu_pos == position(ii))=ii;
end
n = size(permu,1);

if strcmp(batch,'all')
    batch = n;
end
if 0<topN && topN<1
    topN = int32(n*topN);
end

if batch == 1
    EE = zeros(n,1);
for ii = 1:n
    if mod(ii,50)==0
        fprintf('%d out of %d\n',ii,n);
    end
    EE(ii) = energy(to_para(permu(ii,:),N)) ;
    
end

end
if batch >1
    topN = min(topN,n);
    EE=ones(batch+topN,1)*Inf;
    IE = zeros(batch+topN,1);
    start=1;
    endding =min(start+batch-1,n);
    while(start <= endding)
        fprintf('%d - %d \n',start,endding);
        x = to_para(permu(start:endding,:),N);
        EE(topN+1:endding-start+1+topN) = energy(x);
        IE(topN+1:endding-start+1+topN) = start:endding;
        [EE,I] = sort(EE);
        IE = IE(I);
        start = endding +1;
        endding = min(start+batch-1,n);
    end
    fprintf('done\n');
    E = EE(1:topN);
    fp = permu(IE(1:topN),:)';
    return
end
% return top 10 
[E,I]=sort(EE);
E = E(1:topN);
fp = permu(I(1:topN),:)';
end

function x=to_para(fp,N)
x = (double(fp)-1)/N;
end
function permu = crossing_pairs(points,N)
[pairs,d] = get_pairs(points,N,0.25);
pairs = sort(pairs,2); % just in case ...
n = size(pairs,1);
pp = nchoosek(1:n,2);
% check crossing
p1 = pairs(pp(:,1),:);% x1< x2
p2 = pairs(pp(:,2),:);% x3< x4
x1 = p1(:,1);
x2 = p1(:,2);
x3 = p2(:,1);
x4 = p2(:,2);
x3in = (x3-x1).*(x2-x3);
x4in = (x4-x1).*(x2-x4);

crossing = x3in.*(-x4in)>0;
p1 = p1(crossing,:);
p2 = p2(crossing,:);
permu = sort([p1,p2],2);



end

function [pairs,d] = get_pairs(points,N,a)
pos = double(sort(points)-1)/N;
p = nchoosek(pos,2);
p0 = nchoosek(sort(points),2);
d1 = p(:,2)-p(:,1);
d = min(d1,1-d1);
l = d>a;
pairs = p0(l,:);
d = d(l);

end