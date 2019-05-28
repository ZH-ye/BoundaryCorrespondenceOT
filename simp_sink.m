function  [fp,selectable]=simp_sink(data,p,lambda,use_local_search,selectable_str)
% simplify the polygon to find initial values for corner selection
% then run the local search step
% data: 2*N , the polygon
% p : p,[tau]
% lambda: 2
% use_local_search: bool
% selectable_str: char array, influence the post process.
%                possible values : 'all'  / 'positive curvature' / 
%                                   'feature points'
%                
if nargin <3
    lambda = [0.1,0.1];
end
if nargin < 4
    use_local_search = true;
end
if nargin < 5
selectable_str='all'; 
end
global g_time;
g_time=[];
g_time.time = zeros(1,2);
g_time.size = zeros(1,2);
t1 = tic;
curv_length = to_curvature_length(data);
curvature = curv_length(1,:);
edge_length = curv_length(2,:);
position = cumsum(edge_length);
position = position/position(end);
mu = curvature;
mu = mu ./ sum(mu);
g_time.size(1) = nnz(mu>1e-10);
g_time.size(2) = 4+nnz(mu<-1e-10);

 energy = @(x) object(x,mu,p,lambda);

%%%%%%%%%%%%%%%%% get initial values of corners %%%%%%%%%%%%
 [fp,selectable,E0] = feature_and_filter(data,@(d)get_feature(d),energy,1024,10);
 E0
g_time.time(1) = toc(t1);
t1 = tic;
if ~ use_local_search
    return;
end

%%%%%%%%%%%%%%%%% further optimization %%%%%%%%%%%%%%%


%%%%%%%% a switch for post process %%%%%%%
N = size(data,2);
switch selectable_str(1)
    case 'a' % all points can be selected as corners
        sss = 1:N;
    case 'p' % only points with positive curvatures
        sss = find(curvature>0);
    case 'f' % only feature points (determined by polygon decimate altorithm)
        sss = [];
end
% local search start with initial values in fp 
% start with 10 different initial values
[fp,E]=local_search_mex(data,fp,p,lambda,sss); 
g_time.time(2)= toc(t1);
end


function features = get_feature(data,Np)
if nargin <2
    Np=16;
end
N = size(data,2);
% oldp = path();
% addpath('../');
C = data';
C = [C;data(:,1)'];

[C_out,i_rem,CI]=DecimatePoly(C,[1e-4,1],false);
if numel(C_out)/2 -1 >Np
[C_out,i_rem,CI]=DecimatePoly(C,[Np/(N+1),2],false);
end

f = (1:N);
i_rem = i_rem(1:N);
f = f(i_rem == false);
%features = f;
 features = farest_admissible_alt(data,f,true);
% path(oldp);
end

function x = to_length_para(x,position)
N = numel(position);
x_pos = int32(1+x*(N-1));
x = position(x_pos);
end
function d = get_sink(cost_type,x,c)

% x: n by 4
N = numel(c);
%[~,d]=get_phi([0.25,0.25,0.25,0.25],x,c,1/N*ones(N,1));
% call sinkhorn dual distense
n = size(x,1);
x_pos = int32(1+x*(N-1));
x_pos_linear = x_pos + int32((0:N:(n-1)*N)');
a = zeros(N,1);
b = zeros(N,n);
a(c>0) = c(c>0);
c_neg = -c(c<0);
if size(c_neg,2)>size(c_neg,1)
    c_neg = c_neg'; % to column vector
end
c_neg = repmat(c_neg,1,n);
b(c<0,:) = c_neg;
b(x_pos_linear) = b(x_pos_linear)+0.25;
a = reshape(a,[N,1]);
a = a / sum(a);
%b = reshape(b,[N,n]);
b = b./sum(b,1);
xs = (0:N-1)/N;
M = mod( xs-xs',1);
M = min(M,1-M);
if numel(cost_type)>=2
    thr = cost_type(2);
    if thr >0
        M = min(M,thr);
    end
end
if cost_type(1) ~=1
M = (M+1e-10).^cost_type(1);
end
l = 100/(max(M(:))); % l*max(M)<200 as suggested by http://marcocuturi.net/SI.html
K = exp(-l*M);
% U = K.*M;
a_posi = a>1e-10;
b_posi = any(b>1e-10,2);
ap  = a(a_posi);
bp = b(b_posi,:);
Mp = M(a_posi,b_posi);
% Kp = K(a_posi,b_posi);
lp = 100/(max(Mp(:))); 
Kp = exp(-lp*Mp);
Up = Kp.*Mp;
size(Mp)

t1 = tic;
[D,L]=sinkhornTransport(ap,...
    bp,Kp,Up,lp,[],[],5e-3);
d = D';

end
function d = object(x,c,p,lambda)
emd = get_sink(p,x,c);
if size(c,1)==2
    x = to_length_para(x,c(2,:));
end
E_op=0;
E_pair=0;
if lambda(1)>1e-10
    E_op = get_equal_opposite_side(x);
end

if lambda(2)>1e-10
    E_pair = get_equal_pair(x);
end

d = emd + lambda(1)*E_op +  lambda(2)*E_pair;

end
