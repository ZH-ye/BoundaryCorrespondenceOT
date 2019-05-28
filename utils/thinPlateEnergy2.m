function H = thinPlateEnergy2(knotsx,knotsy,mx,my)
orderx = numel(knotsx) - mx;
ordery = numel(knotsy) - my;

[gx,whx] = getGaussQuadrature(knotsx,orderx);
[gy,why] = getGaussQuadrature(knotsy,ordery);
% x
vx = spcol(knotsx,orderx,brk2knt(gx,3));
vy = spcol(knotsy,ordery,brk2knt(gy,3));

[X,Y] = meshgrid(1:mx,1:my);
X = X';
Y = Y';
[I,J] = meshgrid(1:mx*my,1:mx*my);
L = (knotsx(X(I(:))+orderx)> knotsx(X(J(:))) | knotsx(X(J(:))+orderx)> knotsx(X(I(:))))'  ...
    & (knotsy(Y(I(:))+ordery)> knotsy(Y(J(:))) | knotsy(Y(J(:))+ordery)> knotsy(Y(I(:))))';% ...
    %& I(:) >= J(:);
I1 = I(L);
J1 = J(L);
v=zeros(size(I1));
[intbx0,intbx1,intbx2] = intb(vx,gx,whx,mx,orderx); 
[intby0,intby1,intby2] = intb(vy,gy,why,my,ordery); 
for it=1:length(I1)
    i = X(I1(it));
    j = Y(I1(it));
    k = X(J1(it));
    l = Y(J1(it));
    %swap
    if i>k
        temp = i;
        i = k;
        k = temp;
%        [k,i] = deal(i,k);
    end
    if j>l
        temp = j;
        j = l ; 
        l = temp;
%        [l,j] = deal(j,l);
    end
    
    % matrix element for <B^{ij}, B^{kl}>
%     vvx = vx(:,i).*vx(:,k);
%     intbik0 = whx* vvx(1:3:end);% $\int B^i_{xx} B^k_{xx}$
%     intbik1 = whx* vvx(2:3:end);% $\int B^i_{x} B^k_{x}$
%     intbik2 = whx* vvx(3:3:end);% $\int B^i B^k$
%     vvy = vy(:,j).*vy(:,l);
%     intbjl0 = why*vvy(1:3:end);
%     intbjl1 = why*vvy(2:3:end);
%     intbjl2 = why*vvy(3:3:end);
    intbik0 = intbx0(i,k);% $\int B^i_{xx} B^k_{xx}$
    intbik1 = intbx1(i,k);% $\int B^i_{x} B^k_{x}$
    intbik2 = intbx2(i,k);% $\int B^i B^k$
    intbjl0 = intby0(j,l);
    intbjl1 = intby1(j,l);
    intbjl2 = intby2(j,l);
    v(it) = intbik2*intbjl0 + 2*intbik1*intbjl1 ...
        + intbik0*intbjl2;
end
H=sparse(I1,J1,v);

end

function [intb0,intb1,intb2] = intb(v,x,wh,m,order)
intb0=zeros(m);
intb1=zeros(m);
intb2=zeros(m);
for i = 1:m
    for j = i:min(m,i+order-1)
        vv=v(:,i).*v(:,j);
        intb0(i,j)=wh*vv(1:3:end);
        intb1(i,j)=wh*vv(2:3:end);
        intb2(i,j)=wh*vv(3:3:end);
    end
end
end

function [x,wh]= getGaussQuadrature(knots,order)
% get gauss quadrature for spline f with (knots) of (order)
% so that \int_0^1 f = wh*f(x)
[gx,gw] = thirdParty.lgwt(order,-1,1);
uknots = unique(knots);
h = uknots(2:end)-uknots(1:end-1);
x(order*length(h))=0;
wh(order*length(h))=0;
for i = 1:length(h)
    x(1+(i-1)*order:i*order) = (gx+1)*(h(i)/2)+uknots(i);
    wh(1+(i-1)*order:i*order) = gw*h(i)/2;
end
[x,I]=sort(x);
wh = wh(I);
end