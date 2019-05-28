function coeff = fit_boundary(mx,my,knots_x,knots_y,points,fp)
% points: 2 by M
% fp: index of corners given in order: left bottom, right bottom,
%                                      top right  , top left

if size(points,1)>size(points,2)
    points=points';
end
cornerpoints= points(:,sort(fp))';
cc = cornerpoints;
fixed = [-1,-1;1,-1;1,1;-1,1];
dist_Min= 10000;
for i = 1:4
    c = [cornerpoints(i+1:end,:);cornerpoints(1:i,:)];
    dd = fixed-c;
    dist = sum(dd(:).^2);
    if dist <dist_Min
        dist_Min = dist;
        cc = c;
    end
end
cornerpoints=cc;
curv_length = to_curvature_length(cc');
outerangle = curv_length(1,:);


coeff = zeros(size(points,1),mx,my);
num_pt = size(points,2);
offset = fp(1)-1;
fp = mod(fp-offset-1,num_pt)+1;
fp(5)=num_pt+1;
points = [points(:,(offset+1):end),points(:,1:(offset))];
% indicate the curve that each point should be assigned to
cp.cpts=[];
cp.curvenum=zeros(1,num_pt);
for curveN= 1:4
    cp.curvenum(fp(curveN):(fp(curveN+1)-1))=curveN;
end
cp.parameters=zeros(1,num_pt);
%1
n = get_outer_normal(points,fp(1:2));
[cc,p] = curve_fit(knots_x,mx,points(:,fp(1):fp(2)),n);
coeff(:,:,1)=cc;
cp.parameters(fp(1):fp(2)-1)=p(1:end-1);
%2
n = get_outer_normal(points,fp([2,3]));
[cc,p] = curve_fit(knots_y,my,points(:,fp(2):fp(3)),n);
coeff(:,end,:)=cc;
cp.parameters(fp(2):fp(3)-1)=p(1:end-1);
%3
n = get_outer_normal(points,fp([4,3]));
[cc,p] = curve_fit(knots_x,mx,fliplr(points(:,fp(3):fp(4))),n);
coeff(:,:,end)=cc;
cp.parameters(fp(3):fp(4)-1)=fliplr(p(2:end));
%4
n = get_outer_normal(points,fp([1,4]));
[cc,p] = curve_fit(knots_y,my,fliplr([points(:,fp(4):end),points(:,1)]),n);
coeff(:,1,:)=cc;
cp.parameters(fp(4):fp(5)-1)=fliplr(p(2:end));
Hpe = thinPlateEnergy2(knots_x,knots_y,mx,my);
mu = 0.0001;
%[X,Y] = meshgrid(1:mx,1:my);
% [coeffx,coeffy ]= meshgrid((0:mx-1)/(mx-1),(0:my-1)/(my-1));
% coeffx=coeffx';coeffy=coeffy';
% coeff = reshape([coeffx(:)';coeffy(:)'],[2,size(coeffx)]);%+0.1*rand([2,size(coeffx)]);
% u = coeff(1,:,:);
% v = coeff(2,:,:);
% x = (1-u).*(1-v)*points(1,fp(1))+ ...
%     u.*(1-v)* points(1,fp(2))+...
%     u.*v*points(1,fp(3))+...
%     (1-u).*v * points(1,fp(4));
% y = (1-u).*(1-v)*points(2,fp(1))+ ...
%     u.*(1-v)* points(2,fp(2))+...
%     u.*v*points(2,fp(3))+...
%     (1-u).*v * points(2,fp(4));
% coeff = zeros(2,mx,my);
% coeff(1,:,:)=x;
% coeff(2,:,:)=y;
cb.mx=mx;
cb.my=my;
cb.knotsx=knots_x;
cb.knotsy=knots_y;
opts = optimoptions('quadprog','Display','off');
global use_trans_b;
if ~isempty(use_trans_b)
    if use_trans_b
        if any(outerangle<pi/3)
            tform = fitgeotrans(fixed,fixed,'affine');
        else
            tform = fitgeotrans(cornerpoints,fixed,'affine');
        end
    else
        tform = fitgeotrans(fixed,fixed,'affine');
        
    end
else
    if any(outerangle<pi/3)
        tform = fitgeotrans(fixed,fixed,'affine');
    else
        tform = fitgeotrans(cornerpoints,fixed,'affine');
    end
end

% tform = fitgeotrans(fixed,fixed,'affine');
%points = tform.transformPointsForward(points')';
coeff_origin = coeff;
[x,y] = tform.transformPointsForward(coeff(1,:,:),coeff(2,:,:));
coeff(1,:,:) =x;
coeff(2,:,:) =y;
for k = 1:1
    %cp = closestp(points,cb,coeff,cp,true);
    %[H,g,C]=GeoDist(points,cb,coeff,cp);
    [cc,boundary]= tosolution(coeff);
    beq = cc(boundary>0);
    remove_zero_rows=@(a)a(any(a,2),:);
    Aeq = remove_zero_rows(diag(boundary));
    %dist1 = sum(sum((cp.cpts-pts).^2));
    %cc= tosolution(coeff);
    %dist2 = cc'*H*cc -2*cc'*g +C;
    %display(dist2);
    % sol = H\g;
    % dist3 = sol'*H*sol - 2*sol'*g +C;
    %display(dist3)
    l=mx*my;
    H = sparse(2*l,2*l);
    
    H(1:l,1:l)=Hpe;
    H(l+(1:l),l+(1:l))=Hpe;
    %E_old = cc'*H*cc -2*cc'*g +C;
    % 
    %sol = H\g;
    sol = quadprog(H,zeros(size(cc)),[],[],Aeq,beq,[],[],[],opts);
    %H(1:l,1:l)=H(1:l,1:l)+mu*Hpe;
    %H(l+(1:l),l+(1:l))=H(l+(1:l),l+(1:l))+mu*Hpe;
    %E_new = sol'*H*sol - 2*sol'*g +C;
%     if E_new>E_old - 0.001 
%         break;
%     end
    %sol= 0.8*sol+0.2*cc;
    % coeff2=zeros(size(coeff));
    % coeff2(1,:,:)=reshape(sol(1:l),mx,my);
    % coeff2(2,:,:) = reshape(sol(l+(1:l)),mx,my);
    coeff = tocoeff(sol,mx,my,2);
end

[x,y] = tform.transformPointsInverse(coeff(1,:,:),coeff(2,:,:));
coeff(1,:,:) =x;
coeff(2,:,:) =y;
coeff(:,1,:) = coeff_origin(:,1,:);
coeff(:,end,:) = coeff_origin(:,end,:);
coeff(:,:,1) = coeff_origin(:,:,1);
coeff(:,:,end) = coeff_origin(:,:,end);


end

function n = get_outer_normal(points,id)
dim = 2;
num = numel(id);
total = size(points,2);
n = zeros(dim,num);
for ii = 1:num
    T = points(:,mod(id(ii)+1-1,total)+1)-points(:,mod(id(ii)-1-1,total)+1);
%     cT = T(1)+1i*T(2);
%     cN = cT*(-1i);
    N = [T(2);-T(1)];
    n(:,ii)=N/norm(N);
end
end
