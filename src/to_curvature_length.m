function mat_ret = to_curvature_length(mat)
% convert x y coordinates to curvature and length of each edge
% mat:2*1024
N = size(mat,2);
%m = size(m,3);
next = [2:N,1];
prev = [N,1:(N-1)];
fwd = mat(:,next,:)-mat;
bak = mat - mat(:,prev,:);
% x1 = sum(fwd.*bak,1); % cos ...
% fwd[:,:,1]*bak[:,:,0]-fwd[:,:,0]*bak[:,:,1]
x2 = fwd(2,:,:).*bak(1,:,:)-fwd(1,:,:).*bak(2,:,:);% sin, for the sign of the angle
mat_ret = zeros(size(mat));
%mat_ret(1,:,:)=atan2(x2,x1); % unstable
mat_ret(1,:,:)=sign(x2).* stable_ang(bak,fwd);
mat_ret(2,:,:) = vecnorm(fwd,2,1);

end

function ang = stable_ang(v1,v2)
% https://scicomp.stackexchange.com/questions/27689/numerically-stable-way-of-computing-angles-between-vectors
% v1, v2 : dim * N
a = vecnorm(v1,2,1);
b = vecnorm(v2,2,1);
c = vecnorm(v1-v2,2,1);
I = a >=b;
b2 = b(~I);
b(~I) = a(~I);
a(~I) = b2;
% exchange values so that a >=b
II = b>=c;
mu = zeros(size(a));
mu(II) = c(II)- (a(II)-b(II));
mu(~II) = b(~II)-(a(~II)-c(~II));
ang = 2* real(atan(sqrt(mu.*((a-b)+c)./((a+(b+c)).*((a-c)+b)))));
end