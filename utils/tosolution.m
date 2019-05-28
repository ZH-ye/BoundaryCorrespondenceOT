function [solution,boundary]=tosolution(coeff)
[dim,mx,my] = size(coeff);
b = zeros(size(coeff));
b(:,:,1)=1;
b(:,end,:)=1;
b(:,1,:)=1;
b(:,:,end)=1;
solution = zeros([dim*mx*my,1]);
boundary = zeros([dim*mx*my,1]);
for i = 1:dim
    tmp = reshape(b(i,:,:),mx,my);
    boundary((i-1)*mx*my + (1:mx*my) ) = tmp(:);
end
for i = 1:dim
    tmp = reshape(coeff(i,:,:),mx,my);
    solution((i-1)*mx*my + (1:mx*my) ) = tmp(:);
end
end