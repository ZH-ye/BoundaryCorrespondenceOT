function coeff = tocoeff(solution,mx,my,dim)

l= mx*my;
coeff = zeros(dim,mx,my);
for i = 1:dim
    coeff(i,:,:) = reshape(solution((i-1)*l + (1:l)) ,mx,my);
end

end