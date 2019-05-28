function mu = spline_distortion_mu(sp,X)
% sp: spline object generated by spmak
% X: 2*N
% outputs:
% mu: 1*N complex
dx = fnval(fnder(sp,[1,0]),X);
dy = fnval(fnder(sp,[0,1]),X);
E = dot(dx,dx,1);
F = dot(dx,dy,1);
G = dot(dy,dy,1);
mu = (E-G+2*1i*F)./(E+G+2*sqrt(E.*G-F.*F));
end