function coeff=spline_parametrization(points,c,iter)
sh = spline_helper;
sh.set_points(points);
if nargin < 3
    iter = 10;
end
sh.try_update_spline(c,iter)
% sh.current_fp=c;
% sh.fit_four_boundary();
% sh.refine(10);
coeff = sh.coeff;
end