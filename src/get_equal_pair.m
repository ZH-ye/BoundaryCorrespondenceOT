function f = get_equal_pair(x)
n = size(x,1);
el = [x(:,2)-x(:,1),x(:,3)-x(:,2),x(:,4)-x(:,3),1+x(:,1)-x(:,4)];

f = 0.5*(el(:,1)+el(:,3)-el(:,2)-el(:,4)).^2;
end