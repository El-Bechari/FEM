
function [lambda,gradx,grady]=ref_base_func_P1(x)
lambda=[1-x(:,1)-x(:,2)  x(:,1)  x(:,2)];
gradx=[-1 1 0];
grady=[-1 0 1];
end