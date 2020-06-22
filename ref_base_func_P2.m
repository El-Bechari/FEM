function [lambda,gradx,grady]=ref_base_func_P2(xx)
A=[ 1     0     0     0     0     0
    -3    -1     0     4     0     0
    2     2     0    -4     0     0
    -3     0    -1     0     0     4
    2     0     2     0     0    -4
    4     0     0    -4     4    -4];
x=xx(:,1);y=xx(:,2);
zx=0*x;
lambda= [1+zx x x.^2 y y.^2 x.*y]*A;

gradx=  [zx 1+zx 2*x zx zx y]*A;
grady=  [zx zx zx 1+zx 2*y x]*A;
end