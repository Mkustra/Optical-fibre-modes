function F = fun(x)

global m
global v

F(1)=(x(1).*( besselj(m-1,x(1))./besselj(m,x(1)) ) ) + (  x(2).* (besselk(m-1,x(2))./besselk(m,x(2)))  );
F(2)=v.^2-x(1).^2-x(2).^2;