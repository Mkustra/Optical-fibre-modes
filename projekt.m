clear all
close all 
clc

global v;
global m;

% Dane wejsciowe
wektror=zeros(1,6);
('WprowadŸ wektor [d³. fali [nm], promieñ rdzenia [um], wsp. za³amania p³aszcza, wsp. za³amania rdzenia, wartoœæ m, wartoœæ p]')

wektor  = input('');
lambda1 = wektor(1);
a1      = wektor(2);
n       = wektor(3);
nc      = wektor(4);
m       = wektor(5);
p       = wektor(6);

a=a1*1e-6;
lambda=lambda1.*1e-9;

if n>=nc
    ('Niepoprawne dane wejœciowe')
    return
end

% Czestotliwosc znormalizowana
v=((2.*pi.*a)./lambda).*sqrt(nc.^2-n.^2);

%Wyznaczenie rozwiazan ukladu rownan
u0=linspace(1e-3.*v,v-1e-3.*v, ceil(v));
w0=sqrt(v.^2-u0.^2);
n=1;

for i=1:length(u0)
        xprim = fsolve(@fun, [u0(i), w0(i)]);
        if imag(xprim(1))==0 && imag(xprim(2))==0 && xprim(1)>=0 && xprim(2)>=0
            x(:,n) = xprim; 
            n=n+1;
        end
end

if x == -3 
    ('Numeryczny b³¹d obliczeñ')
    return
end

x=round(10.*x)/10;
x=unique(x','rows');
u=x(:,1);
w=x(:,2);

%Obliczanie pola elektrycznego
if length(u)>=p
x=linspace(-2*a,2*a,150);

for i=1:length(x)
    for j=1:length(x)
       [fi,r]=cart2pol(x(i),x(j));
        if r<=a
            Ey(i,j)=( besselj(m, (u(p).*r) ./a) ./ ( besselj(m,u(p))  ) ) .* cos(m.*fi);

        else
            Ey(i,j)=( besselk(m, (w(p).*r) ./a) ./ ( besselk(m,w(p))  ) ) .* cos(m.*fi);
        end
    end
end

figure(1)
surf(x,x,Ey)
figure(2)
surf(x,x,Ey)
view(0,90)

else
 ('Podany mod nie rozchodzi sie w œwiat³owodzie o zadanych parametrach')
end