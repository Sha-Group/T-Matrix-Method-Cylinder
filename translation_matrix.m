% Translation matrix£¬from Hankel to Bessel function

function output=translation_matrix(xi,yi,xj,yj,MAX,k_rho)
d=sqrt((xi-xj)^2+(yi-yj)^2);
phi=atan2(yi-yj,xi-xj);
%phi=atan2(yj-yi,xj-xi);
Alpha_ij=zeros(2*MAX+1,2*MAX+1);
for m=1:2*MAX+1
    for n=1:2*MAX+1
        Alpha_ij(m,n)=cb2(n-m,k_rho,d)*exp(1j*(n-m)*phi);%
    end
end
Alpha=[Alpha_ij,zeros(2*MAX+1,2*MAX+1);zeros(2*MAX+1,2*MAX+1),Alpha_ij];
output=Alpha; 

% Cylindrical Bessel 2 (H^1_n)
function output=cb2(degree,k_rho,rho)
output=besselh(degree,1,k_rho*rho);  