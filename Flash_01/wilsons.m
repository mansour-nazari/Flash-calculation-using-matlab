function K=wilsons(T,P,Tc,Pc,omega)
Tc=Tc+460;% Convet to R
% initial geuss using wilsons equation
K=Pc/P.*exp(5.37*(1+omega).*(1-Tc./T));

end