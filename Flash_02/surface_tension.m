function S=surface_tension(x,y,MW,po,pg,pchi)

%Macleod-Sugden correlation

Mo=sum(x.*MW);
Mg=sum(y.*MW);

A=po*0.0160185/(Mo);
B=pg*0.0160185/(Mg);



S=sum(pchi.*(A*x-B*y));
S=S^4;