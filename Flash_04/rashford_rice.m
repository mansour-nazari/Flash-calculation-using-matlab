function nv=rashford_rice(K,zi,nv)

n=length(K);

            tol=10^-8;
            f=1;
            it=1;
            while abs(f)>tol
                it=it+1;
                f=0;
                df=0;
                for t=1:n
                    f=zi(t,1)*(K(t)-1)/(nv*(K(t)-1)+1)+f;
                    df=-(zi(t,1)*(K(t)-1)^2)/((nv*(K(t)-1)+1)^2)+df;
                end
                nv=nv-f/df;
            end
            
end