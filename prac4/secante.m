%Metodo de la secante
function sol=secante(a,b,tol,imax,fun,params)

for i=1:imax
    fa=fun(a,params);
    fb=fun(b,params);
    x=b-fb*(b-a)/(fb-fa);
    a=b;
    b=x;
    if abs(b-a)<tol;         
        break
    end
end
sol=x;
return


