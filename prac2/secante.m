%Metodo de la secante
function sol=secante(a,b,tol,imax,fun)

for i=1:imax
    fa=fun(a);
    fb=fun(b);
    x=b-fb*(b-a)/(fb-fa);
    a=b;
    b=x;
    if abs(b-a)<tol;         
        break
    end
end
sol=x;
return


