function Ra= minQ3a_ns(k)
global Raoldns
ns=@(Ra)funQ3a_ns(Ra,k);
Raoldns=fsolve(ns,Raoldns);
Ra=Raoldns;
end

