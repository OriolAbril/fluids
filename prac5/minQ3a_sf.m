function Ra= minQ3a_sf(k)
global Raoldsf
sf=@(Ra)funQ3a_sf(Ra,k);
Raoldsf=fsolve(sf,Raoldsf);
Ra=Raoldsf;
end

