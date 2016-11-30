%codi per anar baixnt el reynolds poc a poc fins a trobar el mínim 
%reynolds on encara es veuen TSW


L0=6.1570;
L0name=num2str(L0+1e-12,'%5.2f');
ii=find(L0name == '.') ; L0name(ii)='p';
Reold=4500;
Restep=50;
Re=4450;
Tintegration=20000; 
for ii=1:50
    Re
    Reold
    filestate = ['A_L' L0name '_Re' int2str(Reold) '_T' int2str(Tintegration)];
    Tintegration=20000;
    ppftstep_opt_save(filestate,1,[0 Tintegration],L0,Re)
    fileampli = ['SerTemp_L' L0name 'Re'  int2str(Re)];
    load(fileampli)
    if mean(avec(end-100:end))<0.01
        if Restep~=1
            Restep=max(round(Restep/2),1);
            Re=Reold-Restep;
        else
            disp('Ha relaminaritzat baixant de un en un. Abortar')
            break
        end
    else
        Reold=Re;
        Re=Reold-Restep;
    end
end
