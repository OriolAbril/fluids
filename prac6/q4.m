% If we keep constant t_max and dt:
% then if dt=C/1.5 then C=min(dx,dy)/c0
% min(dx,dy)/c0
clear all;
close all;
% Now we are going to keep constant M=60; N=50
% We vary the friction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PHYSICAL PARAMETERS
% Coriolis parameter
f = 10^(-4);
% Gravity
g = 10;
% Reference depth
h0 = 50;
%wave speed
c0=sqrt(g*h0);
% Deformation radius
R = c0/f; % ??? -> WRITE THE FORMULA FOR THE ROSSBY RADIUS OF DEFORMATION
% Maximum free surface for initial condition
ampli = 0.001*h0; 
% (Use small amplitude with respect to depth because we solve a linearized model)

% NUMERICAL SPATIAL PARAMETERS
% Grid size
M = 60;
N = 50;
% Domain size
Lx = 12*R;
Ly = 10*R;
% Distance between grid points in x and y
dx = Lx/(M-2);
dy = Ly/(N-2);
% (Boundary conditions in x and y require ghost points at i=1 and j=1
% Jumps (in number of points) in x and y to plot u and v
dxp = 3;
dyp = 3;

% NUMERICAL TEMPORAL PARAMETERS
% Final time
TF =10* Lx/c0;
% Time step as a function of Courant parameter
C=min(dx,dy)/c0;
% S'ha de trobar el dt maxim possible perque hi hagi estabilitat
% Es pregunta pel C maxim;
%C=0.3; 
r=1e-3; % per a l'apartat quatre utilitzarem r=[1e-5, 1e-4, 1e-3]
dt = C/2; % ??? -> WRITE IT AS A FUNCTION OF THE CFL CONDITION, USING min(dx,dy) TO REPRESENT THE SPATIAL GRID SIZE
% Number of time steps as a function of time step and final time
Nsteps = floor(TF/dt);
% Time step for ploting 
dtp = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DECLARE ARRAYS
% Coordinates
x = zeros(1,M);
y = zeros(1,N);
% Velocities
u = zeros(M,N);
v = zeros(M,N);
% New velocities
un=u;
vn=v;
% Sea surface elevation
zeta = zeros(M,N);
% New elevation
zetan = zeta;
% Unpertubed depth
h = zeta;
% Variables for plotting
Mp = floor(M/dxp);
Np = floor(N/dyp);
xp = zeros(1,Mp);
yp = zeros(1,Np);
% Land-sea mask to apply boundary conditions on the velocities
issea=ones(M,N);

% Coordinates of cell centers
% Zero in i=3/2. First column (i=1) is a ghost point for BC
for i=1:M
    x(i) = (i-3/2)*dx;
end
% Zero in j=3/2. First line (j=1) is a ghost point for BC
for j=1:N
    y(j) = (j-3/2)*dy;
end

% Coordinates of cell centers for plotting velocity vectors (also avoiding ghost points)
xp(1,:) = x(1,2:dxp:M-1);
yp(1,:) = y(1,2:dyp:N-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Staggered C-grid:
% zeta(i,j) is at node i,j
% h(i,j) is at node i,j
% u(i,j) is in i-1/2
% v(i,j) is in j-1/2

t = 0;

% Initial condition: a sinusoidal Kelvin wave
kx = 2*pi/(Lx); % ??? -> WRITE THE WAVELENGTH
for i=1:M
    for j=1:N
        F1 = sin(kx*x(i)); % ??? -> WRITE THE F1 FUNCTION
        zeta(i,j) = ampli*exp(-y(j)/R)*F1; % ??? -> WRITE THE eta OF A KELVIN WAVE
    end
end

for i=1:M
    for j=1:N
        F1 = sin(kx*x(i)); % ??? -> WRITE THE F1 FUNCTION
        u(i,j) =ampli*c0/h0*exp(-y(j)/R)*F1; % ??? -> WRITE THE u OF A KELVIN WAVE
    end
end

% Bathymetry of the sea
h(:,:) = h0;

% Definition of the land regions with the land-sea mask: here land points in j=1 and j=N 
issea(:,1) = 0;
issea(:,N) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPHIC INTERFASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables for plotting
zetap = zeros(N,M); zetap(:,:) = NaN;
up = zeros(Np,Mp); up(:,:) = NaN;
vp = zeros(Np,Mp); vp(:,:) = NaN;
% Plot zeta everywhere except at the ghost points
for i=2:M-1
    zetap(2:N-1,i)=zeta(i,2:N-1);
end
% Calculate uc vc at the cell centers and only every dp points
for i=2:dxp:M-1
    ii=1+(i-2)/dxp;
    up(1:Np,ii)=(u(i,2:dyp:N-1)+u(i+1,2:dyp:N-1))/2;
end
for j=2:dyp:N-1
    jj=1+(j-2)/dyp;
    vp(jj,1:Mp)=(v(2:dxp:M-1,j)+v(2:dxp:M-1,j+1))/2;
end

% Show initial free surface and velocity in a contour plot
figure
pcolor(x,y,zetap)
shading flat
colorbar
hold on
curr=quiver(xp,yp,up,vp,0.8,'k');
hold off
xlabel('Alongshore coordinate (m)')
ylabel('Cross-shore coordinate (m)')
title ('To start press any key')
pause
clear zetap up vp


tic %To know the elapsed time for the resolution and plot 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESOLUTION AND PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Staggered C-grid:
% zeta(i,j) is at node i,j
% h(i,j) is at node i,j
% u(i,j) is in fact at node i-1/2,j
% v(i,j) is in fact at node i,j-1/2

% Make Nsteps time steps
zetavec3=[];
indexsum=0;
index=0;
for nn=1:Nsteps
    
    % EQUATION 1: MASS CONSERVATION -> provides zeta
    
    % New free surface elevation
    for i=2:M-1
        for j=2:N-1
            %zetan(i,j) = zeta(i,j)-dt*h0*((u(i+1,j)-u(i,j))/dx+(v(i,j+1)-v(i,j))/dy); % if h ct MASS CONSERVATION EQUATION
            dxu=(u(i+1,j)-u(i,j))/dx;
            dyv=(v(i,j+1)-v(i,j))/dy;
            dhx=(h(i+1,j)-h(i-1,j))/(2*dx);
            dhy=(h(i,j+1)-h(i,j-1))/(2*dy);
            uij=(u(i+1,j)+u(i,j))/2;
            vij=(v(i,j+1)+v(i,j))/2;
            zetan(i,j) = zeta(i,j)-dt*(h(i,j)*(dxu+dyv)+uij*dhx+vij*dhy);
        end
    end
    % Boundary condition at i=1 (ghost node)
    zetan(1,:) = zetan(M-1,:);
    % Boundary condition at i=M (real node)
    zetan(M,:) = zetan(2,:);
    % Boundary condition at j=1 (ghost node)
    zetan(:,1) = zeta(:,1); % ??? -> WRITE THE zeta BC AT GHOST zeta NODE j=1
    % Boundary condition at j=N (real node)
    zetan(:,N) = zeta(:,N); % ??? -> WRITE THE zeta BC AT REAL zeta NODE j=N
        
    % Update zeta
    zeta = zetan;
    
    indexin=index;
    [etamax,index]=max(zeta(:,3));
    index1=index-indexin;
    if abs(index1)>1
        index1=1;
    end
    if index1<0
        index1=1;
    end
    indexsum=indexsum+index1;
    zetavec3=[zetavec3 max(zeta(:,3))];
    
    % EQUATION 2: NAVIER STOKES -> provides u and v
    
    % New velocity values with fractional-step approach on Coriolis
    if mod(nn,2) == 1 % Odd values of nn, including the first time step 
        multipas = [1 2];
    elseif mod(nn,2) == 0 % Even values of nn
        multipas = [2 1];
    end
    
    for k=1:2
        if multipas(k) == 1 %This is evaluated first if nn odd
            
            for i=2:M-1
                for j=1:N-1
                    dxz=-g*(zeta(i,j)-zeta(i-1,j))/dx;
                    vi12j=f/4*(v(i,j)+v(i,j+1)+v(i-1,j+1)+v(i-1,j));
                    fricu=-(r/h0)*u(i,j);
                    un(i,j) = u(i,j)+dt*(dxz+vi12j+fricu); % ??? -> FILL IT IN FROM THE x-MOMETUM EQUATION
                    % Applying a mask on sea-land interfaces: it makes un=0 if there  
                    % is a land point on the right or on the left of this u node
                    un(i,j)=un(i,j)*issea(i,j)*issea(i-1,j);
                end
            end
            % Boundary condition at i=1 (ghost node)
            un(1,:) = un(M-1,:); % ??? -> WRITE THE u BC AT GHOST u NODE i=1
            % Boundary condition at i=M (real node)
            un(M,:) = un(2,:); % ??? -> WRITE THE u BC AT REAL u NODE i=M
            % Boundary condition at j=N (real node)
            un(:,N) = u(:,N); % ??? -> WRITE THE u BC AT REAL u NODE j=N
            
            % Update u
            u = un;
            
        elseif multipas(k) == 2 %This is evaluated first if nn even
            
            for i=2:M-1
                for j=2:N
                    dyz=-g*(zeta(i,j)-zeta(i,j-1))/dy;
                    uij12=-f/4*(u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1));
                    fricv=-(r/h0)*v(i,j);
                    vn(i,j) = v(i,j)+dt*(dyz+uij12+fricv); % ??? -> FILL IT IN FROM THE y-MOMETUM EQUATION
                    % Applying a mask on sea-land interfaces: it makes vn=0 if there  
                    % is a land point on the top or on the bottom of this v node
                    vn(i,j)=vn(i,j)*issea(i,j)*issea(i,j-1);
                end
            end
            % Boundary condition at i=1 (ghost node)
            vn(1,:) = vn(M-1,:); % ??? -> WRITE THE v BC AT GHOST v NODE i=1
            % Boundary condition at i=M (real node)
            vn(M,:) = vn(2,:); % ??? -> WRITE THE v BC AT GHOST v NODE i=M
            % Boundary condition at j=1 (ghost node)
            vn(:,1) = 0; % ??? -> WRITE THE v BC AT GHOST v NODE j=1
            
            % Update v
            v = vn;
            
        end
    end
    
    t = t+dt;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% GRAPHIC INTERFASE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Every once in a while make a plot
    if mod(nn,dtp)==0 %Plot every dtp time steps
        
        % Variables for plotting
        zetap = zeros(N,M); zetap(:,:) = NaN;
        up = zeros(Np,Mp); up(:,:) = NaN;
        vp = zeros(Np,Mp); vp(:,:) = NaN;
        % Plot zeta everywhere except at the ghost points
        for i=2:M-1
            zetap(2:N-1,i)=zeta(i,2:N-1);
        end
        % Calculate uc vc at the cells center and only every dp points
        for i=2:dxp:M-1
            ii=1+(i-2)/dxp;
            up(1:Np,ii)=(u(i,2:dyp:N-1)+u(i+1,2:dyp:N-1))/2;
        end
        for j=2:dyp:N-1
            jj=1+(j-2)/dyp;
            vp(jj,1:Mp)=(v(2:dxp:M-1,j)+v(2:dxp:M-1,j+1))/2;
        end
        
        % Show initial free surface and velocity in a contour plot
        pcolor(x,y,zetap)
        shading flat
        colorbar
        hold on
        curr=quiver(xp,yp,up,vp,0.8,'k');
        hold off
        xlabel('Alongshore coordinate (m)')
        ylabel('Cross-shore coordinate (m)')
        title(['Elevation at time: ',num2str(t)])
        pause(0.0001);
        clear zetap up vp
        
    end
    
end

toc %To know the elapsed time for the resolution and plot 
tvec=linspace(0,TF,TF/dt);
figure(2)
plot(tvec,zetavec3)
xlabel('Time')
ylabel(' max(eta) ; j=3')
grid on
cmean=indexsum*dx/TF
relerr=100*abs(cmean-c0)/c0