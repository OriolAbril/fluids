clear all;
close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PHYSICAL PARAMETERS
% Coriolis parameter % Gravity % Reference depth % wave speed
f = 10^(-4);           g = 10;   h0 = 50;          c0=sqrt(g*h0);
% Deformation radius % Maximum free surface for initial condition
R = c0/f;              ampli = 0.001*h0; 
% (Use small amplitude with respect to depth because we solve a linearized model)

% NUMERICAL SPATIAL PARAMETERS
% Grid size     % Domain size           % Distance between grid points in x and y
M = 60; N = 50;   Lx = 12*R; Ly = 10*R;   dx = Lx/(M-2); dy = Ly/(N-2);
% (Boundary conditions in x and y require ghost points at i=1 and j=1
% Jumps (in number of points) in x and y to plot u and v
dxp = 3; dyp = 3;

% NUMERICAL TEMPORAL PARAMETERS
% Final time                % Time step for ploting 
TF = 2*pi*Lx/sqrt(g*h0);      dtp = 10;
% Time step as a function of Courant parameter
dt = min(dx,dy)/c0/2; % ??? -> WRITE IT AS A FUNCTION OF THE CFL CONDITION, USING min(dx,dy) TO REPRESENT THE SPATIAL GRID SIZE
% Number of time steps as a function of time step and final time
Nsteps = floor(TF/dt);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DECLARE ARRAYS
% Coordinates
x = (((1:M)-3/2)*dx); y = (((1:N)-3/2)*dy);
% Zero in j=3/2. First line (j=1) is a ghost point for BC
% Zero in i=3/2. First column (i=1) is a ghost point for BC
% Velocities
u = zeros(M,N); v = zeros(M,N);
% New velocities
un=u; vn=v;
% Sea surface elevation % New elevation % Unpertubed depth
zeta = zeros(M,N);        zetan = zeta;   h = zeta;
% Variables for plotting
Mp = floor(M/dxp); Np = floor(N/dyp);
xp = zeros(1,Mp); yp = zeros(1,Np);
% Land-sea mask to apply boundary conditions on the velocities
issea=ones(M,N);

% Coordinates of cell centers for plotting velocity vectors (also avoiding ghost points)
xp(1,:) = x(1,2:dxp:M-1);
yp(1,:) = y(1,2:dyp:N-1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Staggered C-grid:
% zeta(i,j) is at node i,j
% h(i,j) is at node i,j
% u(i,j) is in i-1/2
% v(i,j) is in j-1/2

t = 0;

% Initial condition: a sinusoidal Kelvin wave
kx = 2*pi/(Lx); % ??? -> WRITE THE WAVELENGTH
[YY,XX]=meshgrid(y,x);
F1=sin(kx*XX);
zeta=ampli*exp(-YY/R).*F1;
u=c0/h0*zeta;

% Bathymetry of the sea
h(:,:) = h0;

% Definition of the land regions with the land-sea mask: here land points in j=1 and j=N 
issea(:,1) = 0;
issea(:,N) = 0;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHIC INTERFASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kelvin_graphic(zeta,u,v,x,y,xp,yp,dxp,dyp,N,M,Np,Mp,'To start press any key')
pause

tic %To know the elapsed time for the resolution and plot 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESOLUTION AND PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Staggered C-grid:
% zeta(i,j) is at node i,j
% h(i,j) is at node i,j
% u(i,j) is in fact at node i-1/2,j
% v(i,j) is in fact at node i,j-1/2

% Make Nsteps time steps
zetavec3=zeros(1,Nsteps);
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
    zetan(:,1) = 0; % ??? -> WRITE THE zeta BC AT GHOST zeta NODE j=1
    % Boundary condition at j=N (real node)
    zetan(:,N) = 0; % ??? -> WRITE THE zeta BC AT REAL zeta NODE j=N
        
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
    zetavec3(nn)=max(zeta(:,3));
    
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
                    un(i,j) = u(i,j)+dt*(dxz+vi12j); % ??? -> FILL IT IN FROM THE x-MOMETUM EQUATION
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
            un(:,N) = 0; % ??? -> WRITE THE u BC AT REAL u NODE j=N
            
            % Update u
            u = un;
            
        elseif multipas(k) == 2 %This is evaluated first if nn even
            
            for i=2:M-1
                for j=2:N
                    dyz=-g*(zeta(i,j)-zeta(i,j-1))/dy;
                    uij12=-f/4*(u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1));
                    vn(i,j) = v(i,j)+dt*(dyz+uij12); % ??? -> FILL IT IN FROM THE y-MOMETUM EQUATION
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRAPHIC INTERFASE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Every once in a while make a plot
    if mod(nn,dtp)==0 %Plot every dtp time steps
        titol=['Elevation at time: ',num2str(t)];
        kelvin_graphic(zeta,u,v,x,y,xp,yp,dxp,dyp,N,M,Np,Mp,titol)
        pause(0.0001);       
    end
    
end

toc %To know the elapsed time for the resolution and plot 
tvec=linspace(0,TF,Nsteps);
figure(2)
plot(tvec,zetavec3)
xlabel('Time')
ylabel(' max(eta) ; j=3')
grid on
indexsum
cmean=indexsum*dx/TF
relerr=100*abs(cmean-c0)/c0