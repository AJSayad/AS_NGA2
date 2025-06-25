% mesh and domain testing and calculations
% this script is intended to be a pre-simulation mesh and domain check
clear; close all; clc; 

%% CASE: AIAAWG shot 6B 100CPD

markersize = 1;
plt_bounds = 1; % flag to plot domain bounds
plt_delta =  1; % flag to plot shock thickness
plt_mesh =   1; % flag to plot mesh

% quickly calculate the extents of the domain and visualize the domain before running a full simulation

%% input file
CFLa = 0.5;            % CFL restriction
tmax = 1e-4;           % max sim time
rho_liq = 1000;        % liquid density
d0 = 0.00127;          % droplet diameter
CPD = 100;             % resolution (cells/d0)
DX = 12;               % # of droplet diameter lengths in x direction
DY = 4;                % # of droplet diameter lengths in y direction (total domain height)
% note: set DZ = 0 for a 2D simulation
DZ = 4;                % # of droplet diameter lengths in z direction (total domain depth)
DX_stretch = 2;        % number of droplet diameter lengths to add at the end of the x direction
DY_stretch = 4;        % number of droplet diameter lengths to add at the top AND bottom of y direction
DZ_stretch = 4;        % number of droplet diameter lengths to add at the top AND bottom of z direction
dctrx = 1.5*d0;        % droplet center in x
dctry = 0.0;           % droplet center in y 
dctrz = 0.0;           % droplet center in z
xshock = dctrx-0.6*d0; % location of the shock

[x,y,z,nx,nx_stretch,ny,ny_stretch,nz,nz_stretch,dx_uni,dy_uni,dz_uni] = generate_mesh(d0,CPD,DX,DY,DZ,DX_stretch,DY_stretch,DZ_stretch);
cell_total = (nx+nx_stretch)*(ny+2*ny_stretch)*(nz+2*nz_stretch)

%% calculate shock thickness 
nshock = 3;          % cells for shock thickening (nshock to the left and right of center)
dshock = 2*3*dx_uni; % full shock thickness
delta = linspace(xshock - dshock/2, xshock + dshock/2, 2*nshock); % extent of shock for plotting
shock_height = zeros(2*nshock); % plot shock at centerline

%% normal shock calculations
% verified with https://devenport.aoe.vt.edu/aoe3114/calc.html
gamma = 1.3; % ratio of specific heats
% pre shock conditions %https://www.calctool.org/atmospheric-thermodynamics/air-density
rho1 = 1.17134;           % pre shock density
P1   = 101000;            % pre shock pressure
T1   = 299;               % pre shock temperature
u_s  = 1052;              % shock velocity
c1 = sqrt(gamma*P1/rho1); % pre shock sound speed
M1 = u_s/c1;              % shock mach number
% post shock conditions from normal shock relations
rho2 = (rho1*(gamma + 1)*M1^2)/(2 + (gamma - 1)*M1^2)
P2   = P1*(1 + (2*gamma/(gamma + 1))*(M1^2 - 1))
T2   = T1*(1 + (2*gamma/(gamma + 1)*(M1^2 - 1)))*((2 + (gamma - 1)*M1^2)/((gamma + 1)*M1^2))
c2   = sqrt(gamma*P2/rho2)
M2   = sqrt((1 + 0.5*(gamma - 1)*M1^2)/((gamma*M1^2) - 0.5*(gamma-1))) 
% gas velocity in shock frame
u_2SF= c2*M2
% gas velocity behind the shock in lab frame
u_2LF= u_s - u_2SF
M2LF = u_2LF/c2

%% time shock hits droplet
edgedistance = (dctrx - 0.5*d0) - xshock;
t_impact = edgedistance / u_2LF
dt = (dx_uni*CFLa)/(u_2LF+c1) % calculate timestep in [s]

%% non dimensional time approximation
tc2 = sqrt(rho_liq/rho2)*(d0/u_2LF) % characteristic breakup time
tau_max = tmax/tc2                  % max non-dim time

%% approximate bow shock stand off distance 
% based on relation from https://doi.org/10.1063/1.323828this stand off distance from solid cylinder
xstandoff = d0*(1 - (sqrt(1 - 4/(gamma+1)^2)*((1/M2LF) + (gamma-1)/2))) % standoff distance from the droplet edge to bowshock
xbow = (dctrx - 0.5*d0) - xstandoff
ybow = 0.5*u_2LF*tmax % (over) approximate how far the bow shock will travel in y 

%% generate a circle to represent the droplet
npts = 100; % set number of points for circle
theta = linspace(0,2*pi,npts); % define angle for circle
dropletx = dctrx + (d0/2)*cos(theta); % x coordinates of droplet from parametric equation
droplety = dctry + (d0/2)*sin(theta); % y coordinates of droplet from parametric equation
dropletz = dctrz + (d0/2)*sin(theta); % z coordinates of droplet from parametric equation

%% generate uniform domain boundaries
left_boundx = zeros(1,ny);
left_boundy = linspace(-(DY*d0)/2,(DY*d0)/2,ny);
right_boundx = (DX*d0)*ones(1,ny);
right_boundy = left_boundy;
bot_boundx = linspace(0,DX*d0,nx);
bot_boundy = (-(DY*d0)/2)*ones(1,nx);
top_boundx = bot_boundx;
top_boundy = ((DY*d0)/2)*ones(1,nx);

if plt_mesh == 1
    %% plot mesh and domain
    [X,Y] = meshgrid(x,y);
    figure
    plot(X,Y,'bo','MarkerSize',markersize)
    hold on
    xline(xshock,'r','LineWidth',3)
    xline(xbow,'g','LineWidth',3) % plot (over)approximation of standoff distance
    yline(ybow,'g','LineWidth',3)
    yline(-ybow,'g','LineWidth',3)

    plot(dropletx,droplety,'c','LineWidth',2)
    if plt_delta==1
        plot(delta,shock_height,'y','linewidth',5)
    end
    if plt_bounds == 1
        plot(left_boundx,left_boundy,'k','linewidth',3)   % left bound of uni region
        plot(right_boundx,right_boundy,'k','linewidth',3) % right bound of uni region
        plot(bot_boundx,bot_boundy,'k','linewidth',3)     % bottom bound of uni region
        plot(top_boundx,top_boundy,'k','linewidth',3)     % top bound of uni region
    end
    axis equal
    title('Mesh check in X-Y plane'); xlabel('Mesh in X'); ylabel('Mesh in Y')
    text(xshock,0.25*((DY*d0)/2),'\leftarrow initial shock','FontSize',10,'FontWeight','bold')
    text(xbow,0.25*((DY*d0)/2),'\leftarrow approx bow shock','FontSize',10,'FontWeight','bold')

    if nz>1
        [X,Z] = meshgrid(x,z);
        figure
        plot(X,Z,'bo','MarkerSize',markersize)
        hold on
        xline(xshock,'r','LineWidth',3)
        plot(dropletx,droplety,'c','LineWidth',2)
        axis equal
        title('Mesh check in X-Z plane'); xlabel('Mesh in X'); ylabel('Mesh in Z')

        [Y,Z] = meshgrid(y,z);
        figure
        plot(Y,Z,'bo','MarkerSize',markersize)
        hold on
        xline(xshock,'r','LineWidth',3)
        plot(dropletx,droplety,'c','LineWidth',2)
        axis equal
        title('Mesh check in Y-Z plane'); xlabel('Mesh in Y'); ylabel('Mesh in Z')
    end
end