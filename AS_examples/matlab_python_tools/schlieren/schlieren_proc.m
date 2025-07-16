%% Post processing example for matlab
clear all; close all; clc; plot_formatting; format long

time_flag = 1; % compute non-dimensional time?

%% mesh grid
d0 = 0.00187;
alpha = 1.03;          % stretching ratio
Lx =           12*d0    % extent of uniform region
Lx_ref =       0;      % start of uniform region
nx =           1200;    % number of cells in uniform region
nx_stretchL =  0;      % cells to add at left
nx_stretchR =  75;     % cells to add at right
dx_uni =       Lx/nx;  % uniform region spacing

Ly =         4*d0      % absolute domain height (top = 0.5*Ly, bottom = -0.5*Ly)
ny =         400;      % number of cells in uniform region
ny_stretch = 175;      % number of cells to add to top and bottom (total added = 2*ny_stretch)
dy_uni = Ly/ny;        % uniform region spacing

Lz = dx_uni;
nz = 1;
dz = Lz/nz;

[x,y] = genmesh(Lx_ref,nx,nx_stretchL,nx_stretchR,dx_uni,ny,ny_stretch,dy_uni,alpha);

x = x(1:end-1); % trim off last point (since we define cell edges)
y = y(1:end-1);
[YY,XX] = meshgrid(y,x);
XX = XX/d0;
YY = YY/d0;

%% schileren constants
Kliq = 5;
Kgas = 15;
const = 2.5e5;

%% characteristic breakup time
if time_flag == 1
    % normal shock calculations
    % verified with https://devenport.aoe.vt.edu/aoe3114/calc.html
    gamma = 1.3; % ratio of specific heats
    % pre shock conditions
    % rho1 = 1.176;
    rho_liq = 1000; % liquid density
    rho1 = 1.17617; %https://www.calctool.org/atmospheric-thermodynamics/air-density
    P1   = 101000;
    T1   = 297;
    u_s  = 1775; % shock velocity
    c1 = sqrt(gamma*P1/rho1); % pre shock sound speed
    M1 = u_s/c1; % shock mach number
    % post shock conditions
    rho2 = (rho1*(gamma + 1)*M1^2)/(2 + (gamma - 1)*M1^2);
    P2   = P1*(1 + (2*gamma/(gamma + 1))*(M1^2 - 1));
    T2   = T1*(1 + (2*gamma/(gamma + 1)*(M1^2 - 1)))*((2 + (gamma - 1)*M1^2)/((gamma + 1)*M1^2));
    c2 = sqrt(gamma*P2/rho2);
    M2 = sqrt((1 + 0.5*(gamma - 1)*M1^2)/((gamma*M1^2) - 0.5*(gamma-1)));
    % gas velocity in shock frame
    u_2SF = c2*M2;
    % gas velocity in lab frame
    u_2LF = u_s - u_2SF;
    M2LF = u_2LF/c2;
    tc2 = sqrt(rho_liq/rho2)*(d0/u_2LF); % characteristic breakup time
end

%% file info
folder = 'data/';
filepattern = fullfile(folder, 'AIAAWG_shot3_100cpd*.csv');
files = dir(filepattern); % access names by files.name
files = natsortfiles(files);

% allocate arrays for data and processing results
nfiles = length(files);
dataMatrices = cell(1, nfiles);
%resultMatrices = cell(1, nfiles);

%% notes
% overlay VOF with a threshold using sky colormap

for q=1:nfiles 
    filename = fullfile(folder, files(q).name)
    
    % read in the data as a cell array
    dataMatrices{q} = import_schlieren_props(filename);

    % access each matrix in the cell array as 
    %test = dataMatrices{q} % where q represents each data file

    %NEED TO GENERALIZE FOR 3D
    work_array = dataMatrices{q}; % access data set
    
    % non-dim time
    if time_flag == 1
        time = work_array(1,1); % physical time
        tau = time/tc2;
    end

    %schileren
    mixrho_work = work_array(:,3); % store mixrho in work array
    mixrho_work = reshape(mixrho_work,nx+nx_stretchL+nx_stretchR,ny+2*ny_stretch); % rearrange array
    VOF_work = work_array(:,4); % store VOF in work array
    VOF_work = reshape(VOF_work,nx+nx_stretchL+nx_stretchR,ny+2*ny_stretch);

    beta = Kliq*VOF_work + Kgas*(1-VOF_work); % compute beta function
    [grad_mixrhox, grad_mixrhoy, grad_mixrhoz] = fastgrad(mixrho_work,x,y); % compute gradient
    mag_grad_mixrho = sqrt(grad_mixrhox.^2 + grad_mixrhoy.^2); % compute magnitude needs to be generalized for 3D

    phi = exp(-beta.*(mag_grad_mixrho./const)); % compute schlieren

    % generate images
    % needs to be generalized for 3D
    %% plot schlieren
    figure('Position',[500 500 1000 750]);
    ax1 = gca;
    axis equal
    h1 = contourf(XX,YY,phi,200,'linecolor','none');
    colormap(ax1,gray(256))
    C1 = colorbar(ax1); clim([0 1])
    C1.Label.Interpreter = 'latex'; C1.Label.String = '$\phi$'; C1.Label.FontSize = 16;
    axis equal
    hold on

    if time_flag == 1
        title(['Numerical Schlieren at $\tau =$ ',num2str(tau)]); xlabel('x/$d_0$'); ylabel('y/$d_0$')
    else
        title('Numerical Schlieren'); xlabel('x/$d_0$'); ylabel('y/$d_0$')
    end
    xlim([0 Lx/d0]); ylim([-Ly/(2*d0) Ly/(2*d0)])
    text(mean(xlim), 2*min(ylim), '$\phi = exp \left(-\beta(\alpha) \frac{|\nabla\rho|}{C} \right)$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center','FontSize',16);
    text(mean(xlim), 2.5*min(ylim), '$\beta(\alpha) = K_{liq}*\alpha + K_{gas} \cdot (1-\alpha)$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center','FontSize',16);
    text(mean(xlim), 2.75*min(ylim),sprintf('$K_{liq} = %.1f,\\quad K_{gas} = %.1f,\\quad C = %.1e$', Kliq, Kgas, const),'Interpreter','latex','HorizontalAlignment','center','FontSize', 16);
    ax = gca; set(ax, 'Position', [0.13 0.25 0.75 0.68]);  % [left bottom width height] adjust whitespaces

    % save image
    plot_name = fullfile(folder,sprintf('schlieren_%03d.png', q));
    saveas(gcf,plot_name);
    close(gcf);
end

function [x,y] = genmesh(Lx_ref,nx,nx_stretchL,nx_stretchR,dx_uni,ny,ny_stretch,dy_uni,alpha)

% allocate x and y array
x = zeros(nx+nx_stretchL+nx_stretchR+1,1);
y = zeros(ny+2*ny_stretch+1,1);

% generate uniform x
for i = nx_stretchL+1:nx_stretchL+1+nx
    x(i) = Lx_ref + (i-nx_stretchL-1)*dx_uni;
end

% stretch left
for i = nx_stretchL:-1:1
    dx_new = alpha*(abs(x(i+2)-x(i+1))); % new dx for stretching
    x(i) = x(i+1) - dx_new;
end

% stretch right
for i = nx_stretchL+nx+2:nx_stretchL+nx+nx_stretchR+1
    dx_new = alpha*(x(i-1)-x(i-2));
    x(i) = x(i-1) + dx_new;
end

% y mesh
y(ny/2 + ny_stretch + 1) = 0; % define centerline of domain

% generate uniform y
for j = ny/2+ny_stretch+2:ny+ny_stretch+1
    y(j) = y(j-1) + dy_uni;
end

% stretch top
for j = ny+ny_stretch+2:ny+2*ny_stretch+1
    dy_new = alpha*(y(j-2)-y(j-3));
    y(j) = y(j-1) + dy_new;
end

% mirror y mesh over centerline
for j = 1:ny/2+ny_stretch
    y(j) = -y(ny-j+(2*ny_stretch)+2);
end

end
