%% Post processing example for matlab
clear all; close all; clc; plot_formatting; format long

% parallelization
% data distribution is handled intrinsically by matlab client using a round robin distribution
poolobj = gcp('nocreate'); % get current pool object
if ~isempty(poolobj) % if pool object is not empty, delete it so we can create a new one
    delete(poolobj)
end

np = 4; % number of tasks for processing data
parpool('local',np) % initialize parallel pool
set(0, 'DefaultFigureVisible', 'off'); % render graphics off-screen

time_flag = 1; % compute non-dimensional time?
%% mesh grid
d0 = 0.00187;    % droplet diameter
dctrx = 0.0028;  % droplet location in x
dctry = 0.0;     % droplet location in y
dctrz = 0.0;     % droplet location in z
CPD = 500;       % resolution (cells/d0)
DX = 12;         % # of droplet diameter lengths in x direction
DY = 8;          % # of droplet diameter lengths in y direction (total domain height)
DZ = 0;          % # of droplet diameter lengths in z direction (total domain depth)

DX_stretch = 2; % number of droplet diameter lengths to add at the end of the x direction
DY_stretch = 2; % number of droplet diameter lengths to add at the top AND bottom of y direction
DZ_stretch = 0; % number of droplet diameter lengths to add at the top AND bottom of z direction

[x,y,z,nx,nx_stretch,ny,ny_stretch,nz,nz_stretch,dx_uni,dy_uni,dz_uni] = generate_mesh(d0,CPD,DX,DY,DZ,DX_stretch,DY_stretch,DZ_stretch);

x = x(1:end-1); % trim off last point (since we define cell edges)
y = y(1:end-1);
if nz>1
    [YY,XX,ZZ] = meshgrid(x,y,z)
else
    [YY,XX] = meshgrid(y,x);
end
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
filepattern = fullfile(folder, 'schlieren_data_*.csv');
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
    dataMatrices{q} = importschlierenvars(filename);

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
    mixrho_work = reshape(mixrho_work,nx+nx_stretch,ny+2*ny_stretch); % rearrange array
    mixrho_work = distributed(mixrho_work); % AS: distribute mixrho work array among workers
    VOF_work = work_array(:,4); % store VOF in work array
    VOF_work = reshape(VOF_work,nx+nx_stretch,ny+2*ny_stretch);
    VOF_work = distributed(VOF_work); % AS: distribute VOF work array among tasks
    spmd % single program multiple data
        beta = Kliq*VOF_work + Kgas*(1-VOF_work); % compute beta function on multiple tasks
        [grad_mixrhox, grad_mixrhoy, grad_mixrhoz] = fastgrad(mixrho_work,x,y); % compute gradient on multiple tasks
        mag_grad_mixrho = sqrt(grad_mixrhox.^2 + grad_mixrhoy.^2); % compute magnitude needs to be generalized for 3D
        phi = exp(-beta.*(mag_grad_mixrho./const)); % compute schlieren
    end

    % communicate arrays to matlab client
    phi = gather(phi);
    VOF_work = gather(VOF_work);

    % generate images
    % needs to be generalized for 3D
    %% plot schlieren

    % JR updated here
    figure_handle = figure('Position',[500 500 1000 750], 'Units','pixels');
    threshold = 0.95; % VOF threshold for blue overlay
    translucent_level = 0.9; % setting facealpha value = 1, blue is opaque, 0 blue is fully translucent
    [new_figure] = blueOverlay(figure_handle,XX,YY, phi, time_flag,tau, DX, DY, d0, Kliq, Kgas, const, VOF_work, threshold, translucent_level);

    % save image
    plot_name = fullfile(folder,sprintf('schlieren_%03d.png', q));
    %saveas(new_figure,plot_name)
    exportgraphics(new_figure,plot_name)

    close(new_figure);
    % JR ended here
end
