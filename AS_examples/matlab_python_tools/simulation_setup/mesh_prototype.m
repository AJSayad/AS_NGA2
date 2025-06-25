% mesh generation prototyping
clear; close all; clc

%% input file
d0 = 1;    % droplet diameter
CPD = 100; % resolution (cells/d0)
DX = 2;    % # of droplet diameter lengths in x direction
DY = 2;    % # of droplet diameter lengths in y direction (total domain height)
DZ = 2;    % # of droplet diameter lengths in z direction (total domain depth)

DX_stretch = 1; % number of droplet diameter lengths to add at the end of the x direction
DY_stretch = 1; % number of droplet diameter lengths to add at the top AND bottom of y direction
DZ_stretch = 1; % number of droplet diameter lengths to add at the top AND bottom of z direction

%% domain setup and mesh generation
Lx = DX*d0; % compute domain length in x
Ly = DY*d0; % compute domain length in y
Lz = DZ*d0; % compute domain length in z

Lx_stretch = DX_stretch*d0; % compute length of stretched region in x
Ly_stretch = DY_stretch*d0; % compute length of stretched region in y
Lz_stretch = DZ_stretch*d0; % compute length of stretched region in z

nx = ceil((CPD*Lx)/d0); % compute # of uniform region cells in x
ny = ceil((CPD*Ly)/d0); % compute # of uniform region cells in y
nz = ceil((CPD*Lz)/d0); % compute # of uniform region cells in z

dx_uni = Lx/nx; % compute uniform mesh spacing in x 
dy_uni = Ly/ny; % compute uniform mesh spacing in y
dz_uni = Lz/nz; % compute uniform mesh spacing in z

if DZ==0 % if 2D
    nz = 1;
    Lz = dx_uni;
    dz_uni = dx_uni;
end

alpha = 1.03; % stretching ratio

% !> compute nx_stretch
nx_stretch_max = 1000; nx_stretch = 1; dx_test = dx_uni; Lx_test = 0;
while nx_stretch < nx_stretch_max
    dx_test = dx_uni*alpha^nx_stretch;
    Lx_test = Lx_test + dx_test;
    nx_stretch = nx_stretch + 1;
    if Lx_test >= Lx_stretch
        break % in fortran this will be exit
    end
end
x = zeros(nx+nx_stretch+1,1); % allocate x mesh

% !> generate uniform mesh cell edges in x
for i=1:nx+1
    x(i) = i*dx_uni;
end
% !> generate stretched mesh in x
for i=nx+2:nx+nx_stretch+1
    dx_stretch = alpha*(x(i-1) - x(i-2));
    x(i) = x(i-1) + dx_stretch;
end

% !> compute ny_stretch
ny_stretch_max = 1000; ny_stretch = 1; dy_test = dy_uni; Ly_test = 0;
while ny_stretch < ny_stretch_max
    dy_test = dy_uni*alpha^ny_stretch;
    Ly_test = Ly_test + dy_test;
    ny_stretch = ny_stretch + 1;
    if Ly_test >= Ly_stretch
        break % in fortran this will be exit
    end
end
if (mod(ny_stretch,2)~=0)
    ny_stretch = ny_stretch +1; % add 1 to ny_stretch to ensure divisibility by 2 for mirroring
end
y = zeros(ny+2*ny_stretch+1,1); % allocate y array
y(ny/2 + ny_stretch + 1) = 0;   % define centerline to be zero
if (mod(ny+2*ny_stretch,2)~=0)
    ny = ny +1;
    warning('[ShockDrop Grid init] Number of cells in y must be divisible by 2 for grid mirroring. Adding 1 cell to ny.')
end

% !> uniform mesh cell edges in y 
for j=(ny/2 + ny_stretch + 2):ny+ny_stretch+1
    y(j) = y(j-1) + dy_uni;
end

% !> stretched region in y
for j=ny+ny_stretch+2:ny+(2*ny_stretch)+1
    dy_stretch = alpha*(y(j-1)-y(j-2));
    y(j) = y(j-1) + dy_stretch;
end

% !> mirror y across y=0 centerline 
for j=1:ny/2+ny_stretch
    y(j) = -y(ny-j+(2*ny_stretch)+2);
end

if nz>1
    % !> compute nz_stretch
    nz_stretch_max = 1000; nz_stretch = 1; dz_test = dz_uni; Lz_test = 0;
    while nz_stretch < nz_stretch_max
        dz_test = dz_uni*alpha^nz_stretch;
        Lz_test = Lz_test + dz_test;
        nz_stretch = nz_stretch + 1;
        if Lz_test >= Lz_stretch
            break % in fortran this will be exit
        end
    end
    if (mod(nz_stretch,2)~=0)
        nz_stretch = nz_stretch +1; % add 1 to nz_stretch to ensure divisibilitz by 2 for mirroring
    end
    z = zeros(nz+2*nz_stretch+1,1); % allocate z array
    z(nz/2 + nz_stretch + 1) = 0;   % define centerline to be zero
    if (mod(nz+2*nz_stretch,2)~=0)
        nz = nz +1;
        warning('[ShockDrop Grid init] Number of cells in z must be divisible by 2 for grid mirroring. Adding 1 cell to nz.')
    end

    % !> uniform mesh cell edges in z
    for k=(nz/2 + nz_stretch + 2):nz+nz_stretch+1
        z(k) = z(k-1) + dz_uni;
    end

    % !> stretched region in z
    for k=nz+nz_stretch+2:nz+(2*nz_stretch)+1
        dz_stretch = alpha*(z(k-1)-z(k-2));
        z(k) = z(k-1) + dz_stretch;
    end

    % !> mirror z across z=0 centerline
    for k=1:nz/2+nz_stretch
        z(k) = -z(nz-k+(2*nz_stretch)+2);
    end
else
    z = []; % if 2D leave z empty
end

%% test plotting 
[XX,YY,ZZ] = meshgrid(y,x,z);
figure('Color','w'); tiledlayout(1,3); title('Generated mesh pints (cell edges)')
nexttile; plot(x,'ro--'); title('x'); xlabel('Index'); ylabel('Physical location')
nexttile; plot(y,'bo--'); title('y'); xlabel('Index'); ylabel('Physical location')
nexttile; plot(z,'go--'); title('z'); xlabel('Index'); ylabel('Physical location')
