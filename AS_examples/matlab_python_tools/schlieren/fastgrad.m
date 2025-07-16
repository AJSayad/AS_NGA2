function [gx, gy, gz] = fastgrad(phi, x, y, z)
%GRAD   Vectorized centered finite-difference gradient with one-sided edges.
%   [GX, GY, GZ] = GRAD(PHI, X, Y, Z) returns components of ∇PHI on a
%   possibly nonuniform grid defined by vectors X, Y, Z.
%
%   Works for 2D (size(phi) = [nx, ny]) or 3D ([nx, ny, nz]). For 2D, GZ
%   is returned empty.

% Determine dimensions
sz = size(phi);
dims = numel(sz);
nx = sz(1); ny = sz(2);
if dims == 3, nz = sz(3); else nz = 1; end

% Preallocate
gx = zeros(sz); 
gy = zeros(sz);
if dims==3, gz = zeros(sz); else gz = []; end

%% Interior centered differences
% --- X-direction ---
% compute delta x for each i in 2..nx-1
dx = x(3:end) - x(1:end-2);             % length nx-2
% use broadcasting: reshape to [nx-2,1,1] (or [nx-2,1] for 2D)
Dx = reshape(2*dx, [numel(dx), ones(1,dims-1)]);
% slice phi at i+1 and i-1
gx(2:end-1,:,:) = (phi(3:end,:,:) - phi(1:end-2,:,:)) ./ Dx;

% --- Y-direction ---
dy = y(3:end) - y(1:end-2);
Dy = reshape(2*dy, [1, numel(dy), ones(1,dims-2)]);
gy(:,2:end-1,:) = (phi(:,3:end,:) - phi(:,1:end-2,:)) ./ Dy;

% --- Z-direction (3D only) ---
if dims==3
    dz = z(3:end) - z(1:end-2);
    Dz = reshape(2*dz, [1,1,numel(dz)]);
    gz(:,:,2:end-1) = (phi(:,:,3:end) - phi(:,:,1:end-2)) ./ Dz;
end

%% Boundaries: one-sided
% X boundaries for all j,k
dxL = x(2) - x(1);
dxR = x(end) - x(end-1);
gx(1,:,:)   = (phi(2,:,:)   - phi(1,:,:))   / dxL;
gx(end,:,:) = (phi(end,:,:)- phi(end-1,:,:)) / dxR;

% Y boundaries for all i,k
dyB = y(2) - y(1);
dyT = y(end) - y(end-1);
gy(:,1,:)   = (phi(:,2,:)   - phi(:,1,:))   / dyB;
gy(:,end,:) = (phi(:,end,:) - phi(:,end-1,:)) / dyT;

% Z boundaries for all i,j (3D only)
if dims==3
    dzF = z(2) - z(1);
    dzBa= z(end)- z(end-1);
    gz(:,:,1)   = (phi(:,:,2)   - phi(:,:,1))   / dzF;
    gz(:,:,end) = (phi(:,:,end) - phi(:,:,end-1)) / dzBa;
end

end
