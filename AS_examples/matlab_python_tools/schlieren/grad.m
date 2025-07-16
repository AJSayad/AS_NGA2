function [grad_phix, grad_phiy, grad_phiz] = grad(phi,x,y,z)
%% description
% numerically comput the gradient of a scalar field using centered finite
% differences. Applies one sided finite differences on the boundaries

%% input
% phi=scalar field
% dx, dy, dz = mesh spacing in each direction
%
%% output
% grad_phix,y,z = gradient of scalar field in each direction

dims = ndims(phi); % determine if 2D or 3D
[nx, ny, nz] = size(phi); 

% initialize gradient arrays 
grad_phix = zeros(size(phi));
grad_phiy = zeros(size(phi));
if dims == 3 
    grad_phiz = zeros(size(phi));
else
    grad_phiz = []; % leave empty if 2D
end

% compute the centered finite differences in the interior of the domain
if dims == 3 % if 3D
    for i = 2:nx-1
        for j = 2:ny-1
            for k = 2:nz-1
                dx = x(i+1) - x(i-1);
                dy = y(j+1) - y(j-1);
                dz = z(k+1) - z(k-1);
                grad_phix(i,j,k) = (phi(i+1,j,k) - phi(i-1,j,k)) / (2*dx);
                grad_phiy(i,j,k) = (phi(i,j+1,k) - phi(i,j-1,k)) / (2*dy);
                grad_phiz(i,j,k) = (phi(i,j,k+1) - phi(i,j,k-1)) / (2*dz);
            end
        end
    end
else % if 2D
    for i = 2:nx-1
        for j = 2:ny-1
            dx = x(i+1) - x(i-1);
            dy = y(j+1) - y(j-1);
            grad_phix(i,j,:) = (phi(i+1,j,:) - phi(i-1,j,:)) / (2*dx);
            grad_phiy(i,j,:) = (phi(i,j+1,:) - phi(i,j-1,:)) / (2*dy);
        end
    end
end % if

% handle boundary conditions with forward/backward differences
% right and left sides
for j = 1:ny
    for k = 1:nz
        dx_left = x(2) - x(1);
        dx_right = x(end) - x(end-1);
        grad_phix(1,j,k) = (phi(2,j,k) - phi(1,j,k)) / dx_left; % left side
        grad_phix(nx,j,k) = (phi(nx,j,k) - phi(nx-1,j,k)) / dx_right; % right side
    end 
end

% top and bottom (y boundaries)
for i = 1:nx
    for k = 1:nz
        dy = y(end) - y(end-1);
        grad_phiy(i,1,k) = (phi(i,2,k) - phi(i,1,k)) / dy; % bottom
        grad_phiy(i,ny,k) = (phi(i,ny,k) - phi(i,ny-1,k)) / dy; % top
    end 
end


% front and back 
if dims == 3 
    for i = 1:nx
        for j = 1:ny
            dz = z(end) - z(end-1);
            grad_phiz(i,j,1) = (phi(i,j,2) - phi(i,j,1)) / dz; % front
            grad_phiz(i,j,nz) = (phi(i,j,nz) - phi(i,j,nz-1)) / dz; % back
        end
    end
end % if

end % function