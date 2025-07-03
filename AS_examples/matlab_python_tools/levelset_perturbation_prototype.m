% level set and initial perturbation prototyping
% AS 7/3/2025

% A level set function is a scalar field where the interface is represented as a contour of the field.
% the sign of the funciton usually gives info about the locaiton of the point relative to the interface/contour
% usually initialized as a signed distance funciton

% be careful when defining wavelength, if it is smaller than the mesh spacing, the wave will not be resolved

%% define grid
dctrx = 0;    % x-center coordinate
dctry = 0;    % y-center coordinate
dctrz = 0;    % z-center coordinate
d0 = 0.00187; % diameter
r0 = d0/2;    % radius
nx = 200;
ny = 200;
x = linspace(-1.25*d0,1.25*d0,nx);
y = linspace(-1.25*d0,1.25*d0,ny);
[XX, YY] = meshgrid(x,y);

%% geometry definition
lambda = (1e-6)*[280, 200, 150, 100, 50];              % wavelengths in meters
epsilon = lambda.*[0.025, 0.025, 0.025, 0.025, 0.025]; % amplitude (2.5% of wavelength)
theta = atan2(YY-dctry,XX-dctrx);                      % compute polar angle for each point
r = sqrt((XX - dctrx).^2 + (YY - dctry).^2);           % radial distance

% compute total perturbation
perturbation = zeros(size(theta));
for i = 1:length(lambda)
    perturbation = perturbation + epsilon(i)*sin(2*pi*r0*theta/lambda(i));
end

% define level sets (unperturbed and perturbed)
phi_unperturbed = r - r0;
phi_perturbed   = r - (r0 + perturbation);
max_perturbation = 100*(max(perturbation(:))/d0) % check max perturbation

% overlay the 0 level to compare
figure
contour(XX,YY,phi_unperturbed,[0 0],'k','LineWidth',2) % plot level set zero
hold on
contour(XX, YY, phi_perturbed, [0 0], 'r', 'LineWidth', 2); % perturbed interface
axis equal; grid on;
xlabel('x'); ylabel('y'); legend('Unperturbed circle','perturbed circle')
title('Comparison of level sets')


% % plot unperturbed level set
% figure;
% contourf(XX,YY,phi_unperturbed,50,'linecolor','none')  % plot filled contours
% hold on
% contour(XX,YY,phi_unperturbed,[0 0],'k','LineWidth',2) % plot level set zero
% axis equal; colorbar; colormap('turbo');
% title('2D level set of an unperturbed circle'); xlabel('x'); ylabel('y')
% 
% % plot perturbed level set
% figure;
% contourf(XX, YY, phi_perturbed, 50, 'LineColor', 'none');   % plot filled contours
% hold on;
% contour(XX, YY, phi_perturbed, [0 0], 'k', 'LineWidth', 2); % plot level set zero
% axis equal; colorbar; colormap('turbo')
% title('2D Multi-Mode Perturbed Circle (Level Set)'); xlabel('x'); ylabel('y')
