% level set and initial perturbation prototyping
clear; close all; clc

% A level set function is a scalar field where the interface is represented as a contour of the field.
% the sign of the funciton usually gives info about the locaiton of the point relative to the interface/contour
% usually initialized as a signed distance funciton

%% define grid
nx = 100;
ny = 100;
x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
[XX, YY] = meshgrid(x,y);

%% geometry definition
dctrx = 0; % x-center coordinate
dctry = 0; % y-center coordinate
dctrz = 0; % z-cetner coordinate
d0 = 1.0;  % diameter
r0 = d0/2; % radius
epsilon = d0*[0.01, 0.007, 0.005]; % amplitude of our perturbation
mode = [20, 11, 10];                 % wave number modes for our perturbation (spatial frequency of wave, 2pi/lambda) lambda=linear wavelength
phase = [pi,pi,0];             % phase shifts (offset of waves)
theta = atan2(YY-dctry,XX-dctrx);  % compute polar angle for each point (we use atan2 to keep track of the quadrant (matlab funciton) instead of just atan)

% compute total perturbation based (summation of modes)
perturbation = zeros(size(theta));
for i=1:length(mode)
    perturbation = perturbation + epsilon(i)*sin(mode(i)*theta + phase(i));
end

% define level sets (unperturbed and perturbed)
phi_unperturbed = sqrt(((XX - dctrx).^2) + (YY - dctry).^2) - r0;
phi_perturbed   = sqrt(((XX - dctrx).^2) + (YY - dctry).^2) - (r0 + perturbation);

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

% overlay the 0 level to compare
figure
contour(XX,YY,phi_unperturbed,[0 0],'k','LineWidth',2) % plot level set zero
hold on
contour(XX, YY, phi_perturbed, [0 0], 'r', 'LineWidth', 2); % perturbed interface
axis equal; 
xlabel('x'); ylabel('y'); legend('Unperturbed circle','perturbed circle')
title('Comparison of level sets')

