function [figure_with_overlay] = blueOverlay(figure_handle,XX,YY, phi, time_flag,tau, Lx, Ly, d0, Kliq, Kgas, const, VOF_work, threshold, translucent_level)
                                             
figure_with_overlay = figure_handle;
ax = axes('PositionConstraint','innerposition', 'Units', 'pixels');

% Schlieren Plot
s1  = surf(ax, XX, YY, phi, 'EdgeColor','none');  % Schlieren surface
view(ax, 2)                                       % top-down
axis(ax, 'equal', 'tight')
shading(ax, 'interp')
colormap(ax, gray(256))                           % Schlieren colormap
hold(ax, 'on')

% set limits on axes
xlim([0 Lx/d0]); ylim([-Ly/(2*d0) Ly/(2*d0)])

% Title
if time_flag == 1
    title(['Numerical Schlieren at $\tau =$ ',num2str(tau,'%.2f')]); xlabel('x/$d_0$'); ylabel('y/$d_0$')
else
    title('Numerical Schlieren'); xlabel('x/$d_0$'); ylabel('y/$d_0$')
end

% Schlieren Colorbar
cbGray = colorbar(ax,'Location','eastoutside', 'Units', 'pixels');
cbGray.AxisLocation = 'in';
cbGray.Label.Interpreter = 'latex'; cbGray.Label.String = '$\phi$'; cbGray.Label.FontSize = 16;
shift_right = 100; % pixels
pos            = cbGray.Position;         % [x y w h] of the bar
pos(1)         = pos(1) + pos(3) + shift_right;   % shift right by its own width + gap
cbGray.Position = pos;

% Use sky colormap for liquid
blueMap = flip(sky(256), 1);

% text labels
text(mean(xlim), 2*min(ylim), '$\phi = exp \left(-\beta(\alpha) \frac{|\nabla\rho|}{C} \right)$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center','FontSize',16);
text(mean(xlim), 2.5*min(ylim), '$\beta(\alpha) = K_{liq}*\alpha + K_{gas} \cdot (1-\alpha)$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center','FontSize',16);
text(mean(xlim), 2.75*min(ylim),sprintf('$K_{liq} = %.1f,\\quad K_{gas} = %.1f,\\quad C = %.1e$', Kliq, Kgas, const),'Interpreter','latex','HorizontalAlignment','center','FontSize', 16);
text(1.125*max(xlim), 1.05*min(ylim), 'Gas Phase', ... 
     'HorizontalAlignment','center', ...
     'VerticalAlignment',  'top', ...
     'Interpreter',       'latex', ...
     'FontSize',          14);
text(1.125*max(xlim), 1.25*min(ylim), 'Liquid Phase', ... 
     'HorizontalAlignment','center', ...
     'VerticalAlignment',  'top', ...
     'Interpreter',       'latex', ...
     'FontSize',          14, ...
     'Color',             blueMap(1,:));

% Overlay blue for VOF > threshold
mask = VOF_work > threshold;
blue_values   = 1.0.*mask+max(max(abs(phi)))*10; % just need to make sure these are "above" the other surface for top-down view    

% Look-up table: convert phi value value to an RGB triplet from chosen
% colormap
max_phi_liquid = max(max(phi(mask)));
blue_values_norm = phi ./ max_phi_liquid;
idx = max(1, round(blue_values_norm * (size(blueMap,1)-1)) + 1);   % indices 1…256
dropRGB        = zeros([size(VOF_work) 3]);

% populate appropriate R G B triplets
for color_index = 1:3 % [R G B]
    tmp            = dropRGB(:,:,color_index);
    tmp(mask)    = blueMap(idx(mask),color_index);  % colour only inside the disc
    dropRGB(:,:,color_index) = tmp;
end

% Per-vertex transparency: 0.5 for the droplet, 0 elsewhere
dropAlpha = translucent_level * mask;

% Overlay blue for VOF
s2 = surf(ax, XX, YY, blue_values, dropRGB, ...          % true-colour texture
          'EdgeColor', 'none', ...
          'FaceColor', 'texturemap');
s2.FaceAlpha        = 'texturemap';           % make alpha come from a map
s2.AlphaData        = dropAlpha;              % the map itself
s2.AlphaDataMapping = 'none';                 % use values literally

set( ancestor(ax,'figure'), 'Units','pixels');
set( [ax cbGray],           'Units','pixels');

% Place blue color scale (also phi) to hte right
axBlue = axes('Position', ax.Position, ...
              'Color','none',          ...
              'XTick',[], 'YTick',[],  ...
              'CLim',[0 1], 'Visible','off');

colormap(axBlue, blueMap); % Make sure color bar is consistent with plot 
axBlue.Position = ax.Position;
set(axBlue,'PositionConstraint','innerposition');

% ----- dimensions for the blue colour-bar -----
wBlue = cbGray.Position(3);   % keep same width as the grey bar
gap   = wBlue/15; 
cbBlue_Position = [cbGray.Position(1) + cbGray.Position(3) + gap,  cbGray.Position(2) ,  cbGray.Position(3) ,  cbGray.Position(4)];
cbBlue          = colorbar(axBlue,'Location', 'manual', 'units','pixels');
cbBlue.Position = cbBlue_Position;
cbBlue.Ticks=[];
set(figure_with_overlay, 'Position', [500 500 1100 750], 'Units','pixels');
