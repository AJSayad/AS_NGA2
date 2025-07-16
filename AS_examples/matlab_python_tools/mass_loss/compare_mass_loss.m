%% combined mass loss plot
clear; close all; clc; plot_formatting

%% TO DO: 1) figure out why main drop isn't plotting 2) update plotting colors?

plot_title = 'Unperturbed drop vs. Perturbed drop mass loss';

%% step 0: generate mesh
% we only need to do this once for this comparison since both sims were run on the same mesh
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

[x,y,z,nx,nx_stretch,ny,ny_stretch,nz,nz_stretch,dx_uni,dy_uni,dz_uni] = generate_mesh(d0,CPD,DX,DY,DZ,DX_stretch,DY_stretch,DZ_stretch)
x = x(1:end-1);          % trim off last point (since we define cell edges)
y = y(1:end-1);
[YY,XX] = meshgrid(y,x); % create mesh grid

%% step 1: load data
% perturbed drop data
perturbed_tau_snaps         = cell2mat(struct2cell(load('perturbed_data/tau_snaps.mat')));                % load tau snapshots
perturbed_COMx              = cell2mat(struct2cell(load('perturbed_data/COMx.mat')));                     % load COMx snapshots
perturbed_mass_loss_history = cell2mat(struct2cell(load('perturbed_data/mass_time_history_AIAAWG_shot3_perturbed_500cpd.mat'))); % load mass loss time history data 
perturbed_VOFbin            = struct2cell(load('perturbed_data/VOF_bin_snaps.mat'));                      % load binarized full VOF field data
perturbed_filtered_VOFbin   = struct2cell(load('perturbed_data/filtered_VOF_bin_snaps.mat'));             % load binarized filtered VOF field data
nsnaps_perturbed = numel(perturbed_VOFbin{1,1});                                           % get the number of snapshots 

% unperturbed drop data
unperturbed_tau_snaps         = cell2mat(struct2cell(load('unperturbed_data/tau_snaps.mat')));                % load tau snapshots
unperturbed_COMx              = cell2mat(struct2cell(load('unperturbed_data/COMx.mat')));                     % load COMx snapshots
unperturbed_mass_loss_history = cell2mat(struct2cell(load('unperturbed_data/mass_time_history_AIAAWG_shot3_unperturbed_500cpd.mat'))); % load mass loss time history data 
unperturbed_VOFbin            = struct2cell(load('unperturbed_data/VOF_bin_snaps.mat'));                      % load binarized full VOF field data
unperturbed_filtered_VOFbin   = struct2cell(load('unperturbed_data/filtered_VOF_bin_snaps.mat'));             % load binarized filtered VOF field data
nsnaps_unperturbed = numel(unperturbed_VOFbin{1,1});                                           % get the number of snapshots 


%% step 2: assemble figure
figure('Color','w','Units','pixels','Position',[100 100 1300 1000]);
% Create a 2-row grid: top row has n_snap columns, bottom row spans all
tfig = tiledlayout(3, nsnaps_perturbed, 'TileSpacing','compact', 'Padding','compact'); % change this to 3 rows eventually

unperturbed_m_snaps = zeros(1,nsnaps_unperturbed);
% MIDDLE ROW: unperturbed drop
for ii = nsnaps_unperturbed+1:2*nsnaps_unperturbed
    % Create axes in top row
    ax = nexttile(tfig, ii); % Top row: tile 1 to n_snap
    axes(ax); hold on; axis equal; view(2); % Set context

    % set colors [R G B]
    CO1(size(XX,1),size(YY,2),1:3) = 0;
    %CO1(:,:,1) = 1; % red
    CO1(:,:,1) = 0.610; % R
    CO1(:,:,2) = 0.762; % G
    CO1(:,:,3) = 0.812; % B
    CO1(cell2mat(unperturbed_VOFbin{1,1}(ii-nsnaps_unperturbed))<1) = nan; % filter based on binarized VOF
    surf(XX,YY,cell2mat(unperturbed_VOFbin{1,1}(ii-nsnaps_unperturbed)),CO1,'EdgeColor', 'none'); view(2); axis equal; % plot full liquid field
    CO2(size(XX,1),size(YY,2),1:3) = 0;
    %CO2(:,:,3) = 1; % blue
    CO2(:,:,1) = 0.000;
    CO2(:,:,2) = 0.000;
    CO2(:,:,3) = 0.545;
    CO2(cell2mat(unperturbed_filtered_VOFbin{1,1}(ii-nsnaps_unperturbed))<1) = nan; % filter based on main droplet
    surf(XX,YY,cell2mat(unperturbed_filtered_VOFbin{1,1}(ii-nsnaps_unperturbed))+1,CO2,'EdgeColor','none'); axis equal % plot main droplet mass

    % Set view and axis limits
    view(ax, 2); % 2D top-down view
    xlim(ax,[unperturbed_COMx(ii-nsnaps_unperturbed)-2*d0 unperturbed_COMx(ii-nsnaps_unperturbed)+2*d0]); % track droplet COM in x direction
    ylim(ax,[dctry-2*d0 dctry+2*d0]);

    % Title
    tau_val  = unperturbed_tau_snaps(ii-nsnaps_unperturbed);
    % we need to find the indices closest to tau_val to get the mass value
    temp_tau = abs(unperturbed_mass_loss_history.tau - tau_val);
    idx = find(temp_tau==min(abs(unperturbed_mass_loss_history.tau - tau_val)));
    mass_val = unperturbed_mass_loss_history.mass(idx) / unperturbed_mass_loss_history.mass(1);
    unperturbed_m_snaps(ii-nsnaps_unperturbed) = mass_val; % store mass value for plotting 
    tau_str  = ['$\tau$ = ' num2str(tau_val, 2)];
    mass_str = ['$m/m_0$ = ' num2str(mass_val, 2)];
    title(ax, {tau_str, mass_str}, 'Interpreter', 'latex');
    %axis off % clean up figure
    if ii == nsnaps_unperturbed+1 
        ylabel(ax, 'Unperturbed')
    end
    set(ax, 'XTick', [], 'YTick', []) % Hide ticks but keep label and box
    clearvars CO1 CO2
end

perturbed_m_snaps = zeros(1,nsnaps_perturbed);
% TOP ROW: perturbed drop
for ii = 1:nsnaps_perturbed
    % Create axes in top row
    ax = nexttile(tfig, ii); % Top row: tile 1 to n_snap
    axes(ax); hold on; axis equal; view(2); % Set context

    % set colors [R G B]
    CO1(size(XX,1),size(YY,2),1:3) = 0;
    %CO1(:,:,1) = 1; % red
    CO1(:,:,1) = 0.610; % R
    CO1(:,:,2) = 0.762; % G
    CO1(:,:,3) = 0.812; % B
    CO1(cell2mat(perturbed_VOFbin{1,1}(ii))<1) = nan; % filter based on binarized VOF
    surf(XX,YY,cell2mat(perturbed_VOFbin{1,1}(ii)),CO1,'EdgeColor', 'none'); view(2); axis equal; % plot full liquid field
    CO2(size(XX,1),size(YY,2),1:3) = 0;
    %CO2(:,:,3) = 1; % blue
    CO2(:,:,1) = 0.000;
    CO2(:,:,2) = 0.000;
    CO2(:,:,3) = 0.545;
    CO2(cell2mat(perturbed_filtered_VOFbin{1,1}(ii))<1) = nan; % filter based on main droplet
    surf(XX,YY,cell2mat(perturbed_filtered_VOFbin{1,1}(ii))+1,CO2,'EdgeColor','none'); axis equal % plot main droplet mass

    % Set view and axis limits
    view(ax, 2); % 2D top-down view
    xlim(ax,[perturbed_COMx(ii)-2*d0 perturbed_COMx(ii)+2*d0]); % track droplet COM in x direction
    ylim(ax,[dctry-2*d0 dctry+2*d0]);

    % Title
    tau_val  = perturbed_tau_snaps(ii);
    % we need to find the indices closest to tau_val to get the mass value
    temp_tau = abs(perturbed_mass_loss_history.tau - tau_val);
    idx = find(temp_tau==min(abs(perturbed_mass_loss_history.tau - tau_val)));
    mass_val = perturbed_mass_loss_history.mass(idx) / perturbed_mass_loss_history.mass(1);
    perturbed_m_snaps(ii) = mass_val; % store mass value for plotting 
    tau_str  = ['$\tau$ = ' num2str(tau_val, 2)];
    mass_str = ['$m/m_0$ = ' num2str(mass_val, 2)];
    title(ax, {tau_str, mass_str}, 'Interpreter', 'latex');
    if ii == 1
        ylabel(ax, 'Perturbed')
    end
    set(ax, 'XTick', [], 'YTick', []) % Hide ticks but keep label
    clearvars CO1 CO2
end

% BOTTOM ROW: big line plot
% the bottom row is tiles (n_snap+1) through (2*n_snap)
ax_mass = nexttile(tfig,2*nsnaps_perturbed+1, [1 nsnaps_perturbed] );
plot(ax_mass, unperturbed_mass_loss_history.tau, unperturbed_mass_loss_history.mass ./ unperturbed_mass_loss_history.mass(1), 'o--','Color','k','LineWidth',1.5);
hold on
plot(ax_mass, perturbed_mass_loss_history.tau, perturbed_mass_loss_history.mass ./ perturbed_mass_loss_history.mass(1), 'x--','Color',rgb('cornflowerblue'),'LineWidth',1.5);
xlim([0 max(unperturbed_mass_loss_history.tau)+0.15]); ylim([0 1.1]);
xlabel(ax_mass,'$\tau$');
ylabel(ax_mass,'$m/m_0$');
%title(ax_mass,plot_title);
grid(ax_mass,'on');
axes(ax_mass)
hold(ax_mass,'on')
plot(ax_mass, unperturbed_tau_snaps, unperturbed_m_snaps, 'o','color',rgb('Navy'),'MarkerFaceColor',rgb('Navy'),'MarkerSize',8, 'LineWidth',1.5);
plot(ax_mass, perturbed_tau_snaps, perturbed_m_snaps, 'x','color',rgb('Navy'),'MarkerFaceColor',rgb('Navy'),'MarkerSize',8, 'LineWidth',1.5);
legend('unperturbed','perturbed','','')
