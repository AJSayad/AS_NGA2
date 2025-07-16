%% Droplet mass loss processing
clear all; close all; clc; plot_formatting; format long
time_flag = 1;
plot_mass_snaps = 1;
plot_circles = 1;
mass_file = 'mass_time_history_AIAAWG_shot3_perturbed_500cpd.mat';
plot_title = 'Normalized Shot 3 Main-Drop Mass';
fig_name = 'mass_plot_AIAAWG_shot3_500cpd';
data_file = 'AIAAWG_shot3_perturbed_500cpd_*.csv';

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

[x,y,z,nx,nx_stretch,ny,ny_stretch,nz,nz_stretch,dx_uni,dy_uni,dz_uni] = generate_mesh(d0,CPD,DX,DY,DZ,DX_stretch,DY_stretch,DZ_stretch)

x = x(1:end-1); % trim off last point (since we define cell edges)
y = y(1:end-1);
if nz>1
    [YY,XX,ZZ] = meshgrid(x,y,z)
else
    [YY,XX] = meshgrid(y,x);
end
%XX = XX/d0;
%YY = YY/d0;

nx_total = nx + nx_stretch;
ny_total = ny + 2*ny_stretch;
nz_total = nz + 2*nz_stretch;

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

if (nz>1)
    m0 = (4/3)*rho_liq*pi*((d0/2)^3); % calculate analytical initial mass
else
    m0 = rho_liq*pi*((d0/2)^2);       % calculate analytical initial mass
end

%% file info
folder = 'data/';
filepattern = fullfile(folder, data_file);
files = dir(filepattern); % access names by files.name
files = natsortfiles(files);

% allocate arrays for data and processing results
nfiles = numel(files);
dataMatrices = cell(1, nfiles);
mass_time_history = struct('time',zeros(nfiles,1),'tau',zeros(nfiles,1),'mass',zeros(nfiles,1));

if plot_mass_snaps==1
    % specify tau values for plotting mass contours
    tau_plot = [1 2 3 4 6];
    n_snap = length(tau_plot);
    small_axes = zeros(1,n_snap);
    snapshot_idxs = zeros(1,n_snap); % allocate snapshot index array
    mass_all = cell(1,nfiles); % allocate filtered mass and COM arrays for lookup when plotting
    filtered_mass_all = cell(1,nfiles);
    VOF_bin_all = cell(1,nfiles);
    filtered_VOF_bin_all = cell(1,nfiles);
    COMx_all          = zeros(1,nfiles);
    if nz>1
        vol_cell = dx_uni*dy_uni*dz_uni; % compute uniform region cell volume for COM calculation
    else
        vol_cell = dx_uni*dy_uni;
    end
end

% process data
for q=1:nfiles
    filename = fullfile(folder, files(q).name)

    % read in the data as a cell array
    dataMatrices{q} = import_liqiud_props(filename);

    %NEED TO GENERALIZE FOR 3D
    field_data_array = dataMatrices{q}; % access data set
    
    % access field data
    time = field_data_array(1,1);
    Lrho = field_data_array(:,4);
    VOF = field_data_array(:,6);

    % calc nondim time
    tau = time/tc2;

    if nz>1
        %3D
        % step 0: reshape vectors to matrices
        VOF_field  = reshape(VOF,nx_total,ny_total,nz_total);
        Lrho_field = reshape(Lrho,nx_total,ny_total,nz_total);
        % step 1: binarize VOF
        VOF_bin = zeros(nx_total*ny_total*nz_total,1);
        VOF_bin(VOF>0) = 1;
        VOF_bin_field = reshape(VOF_bin,nx_total,ny_total,nz_total);   % reshape to matrix
        % step 2: connected component labeling 
        CC_VOF_bin = bwconncomp(VOF_bin_field,26);                     % 26 connectivity in 3D 
        % Step 3: Filter out satelite droplets 
        p = regionprops3(CC_VOF_bin,"Volume");                         % calculate max area
        [~,maxidx] = max([p.Volume]);                                  % find region indices connected to max volume (number of cells)
        BWimg = transpose(cc2bw(CC_VOF_bin,ObjectsToKeep=maxidx));     % keep objects with only the max # of cells, take transpose for correct orientation
        % step 4: compute mass
        BWimgvec = reshape(BWimg,nx_total*ny_total*nz_total,1);              % reshape into vector for efficiency
        BWidx = find(BWimgvec>0);                                            % find non-zero cells (cells part of main drop)
        mass = sum(dx_uni*dy_uni*dz_uni*VOF_field.*Lrho_field.*BWimg','all');% integrate mass over filtered drop
        mass_time_history.mass(q) = mass;                                    % log mass
        mass_time_history.time(q) = time;                                    % log time
        mass_time_history.tau(q)  = tau ;                                    % log nondim time
        % save the field + COM for later snapshotting
        if plot_mass_snaps
            mass_field = dx_uni*dy_uni*dz_uni.*VOF_field.*Lrho_field; % compute full mass field for viz
            filtered_mass_field = dx_uni*dy_uni*dz_uni.*VOF_field.*Lrho_field.*BWimg';
            mass_all{q} = mass_field;
            filtered_mass_all{q} = filtered_mass_field;
            VOF_bin_all{q} = VOF_bin_field;
            filtered_VOF_bin_all{q} = BWimg';
            COMx_all(q) = sum(VOF_field.*Lrho_field.*x.*vol_cell)/sum(VOF_field.*Lrho_field.*vol_cell);
        end
    else
        % 2D
        % step 0: reshape vectors to matrices
        VOF_field  = reshape(VOF,nx_total,ny_total);
        Lrho_field = reshape(Lrho,nx_total,ny_total);
        % step 1: binarize VOF
        VOF_bin = zeros(nx_total*ny_total,1);
        VOF_bin(VOF>0) = 1;
        VOF_bin_field = reshape(VOF_bin,nx_total,ny_total);         % reshape to matrix
        % step 2: connected component labeling 
        CC_VOF_bin = bwconncomp(VOF_bin_field,8);                   % 8 connectivity in 2D 
        % Step 3: Filter out satelite droplets 
        p = regionprops(CC_VOF_bin,"Area");                         % calculate max area
        [~,maxidx] = max([p.Area]);                                 % find region indices connected to max area (number of cells)
        BWimg = transpose(cc2bw(CC_VOF_bin,ObjectsToKeep=maxidx));  % keep objects with only the max # of cells, take transpose for correct orientation
        % step 4: compute mass
        BWimgvec = reshape(BWimg,nx_total*ny_total,1);              % reshape into vector for efficiency
        BWidx = find(BWimgvec>0);                                   % find non-zero cells (cells part of main drop)
        mass = sum(dx_uni*dy_uni*VOF_field.*Lrho_field.*BWimg','all');
        mass_time_history.mass(q) = mass;                           % log mass
        mass_time_history.time(q) = time;                           % log time
        mass_time_history.tau(q)  = tau ;                           % log nondim time
        % save the field + COM for later snapshotting
        if plot_mass_snaps
            mass_field = dx_uni*dy_uni.*VOF_field.*Lrho_field; % comute full mass field for viz
            filtered_mass_field = dx_uni*dy_uni.*VOF_field.*Lrho_field.*BWimg';
            mass_all{q} = mass_field;
            filtered_mass_all{q} = filtered_mass_field;
            VOF_bin_all{q} = VOF_bin_field;
            filtered_VOF_bin_all{q} = BWimg';
            COMx_all(q) = sum(VOF_field.*Lrho_field.*x.*vol_cell)/sum(VOF_field.*Lrho_field.*vol_cell);
        end
    end
end

% save computed mass history
save(mass_file,'mass_time_history')
% to access, use data = load(mass_file)
m0_calc = mass_time_history.mass(1); % store the initial calcualted droplet mass for normalization

if (plot_mass_snaps==1)
    %% extract filtered mass
    % determin snapshot indices by matching tau
    for k = 1:n_snap
        [~, idx] = min(abs(mass_time_history.tau - tau_plot(k)));
        snapshot_idxs(k) = idx;
    end
    % now pull out those fields and COMs
    mass_snaps = mass_all(snapshot_idxs); % mass snapshots 
    filtered_mass_snaps = filtered_mass_all(snapshot_idxs); % filtered mass snapshots
    VOF_bin_snaps = VOF_bin_all(snapshot_idxs); % VOF filter for snaps
    filtered_VOF_bin_snaps = filtered_VOF_bin_all(snapshot_idxs);
    COMx               = COMx_all(snapshot_idxs); % center of mass snapshots

    %% assemble final figure
    %figure('Color','w','Units','normalized','Position',[0.1 0.1 0.8 0.8]);
    figure('Color','w','Units','pixels','Position',[100 100 1300 1000]);
    % Create a 2-row grid: top row has n_snap columns, bottom row spans all
    tfig = tiledlayout(2, n_snap, 'TileSpacing','compact', 'Padding','compact');
    for ii = 1:n_snap
        % Create axes in top row
        ax = nexttile(tfig, ii); % Top row: tile 1 to n_snap
        axes(ax); hold on; axis equal; view(2); % Set context

        % set colors [R G B]
        CO1(size(XX,1),size(YY,2),1:3) = 0;
        %CO1(:,:,1) = 1; % red
        CO1(:,:,1) = 0.610; % R
        CO1(:,:,2) = 0.762; % G
        CO1(:,:,3) = 0.812; % B
        CO1(cell2mat(VOF_bin_snaps(ii))<1) = nan; % filter based on binarized VOF
        surf(XX,YY,cell2mat(VOF_bin_snaps(ii)),CO1,'EdgeColor', 'none'); view(2); axis equal; % plot full liquid field
        CO2(size(XX,1),size(YY,2),1:3) = 0;
        %CO2(:,:,3) = 1; % blue
        CO2(:,:,1) = 0.000;
        CO2(:,:,2) = 0.000;
        CO2(:,:,3) = 0.545;
        CO2(cell2mat(filtered_VOF_bin_snaps(ii))<1) = nan; % filter based on main droplet
        surf(XX,YY,cell2mat(filtered_VOF_bin_snaps(ii))+1,CO2,'EdgeColor','none'); axis equal % plot main droplet mass

        % Set view and axis limits
        view(ax, 2); % 2D top-down view
        xlim(ax,[COMx(ii)-2*d0 COMx(ii)+2*d0]); % track droplet COM in x direction
        ylim(ax,[dctry-2*d0 dctry+2*d0]);

        % Title
        tau_val  = mass_time_history.tau(snapshot_idxs(ii));
        mass_val = mass_time_history.mass(snapshot_idxs(ii)) / m0;
        tau_str  = ['$\tau$ = ' num2str(tau_val, 1)];
        mass_str = ['$m/m_0$ = ' num2str(mass_val, 2)];
        title(ax, {tau_str, mass_str}, 'Interpreter', 'latex');
        axis off % clean up figure
        clearvars CO1 CO2
    end

    % BOTTOM ROW: big line plot
    % the bottom row is tiles (n_snap+1) through (2*n_snap)
    ax_mass = nexttile( n_snap+1, [1 n_snap] );
    plot(ax_mass, mass_time_history.tau, mass_time_history.mass ./ mass_time_history.mass(1), 'o--','Color',rgb('cornflowerblue'),'LineWidth',1.5);
    xlim([0 max(mass_time_history.tau)+0.25]); ylim([0 1.1]);
    xlabel(ax_mass,'$\tau$');
    ylabel(ax_mass,'$m/m_0$');
    xlim(ax_mass,[0 5.3]); ylim(ax_mass,[0 1]); 
    title(ax_mass,plot_title);
    grid(ax_mass,'on');

    if plot_circles==1
        ax_mass = nexttile( n_snap+1, [1 n_snap] );
        axes(ax_mass)
        hold(ax_mass,'on')
        tau_snaps = mass_time_history.tau( snapshot_idxs );
        m_snaps = mass_time_history.mass(snapshot_idxs) ./ mass_time_history.mass(1);
        plot(ax_mass, tau_snaps, m_snaps, 'o','color',rgb('Navy'),'MarkerFaceColor',rgb('Navy'),'MarkerSize',8, 'LineWidth',1.5);
    end
    exportgraphics(gcf, fig_name+".png", 'ContentType', 'image');
    exportgraphics(gcf, fig_name+".eps", 'ContentType', 'vector');
    savefig(gcf, fig_name);
else
    figure()
    plot(mass_time_history.tau,mass_time_history.mass/m0_calc,'o--','Color',rgb('cornflowerblue'))
    title(plot_title); xlabel('$\tau$'); ylabel('Normalized Mass m/m$_0$')
    legend('100cpd')
    exportgraphics(gcf, fig_name+".eps", 'ContentType', 'vector');
    savefig(gcf, fig_name);
end

% save mass snapshots and filtered mass snapshots as a file
save('tau_snaps',tau_plot) % save tau snapshots
save('VOF_bin_snaps',VOF_bin_snaps) % save binarized VOF snapshots
save('filtered_VOF_bin_snaps',filtered_VOF_bin_snaps) % save filtered binarized VOF snapshots
