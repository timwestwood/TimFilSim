clear;

%% Specify which simulations to use

scale_by_sync_array = false;

% The parameters aren't automatically stored anywhere, so we must
% hard-code the pairings into the script.
tags_and_wavelengths = {

'0.5L', 0.5;... % The order given here is irrelevant.
'1.0L', 1.0;...
'1.5L', 1.5;...
'2.0L', 2.0;...

};

tags_and_spacings = {

'0.1L', 0.1;... % The order given here is irrelevant too.
% '0.5L', 0.5;...
% '1.0L', 1.0;...
% '1.5L', 1.5;...
% '0.1053L', 0.1053;...
% '0.25L', 0.25;...
% '0.75L', 0.75;...
% '0.375L', 0.375;...
% '0.6L', 0.6;...
'0.11266L', 0.11266;...
'0.129L', 0.129;...
'0.2023L', 0.2023;...
'0.1508L', 0.1508;...

};

name_trunk = 'review_wavelength_and_spacing_fixed_area_sweep_'; % Should include a trailing underscore.
top_dir = '~/TimFilSim/'; % Search for data in this directory and all its subdirectories.

%% Loop over the simulations

% Sort the wavelengths and spacings.
[wavelengths, I] = sort(cell2mat(tags_and_wavelengths(:,2)));
tags_and_wavelengths = tags_and_wavelengths(I,:);

[spacings, I] = sort(cell2mat(tags_and_spacings(:,2)));
tags_and_spacings = tags_and_spacings(I,:);

% Allocate the array to plot. Each column is a given spacing and each row a
% given wavelength; i.e. we will plot spacing on the horizontal axis and
% wavelength on the vertical.
time_avg_R = NaN(size(tags_and_wavelengths,1), size(tags_and_spacings,1));
time_avg_Q = time_avg_R;

for row = 1:size(tags_and_wavelengths,1)

    name = [name_trunk tags_and_wavelengths{row,1} '_and_'];

    for col = 1:size(tags_and_spacings,1)

        full_name = [name tags_and_spacings{col,1}];

        [Rbar, Qbar, L, visc] = read_sim_data(top_dir, full_name);

        time_avg_R(row,col) = Rbar * L^3 / visc; % Provide the factors needed to end up with a dimensionless efficiency.
        time_avg_Q(row,col) = Qbar;

    end

end

pumping_eff = (time_avg_Q.^2)./time_avg_R;

if scale_by_sync_array

    name = [name_trunk 'inf_and_'];

    sync_array_avg_R = NaN(1, size(tags_and_spacings,1));
    sync_array_avg_Q = sync_array_avg_R;

    for col = 1:size(tags_and_spacings,1)

        full_name = [name tags_and_spacings{col,1}];

        [Rbar, Qbar, L, visc] = read_sim_data(top_dir, full_name);

        sync_array_avg_R(col) = Rbar * L^3 / visc; % Provide the factors needed to end up with a dimensionless efficiency.
        sync_array_avg_Q(col) = Qbar;

    end

    sync_array_eff = (sync_array_avg_Q.^2)./sync_array_avg_R;

    pumping_eff = pumping_eff ./ sync_array_eff;

end

contourf(spacings, wavelengths, pumping_eff, 'EdgeColor', 'none');
axis equal;
colormap jet;
bar_handle = colorbar;

if scale_by_sync_array

    colourbar_label = '$\varepsilon\left(\Delta_\theta, \lambda\right)/\varepsilon_\mathrm{sync}\left(\Delta_\theta\right)$';

else

    colourbar_label = '$\varepsilon\left(\Delta_\theta, \lambda\right)$';

end
bar_handle.Label.String = colourbar_label;
bar_handle.Label.Interpreter = 'latex';
bar_handle.Label.FontSize = 24;

box on;
set(gca, 'FontName', 'Times', 'FontSize', 24);
xlabel('$\Delta_\theta/L$', 'Interpreter', 'latex');
ylabel('$\lambda/L$', 'Interpreter', 'latex');

if scale_by_sync_array

    fig_name = fullfile(top_dir, [name_trunk 'scaled_efficiency.fig']);

else

    fig_name = fullfile(top_dir, [name_trunk 'efficiency.fig']);

end

savefig(gcf, fig_name, 'compact');

%% Local functions

function [Rbar, Qbar, L, visc] = read_sim_data(top_dir, full_name)

sim = dir(fullfile(top_dir, '/**/', [full_name '.par']));

assert(numel(sim)==1, ['There are multiple candidate simulations for ' full_name]);

fp = fullfile(sim.folder, sim.name);
fp = fp(1:end-4);

par = read_parameter_file(fp);

if exist([fp '_dissipation_rate.dat'], 'file')

    R = load([fp '_dissipation_rate.dat']);

else

    Fcomplete = load([fp '_seg_forces.dat']);
    F = Fcomplete(:,2:end);

    V = load([fp '_seg_vels.dat']);
    V = V(:,2:end);

    R = sum(V.*F, 2)/par.NFIL;

    dlmwrite([fp '_dissipation_rate.dat'], R, 'precision', '%.10e');

end

if exist([fp '_flow_rate.dat'], 'file')

    Q = load([fp '_flow_rate.dat']);

else

    if ~exist('Fcomplete', 'var')

        Fcomplete = load([fp '_seg_forces.dat']);

    end

    F = Fcomplete(:,3:6:end); % Only want component in flow (i.e. y) direction.

    Z = load([fp '_seg_states.dat']);
    Z = Z(:,4:3:end); % Only want z-components (i.e. heights above the wall).

    Q = sum(F.*Z, 2)/(pi * par.MU * par.NFIL);

    dlmwrite([fp '_flow_rate.dat'], Q, 'precision', '%.10e');

end

Rbar = mean(R);
Qbar = mean(Q);
L = par.FIL_LENGTH;
visc = par.MU;

end

