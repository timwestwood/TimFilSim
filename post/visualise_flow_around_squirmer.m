clear;

%% Load the data

sim_name = 'flow_field_larva_symm_squirmer';

X = load([sim_name '_flow_locations.dat']);
V = load([sim_name '_flow_velocities.dat']);
X = X(:,2:end);
V = V(:,2:end);

%% Find the flow magnitude

% Make the grid
Xbox = X(1,:);
xbox = get_grid(Xbox, 1);
ybox = get_grid(Xbox, 2);
zbox = get_grid(Xbox, 3);

box_count = zeros(numel(xbox), numel(ybox), numel(zbox));
box_sum = zeros(size(box_count));

% Bin the data
for n=1:size(X,1)
    for m=1:size(X,2)/3
        
        xtemp = X(n, 3*m - 2);
        ytemp = X(n, 3*m - 1);
        ztemp = X(n, 3*m);
        
        [~, xid] = min(abs(xtemp - xbox));
        [~, yid] = min(abs(ytemp - ybox));
        [~, zid] = min(abs(ztemp - zbox));
        
        box_count(xid,yid,zid) = box_count(xid,yid,zid) + 1;
        box_sum(xid,yid,zid) = box_sum(xid,yid,zid) + norm(V(n, 3*m - 2 : 3*m));
        
    end
end

vbox = box_sum ./ box_count;

%% Plot the magnitude

% Account for there being a singleton dimension
vbox = squeeze(vbox);

[fig_ax_dirs] = find(size(box_count)-1);
fig_x_lims = [min(Xbox(fig_ax_dirs(1):3:end)) max(Xbox(fig_ax_dirs(1):3:end))];
fig_y_lims = [min(Xbox(fig_ax_dirs(2):3:end)) max(Xbox(fig_ax_dirs(2):3:end))];

figure;
imagesc(fig_y_lims, fig_x_lims, vbox); % x and y are the wrong way round everywhere because of imagesc...
colormap jet;
axis image;

%% Plot the streamlines on top

hold on;
r = load([sim_name '_blob_references.dat']);
scatter(r(fig_ax_dirs(2):3:end), r(fig_ax_dirs(1):3:end), 'k.');

Xp = [X; NaN(1, size(X,2))];
slx = Xp(:, fig_ax_dirs(1):3:end);
sly = Xp(:, fig_ax_dirs(2):3:end);
plot(sly, slx, 'Color', [1 1 1]);

savefig(gcf, [sim_name '_flow_field_with_streamlines.fig'], 'compact');

%% Local functions

function [box] = get_grid(X, dim)

box = unique(X(dim:3:end));
box = linspace(min(box), max(box), numel(box));

end