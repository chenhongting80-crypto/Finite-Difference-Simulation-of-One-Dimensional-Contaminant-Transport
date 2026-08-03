function results = run_method1_distributed(output_root)
%RUN_METHOD1_DISTRIBUTED Simulate distributed biostimulation/bioaugmentation.
%
% The model contains a 35 m source zone followed by a 30 m reactive zone.
% It uses central finite differences for advection and dispersion and
% Monod-type sequential reaction terms for TCE -> cDCE -> VC.
%
% The implementation uses two rolling time layers and records only selected
% snapshots and time series, reducing memory use relative to the archived
% course script.

if nargin < 1 || strlength(string(output_root)) == 0
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    output_root = fullfile(repo_root, 'results', 'generated');
end
output_dir = fullfile(output_root, 'method1');
ensure_directory(output_dir);

%% Active executable parameters
p.alpha_L = 0.015;            % longitudinal dispersivity, m
p.tortuosity = 0.8;           % dimensionless
p.K = 17;                     % hydraulic conductivity, m/day
p.porosity = 0.38;            % dimensionless
p.gradient = 0.0033;          % hydraulic gradient
p.velocity = p.K / p.porosity * p.gradient; % pore-water velocity, m/day

p.D0_TCE = 10.1e-6 * 3600 * 24 * 1e-4;   % m^2/day
p.D0_cDCE = 11.4e-6 * 3600 * 24 * 1e-4;  % m^2/day
p.D0_VC = 13.3e-6 * 3600 * 24 * 1e-4;    % m^2/day
p.Dx_TCE = p.velocity * p.alpha_L + p.porosity * p.tortuosity * p.D0_TCE;
p.Dx_cDCE = p.velocity * p.alpha_L + p.porosity * p.tortuosity * p.D0_cDCE;
p.Dx_VC = p.velocity * p.alpha_L + p.porosity * p.tortuosity * p.D0_VC;

p.dx = 5 / 167;               % m; uniform grid spacing that exactly represents all zone endpoints
p.dt = 0.03;                  % day
p.L1 = 35;                    % source-zone length, m
p.L2 = 30;                    % reactive-zone length, m
p.t_end = 700;                % day
p.KT_A_V = 7.8e-4;            % day^-1
p.TCE_solubility = 1300;      % mg/L
p.initial_TCE = 10;           % mg/L
p.initial_source_mass = 1e9;  % mg

p.mu_TCE = 3.42;              % day^-1
p.mu_cDCE = 0.16;             % day^-1
p.mu_VC = 0.52;               % day^-1
p.cell_density = 2e9;         % cells/L (submitted script value)
p.Ks_TCE = 10 * 131.38;       % mg/L
p.Ks_cDCE = 3.3 * 96.94;      % mg/L
p.Ks_VC = 2.6 * 62.5;         % mg/L
p.yield = 5.2e8 / 35.45;      % cells/mg
p.cross_section = 10.67 * 35; % m^2
p.cell_volume_L = p.porosity * p.cross_section * p.dx * 1000;

validate_stability(p);

courant = p.velocity * p.dt / p.dx;
r_TCE = p.Dx_TCE * p.dt / p.dx^2;
r_cDCE = p.Dx_cDCE * p.dt / p.dx^2;
r_VC = p.Dx_VC * p.dt / p.dx^2;
w_TCE = [r_TCE + courant/2, 1 - 2*r_TCE, r_TCE - courant/2];
w_cDCE = [r_cDCE + courant/2, 1 - 2*r_cDCE, r_cDCE - courant/2];
w_VC = [r_VC + courant/2, 1 - 2*r_VC, r_VC - courant/2];

%% Grids and recording locations
x1 = 0:p.dx:p.L1;
x2 = 0:p.dx:p.L2;
t = 0:p.dt:p.t_end;
Nx1 = numel(x1);
Nx2 = numel(x2);
Nt = numel(t);

plot_times = [10, 50, 150, 200, 250, 300];
plot_indices = zeros(size(plot_times));
for k = 1:numel(plot_times)
    [~, plot_indices(k)] = min(abs(t - plot_times(k)));
end

[~, idx_z1] = min(abs(x1 - p.L1 / 2));
[~, idx_z2_2] = min(abs(x2 - 2));
[~, idx_z2_15] = min(abs(x2 - 15));
[~, idx_z2_28] = min(abs(x2 - 28));

%% Initial conditions, including one downstream ghost node per zone
C1 = zeros(Nx1 + 1, 1);
C1(1:Nx1) = p.initial_TCE;
C1(1) = 0;

C2_TCE = zeros(Nx2 + 1, 1);
C2_cDCE = zeros(Nx2 + 1, 1);
C2_VC = zeros(Nx2 + 1, 1);

source_mass = zeros(Nt, 1);
source_mass(1) = p.initial_source_mass;
source_active = true;
source_stop_time = NaN;

series = zeros(Nt, 13);
series(:, 1) = t(:);
series(1, 2:13) = [ ...
    C1(idx_z1), ...
    C2_TCE(idx_z2_2), C2_cDCE(idx_z2_2), C2_VC(idx_z2_2), ...
    C2_TCE(idx_z2_15), C2_cDCE(idx_z2_15), C2_VC(idx_z2_15), ...
    C2_TCE(idx_z2_28), C2_cDCE(idx_z2_28), C2_VC(idx_z2_28), ...
    source_mass(1), double(source_active)];

snapshots = struct('time', cell(numel(plot_times), 1), ...
                   'x', cell(numel(plot_times), 1), ...
                   'TCE', cell(numel(plot_times), 1), ...
                   'cDCE', cell(numel(plot_times), 1), ...
                   'VC', cell(numel(plot_times), 1));

%% Explicit finite-difference simulation
for step = 1:(Nt - 1)
    % Downstream zero-gradient boundary represented by a ghost node.
    C1(Nx1 + 1) = C1(Nx1 - 1);

    C1_next = C1;
    released_mass = 0;
    for node = 2:Nx1
        if source_active
            source_term = p.KT_A_V * (p.TCE_solubility - C1(node));
        else
            source_term = 0;
        end
        C1_next(node) = w_TCE(1)*C1(node - 1) + w_TCE(2)*C1(node) + ...
            w_TCE(3)*C1(node + 1) + p.dt*source_term;
        if source_active && source_term > 0
            released_mass = released_mass + source_term * p.cell_volume_L * p.dt;
        end
    end
    C1_next(1) = 0;
    C1_next(Nx1 + 1) = C1_next(Nx1 - 1);

    % Spatially accumulated source-mass bookkeeping.
    if source_active
        source_mass(step + 1) = max(0, source_mass(step) - released_mass);
        if source_mass(step + 1) <= 0.001 * p.initial_source_mass
            source_mass(step + 1) = 0;
            source_active = false;
            source_stop_time = t(step + 1);
            fprintf('Method 1 source stops releasing at %.2f days.\n', source_stop_time);
        end
    else
        source_mass(step + 1) = source_mass(step);
    end

    % Interface and downstream ghost-node conditions for reactive Zone 2.
    C2_TCE(1) = C1(Nx1);
    C2_TCE(Nx2 + 1) = C2_TCE(Nx2 - 1);
    C2_cDCE(1) = 0;
    C2_cDCE(Nx2 + 1) = C2_cDCE(Nx2 - 1);
    C2_VC(1) = 0;
    C2_VC(Nx2 + 1) = C2_VC(Nx2 - 1);

    C2_TCE_next = C2_TCE;
    C2_cDCE_next = C2_cDCE;
    C2_VC_next = C2_VC;
    C2_TCE_next(1) = C1_next(Nx1);
    C2_cDCE_next(1) = 0;
    C2_VC_next(1) = 0;

    for node = 2:Nx2
        bio_T = -p.mu_TCE * C2_TCE(node) * p.cell_density / ((p.Ks_TCE + C2_TCE(node)) * p.yield);
        bio_c = -p.mu_cDCE * C2_cDCE(node) * p.cell_density / ((p.Ks_cDCE + C2_cDCE(node)) * p.yield);
        bio_V = -p.mu_VC * C2_VC(node) * p.cell_density / ((p.Ks_VC + C2_VC(node)) * p.yield);

        C2_TCE_next(node) = w_TCE(1)*C2_TCE(node - 1) + w_TCE(2)*C2_TCE(node) + ...
            w_TCE(3)*C2_TCE(node + 1) + p.dt*bio_T;
        C2_cDCE_next(node) = w_cDCE(1)*C2_cDCE(node - 1) + w_cDCE(2)*C2_cDCE(node) + ...
            w_cDCE(3)*C2_cDCE(node + 1) + p.dt*(-bio_T + bio_c);
        C2_VC_next(node) = w_VC(1)*C2_VC(node - 1) + w_VC(2)*C2_VC(node) + ...
            w_VC(3)*C2_VC(node + 1) + p.dt*(-bio_c + bio_V);
    end

    C2_TCE_next(Nx2 + 1) = C2_TCE_next(Nx2 - 1);
    C2_cDCE_next(Nx2 + 1) = C2_cDCE_next(Nx2 - 1);
    C2_VC_next(Nx2 + 1) = C2_VC_next(Nx2 - 1);

    C1 = C1_next;
    C2_TCE = C2_TCE_next;
    C2_cDCE = C2_cDCE_next;
    C2_VC = C2_VC_next;

    series(step + 1, 2:13) = [ ...
        C1(idx_z1), ...
        C2_TCE(idx_z2_2), C2_cDCE(idx_z2_2), C2_VC(idx_z2_2), ...
        C2_TCE(idx_z2_15), C2_cDCE(idx_z2_15), C2_VC(idx_z2_15), ...
        C2_TCE(idx_z2_28), C2_cDCE(idx_z2_28), C2_VC(idx_z2_28), ...
        source_mass(step + 1), double(source_active)];

    snapshot_position = find(plot_indices == step + 1, 1);
    if ~isempty(snapshot_position)
        snapshots(snapshot_position) = combine_zones(t(step + 1), x1, x2, p.L1, ...
            C1, C2_TCE, C2_cDCE, C2_VC, Nx1, Nx2);
    end
end

%% Export outputs
series_table = array2table(series, 'VariableNames', { ...
    'time_day', 'TCE_zone1_x17_5_mg_L', ...
    'TCE_zone2_x2_mg_L', 'cDCE_zone2_x2_mg_L', 'VC_zone2_x2_mg_L', ...
    'TCE_zone2_x15_mg_L', 'cDCE_zone2_x15_mg_L', 'VC_zone2_x15_mg_L', ...
    'TCE_zone2_x28_mg_L', 'cDCE_zone2_x28_mg_L', 'VC_zone2_x28_mg_L', ...
    'source_mass_mg', 'source_active'});
csv_series = series;
concentration_columns = 2:11;
csv_concentrations = csv_series(:, concentration_columns);
csv_concentrations(abs(csv_concentrations) < 1e-12) = 0;
csv_series(:, concentration_columns) = csv_concentrations;
csv_series_table = array2table(csv_series, 'VariableNames', series_table.Properties.VariableNames);
writetable(csv_series_table, fullfile(output_dir, 'method1_timeseries.csv'));

for k = 1:numel(snapshots)
    fig = figure('Visible', 'off');
    plot(snapshots(k).x, snapshots(k).TCE, 'r-', 'LineWidth', 2); hold on;
    plot(snapshots(k).x, snapshots(k).cDCE, 'b--', 'LineWidth', 2);
    plot(snapshots(k).x, snapshots(k).VC, 'g-.', 'LineWidth', 2); hold off;
    xlabel('Distance, x (m)'); ylabel('Concentration (mg/L)');
    title(sprintf('Method 1: concentration profiles at day %.0f', snapshots(k).time));
    legend({'TCE', 'cDCE', 'VC'}, 'Location', 'best'); grid on;
    save_figure(fig, fullfile(output_dir, sprintf('spatial_day_%03.0f.png', snapshots(k).time)));
end

export_time_plot(t, series(:, 2), {'TCE'}, 'Method 1: TCE in source zone at x = 17.5 m', ...
    fullfile(output_dir, 'timeseries_zone1_x17_5m.png'));
export_time_plot(t, series(:, 3:5), {'TCE', 'cDCE', 'VC'}, 'Method 1: reactive zone at x = 2 m', ...
    fullfile(output_dir, 'timeseries_zone2_x2m.png'));
export_time_plot(t, series(:, 6:8), {'TCE', 'cDCE', 'VC'}, 'Method 1: reactive zone at x = 15 m', ...
    fullfile(output_dir, 'timeseries_zone2_x15m.png'));
export_time_plot(t, series(:, 9:11), {'TCE', 'cDCE', 'VC'}, 'Method 1: reactive zone at x = 28 m', ...
    fullfile(output_dir, 'timeseries_zone2_x28m.png'));

fig = figure('Visible', 'off');
plot(t, source_mass / 1000, 'LineWidth', 2);
xlabel('Time (day)'); ylabel('Source mass (g)');
title('Method 1: source-mass depletion'); grid on;
save_figure(fig, fullfile(output_dir, 'source_mass.png'));

results.parameters = p;
results.source_stop_time_day = source_stop_time;
results.timeseries = series_table;
results.snapshots = snapshots;
save(fullfile(output_dir, 'method1_results.mat'), 'results');
end

function snapshot = combine_zones(time_value, x1, x2, L1, C1, C2_TCE, C2_cDCE, C2_VC, Nx1, Nx2)
snapshot.time = time_value;
snapshot.x = [x1, L1 + x2(2:end)]';
snapshot.TCE = [C1(1:Nx1); C2_TCE(2:Nx2)];
snapshot.cDCE = [zeros(Nx1, 1); C2_cDCE(2:Nx2)];
snapshot.VC = [zeros(Nx1, 1); C2_VC(2:Nx2)];
end

function export_time_plot(t, concentrations, labels, plot_title, filename)
fig = figure('Visible', 'off');
plot(t, concentrations, 'LineWidth', 2);
xlabel('Time (day)'); ylabel('Concentration (mg/L)');
title(plot_title); legend(labels, 'Location', 'best'); grid on;
save_figure(fig, filename);
end

function validate_stability(p)
advective_limit = 2 * min([p.Dx_TCE, p.Dx_cDCE, p.Dx_VC]) / p.velocity;
diffusive_limit = p.dx^2 / (2 * max([p.Dx_TCE, p.Dx_cDCE, p.Dx_VC]));
if p.dx >= advective_limit
    warning('Spatial step may violate the report''s advection-dispersion criterion.');
end
if p.dt >= diffusive_limit
    warning('Time step may violate the report''s diffusion criterion.');
end
end
