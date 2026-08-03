function results = run_method2_prb(output_root)
%RUN_METHOD2_PRB Simulate a permeable reactive biobarrier treatment layout.
%
% The model contains a 35 m source zone, a 20 m nonreactive transition
% zone, and a 10 m permeable reactive barrier. Biodegradation is localized
% in the barrier.

if nargin < 1 || strlength(string(output_root)) == 0
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    output_root = fullfile(repo_root, 'results', 'generated');
end
output_dir = fullfile(output_root, 'method2');
ensure_directory(output_dir);

%% Active executable parameters
p.alpha_L = 0.015;
p.tortuosity = 0.8;
p.K = 17;                     % m/day
p.porosity = 0.38;
p.gradient = 0.0033;
p.velocity = p.K / p.porosity * p.gradient;

p.D0_TCE = 10.1e-6 * 3600 * 24 * 1e-4;
p.D0_cDCE = 11.4e-6 * 3600 * 24 * 1e-4;
p.D0_VC = 13.3e-6 * 3600 * 24 * 1e-4;
p.Dx_TCE = p.velocity * p.alpha_L + p.porosity * p.tortuosity * p.D0_TCE;
p.Dx_cDCE = p.velocity * p.alpha_L + p.porosity * p.tortuosity * p.D0_cDCE;
p.Dx_VC = p.velocity * p.alpha_L + p.porosity * p.tortuosity * p.D0_VC;

p.dx = 5 / 167;               % m; uniform grid spacing that exactly represents all zone endpoints
p.dt = 0.05;
p.L1 = 35;
p.L2 = 20;
p.L3 = 10;
p.t_end = 700;                % report figures end at 620 days
p.KT_A_V = 7.8e-4;
p.TCE_solubility = 1300;
p.initial_TCE = 10;
p.initial_source_mass = 1e9;

p.mu_TCE = 3.42;
p.mu_cDCE = 0.16;
p.mu_VC = 0.52;
p.cell_density = 1e10;
p.Ks_TCE = 10 * 131.38;
p.Ks_cDCE = 3.3 * 96.94;
p.Ks_VC = 2.6 * 62.5;
p.yield = 5.2e8 / 35.45;
p.cross_section = 10.67 * 35;
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
x3 = 0:p.dx:p.L3;
t = 0:p.dt:p.t_end;
Nx1 = numel(x1); Nx2 = numel(x2); Nx3 = numel(x3); Nt = numel(t);

plot_times = [100, 200, 300, 400, 500, 620];
plot_indices = zeros(size(plot_times));
for k = 1:numel(plot_times)
    [~, plot_indices(k)] = min(abs(t - plot_times(k)));
end

[~, idx_z1] = min(abs(x1 - p.L1/2));
[~, idx_z2] = min(abs(x2 - p.L2/2));
[~, idx_z3_inlet] = min(abs(x3 - 0.3));
[~, idx_z3_mid] = min(abs(x3 - p.L3/2));

%% Initial conditions
C1 = zeros(Nx1 + 1, 1);
C1(1:Nx1) = p.initial_TCE;
C1(1) = 0;
C2_TCE = zeros(Nx2 + 1, 1);
C3_TCE = zeros(Nx3 + 1, 1);
C3_cDCE = zeros(Nx3 + 1, 1);
C3_VC = zeros(Nx3 + 1, 1);

source_mass = zeros(Nt, 1);
source_mass(1) = p.initial_source_mass;
source_active = true;
source_stop_time = NaN;

series = zeros(Nt, 12);
series(:, 1) = t(:);
series(1, 2:12) = [C1(idx_z1), C2_TCE(idx_z2), ...
    C3_TCE(idx_z3_inlet), C3_cDCE(idx_z3_inlet), C3_VC(idx_z3_inlet), ...
    C3_TCE(idx_z3_mid), C3_cDCE(idx_z3_mid), C3_VC(idx_z3_mid), ...
    source_mass(1), double(source_active), p.velocity];

snapshots = struct('time', cell(numel(plot_times), 1), ...
                   'x', cell(numel(plot_times), 1), ...
                   'TCE', cell(numel(plot_times), 1), ...
                   'cDCE', cell(numel(plot_times), 1), ...
                   'VC', cell(numel(plot_times), 1));

%% Explicit finite-difference simulation
for step = 1:(Nt - 1)
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

    if source_active
        source_mass(step + 1) = max(0, source_mass(step) - released_mass);
        if source_mass(step + 1) <= 0.001 * p.initial_source_mass
            source_mass(step + 1) = 0;
            source_active = false;
            source_stop_time = t(step + 1);
            fprintf('Method 2 source stops releasing at %.2f days.\n', source_stop_time);
        end
    else
        source_mass(step + 1) = source_mass(step);
    end

    % Nonreactive transition zone.
    C2_TCE(1) = C1(Nx1);
    C2_TCE(Nx2 + 1) = C2_TCE(Nx2 - 1);
    C2_TCE_next = C2_TCE;
    C2_TCE_next(1) = C1_next(Nx1);
    for node = 2:Nx2
        C2_TCE_next(node) = w_TCE(1)*C2_TCE(node - 1) + w_TCE(2)*C2_TCE(node) + ...
            w_TCE(3)*C2_TCE(node + 1);
    end
    C2_TCE_next(Nx2 + 1) = C2_TCE_next(Nx2 - 1);

    % Reactive barrier.
    C3_TCE(1) = C2_TCE(Nx2);
    C3_TCE(Nx3 + 1) = C3_TCE(Nx3 - 1);
    C3_cDCE(1) = 0;
    C3_cDCE(Nx3 + 1) = C3_cDCE(Nx3 - 1);
    C3_VC(1) = 0;
    C3_VC(Nx3 + 1) = C3_VC(Nx3 - 1);

    C3_TCE_next = C3_TCE;
    C3_cDCE_next = C3_cDCE;
    C3_VC_next = C3_VC;
    C3_TCE_next(1) = C2_TCE_next(Nx2);
    C3_cDCE_next(1) = 0;
    C3_VC_next(1) = 0;

    for node = 2:Nx3
        bio_T = -p.mu_TCE * C3_TCE(node) * p.cell_density / ((p.Ks_TCE + C3_TCE(node)) * p.yield);
        bio_c = -p.mu_cDCE * C3_cDCE(node) * p.cell_density / ((p.Ks_cDCE + C3_cDCE(node)) * p.yield);
        bio_V = -p.mu_VC * C3_VC(node) * p.cell_density / ((p.Ks_VC + C3_VC(node)) * p.yield);

        C3_TCE_next(node) = w_TCE(1)*C3_TCE(node - 1) + w_TCE(2)*C3_TCE(node) + ...
            w_TCE(3)*C3_TCE(node + 1) + p.dt*bio_T;
        C3_cDCE_next(node) = w_cDCE(1)*C3_cDCE(node - 1) + w_cDCE(2)*C3_cDCE(node) + ...
            w_cDCE(3)*C3_cDCE(node + 1) + p.dt*(-bio_T + bio_c);
        C3_VC_next(node) = w_VC(1)*C3_VC(node - 1) + w_VC(2)*C3_VC(node) + ...
            w_VC(3)*C3_VC(node + 1) + p.dt*(-bio_c + bio_V);
    end
    C3_TCE_next(Nx3 + 1) = C3_TCE_next(Nx3 - 1);
    C3_cDCE_next(Nx3 + 1) = C3_cDCE_next(Nx3 - 1);
    C3_VC_next(Nx3 + 1) = C3_VC_next(Nx3 - 1);

    C1 = C1_next;
    C2_TCE = C2_TCE_next;
    C3_TCE = C3_TCE_next;
    C3_cDCE = C3_cDCE_next;
    C3_VC = C3_VC_next;

    series(step + 1, 2:12) = [C1(idx_z1), C2_TCE(idx_z2), ...
        C3_TCE(idx_z3_inlet), C3_cDCE(idx_z3_inlet), C3_VC(idx_z3_inlet), ...
        C3_TCE(idx_z3_mid), C3_cDCE(idx_z3_mid), C3_VC(idx_z3_mid), ...
        source_mass(step + 1), double(source_active), p.velocity];

    snapshot_position = find(plot_indices == step + 1, 1);
    if ~isempty(snapshot_position)
        snapshots(snapshot_position) = combine_zones(t(step + 1), x1, x2, x3, p.L1, p.L2, ...
            C1, C2_TCE, C3_TCE, C3_cDCE, C3_VC, Nx1, Nx2, Nx3);
    end
end

%% Export outputs
series_table = array2table(series, 'VariableNames', { ...
    'time_day', 'TCE_zone1_mid_mg_L', 'TCE_zone2_mid_mg_L', ...
    'TCE_PRB_inlet_mg_L', 'cDCE_PRB_inlet_mg_L', 'VC_PRB_inlet_mg_L', ...
    'TCE_PRB_mid_mg_L', 'cDCE_PRB_mid_mg_L', 'VC_PRB_mid_mg_L', ...
    'source_mass_mg', 'source_active', 'pore_velocity_m_day'});
csv_series = series;
concentration_columns = 2:9;
csv_concentrations = csv_series(:, concentration_columns);
csv_concentrations(abs(csv_concentrations) < 1e-12) = 0;
csv_series(:, concentration_columns) = csv_concentrations;
csv_series_table = array2table(csv_series, 'VariableNames', series_table.Properties.VariableNames);
writetable(csv_series_table, fullfile(output_dir, 'method2_timeseries.csv'));

for k = 1:numel(snapshots)
    fig = figure('Visible', 'off');
    plot(snapshots(k).x, snapshots(k).TCE, 'r-', 'LineWidth', 2); hold on;
    plot(snapshots(k).x, snapshots(k).cDCE, 'b--', 'LineWidth', 2);
    plot(snapshots(k).x, snapshots(k).VC, 'g-.', 'LineWidth', 2); hold off;
    xlabel('Distance, x (m)'); ylabel('Concentration (mg/L)');
    title(sprintf('Method 2: concentration profiles at day %.0f', snapshots(k).time));
    legend({'TCE', 'cDCE', 'VC'}, 'Location', 'best'); grid on;
    save_figure(fig, fullfile(output_dir, sprintf('spatial_day_%03.0f.png', snapshots(k).time)));
end

export_time_plot(t, series(:, 2), {'TCE'}, 'Method 2: TCE in source-zone midpoint', ...
    fullfile(output_dir, 'timeseries_zone1_midpoint.png'));
export_time_plot(t, series(:, 3), {'TCE'}, 'Method 2: TCE in transition-zone midpoint', ...
    fullfile(output_dir, 'timeseries_zone2_midpoint.png'));
export_time_plot(t, series(:, 4:6), {'TCE', 'cDCE', 'VC'}, 'Method 2: species near PRB inlet', ...
    fullfile(output_dir, 'timeseries_prb_near_inlet.png'));
export_time_plot(t, series(:, 7:9), {'TCE', 'cDCE', 'VC'}, 'Method 2: species at PRB midpoint', ...
    fullfile(output_dir, 'timeseries_prb_midpoint.png'));

fig = figure('Visible', 'off');
plot(t, source_mass / 1000, 'LineWidth', 2);
xlabel('Time (day)'); ylabel('Source mass (g)');
title('Method 2: source-mass depletion'); grid on;
save_figure(fig, fullfile(output_dir, 'source_mass.png'));

results.parameters = p;
results.source_stop_time_day = source_stop_time;
results.timeseries = series_table;
results.snapshots = snapshots;
save(fullfile(output_dir, 'method2_results.mat'), 'results');
end

function snapshot = combine_zones(time_value, x1, x2, x3, L1, L2, C1, C2_TCE, C3_TCE, C3_cDCE, C3_VC, Nx1, Nx2, Nx3)
snapshot.time = time_value;
snapshot.x = [x1, L1 + x2(2:end), L1 + L2 + x3(2:end)]';
snapshot.TCE = [C1(1:Nx1); C2_TCE(2:Nx2); C3_TCE(2:Nx3)];
snapshot.cDCE = [zeros(Nx1 + Nx2 - 1, 1); C3_cDCE(2:Nx3)];
snapshot.VC = [zeros(Nx1 + Nx2 - 1, 1); C3_VC(2:Nx3)];
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
