%==========================================================================
%                  FULL-FLEDGED RCS ANALYSIS TOOL (v3.1)
%==========================================================================
% A robust, high-performance, and feature-rich tool for Radar Cross
% Section (RCS) analysis of complex 3D models.
%
% --- KEY FEATURES (v3.1) ---
% 1.  ROBUST PHYSICS ENGINE (NaN FIX):
%     - The Physical Optics (PO) function now includes a special-case
%       handler for normal incidence, fixing the division-by-zero error
%       that caused 'NaN' results at specular angles.
%
% 2.  FAILSAFE DATA SANITIZATION:
%     - A new post-processing step ensures any residual NaN/Inf values
%       are replaced with a low number, guaranteeing that all final metrics
%       (like frontal_rcs) are always valid numbers.
%
% 3.  USER-SELECTABLE POLARIZATION & OPTIONAL MEC.
% 4.  COMPLETE USER INTERFACE & REPORTING SUITE.
%==========================================================================

clear; clc; close all;

%% =================== 1. USER CONFIGURATION PANEL =====================
% --- Input File and Geometry ---
stl_file = 'Your Stl File.stl';  % The name of your 3D model file
scale_factor = 0.001;          % Use 0.001 if model is in mm, 1.0 if in m

% --- Core Simulation Parameters ---
polarization = 'VV';         % Set polarization: 'VV' or 'HH'
enable_MEC = false;          % Set to true to include edge effects, false for PO only
mec_scaling_factor = 1.0;    % Empirical factor to adjust MEC contribution (1.0 is default)

% --- Material Properties (RAM) ---% --- Material Properties (RAM) ---
% Set eta_ram = 0 for a Perfect Electric Conductor (PEC) to test shaping.
% Use a complex value for Radar Absorbent Material (e.g., 0.43 + 0.05i).
eta_ram = 0.43 + 0.13i; % Set to 0 for a Perfect Electric Conductor (PEC)

% --- Define Major Military Radar Bands for Analysis ---
bands(1).name = 'L-Band';   bands(1).freq = 1.5e9;
bands(2).name = 'S-Band';   bands(2).freq = 3.0e9;
bands(3).name = 'C-Band';   bands(3).freq = 6.0e9;
bands(4).name = 'X-Band';   bands(4).freq = 9.0e9;
bands(5).name = 'Ku-Band';  bands(5).freq = 15.0e9;
bands(6).name = 'Ka-Band';  bands(6).freq = 35.0e9;

% --- Analysis Sweep ---
azimuth_angles_deg = 0:1:359;
elevation_angle_deg = 0;

%% =================== 2. GEOMETRY PRE-PROCESSING =====================
fprintf('--- Initializing Full-Fledged RCS Analysis v3.1 ---\n');
fprintf('1. Loading and pre-processing geometry...\n');
tic_setup = tic;

TR = stlread(stl_file);
faces = TR.ConnectivityList;
vertices = TR.Points * scale_factor;

p1_f=vertices(faces(:,1),:); p2_f=vertices(faces(:,2),:); p3_f=vertices(faces(:,3),:);
face_normals = cross(p2_f - p1_f, p3_f - p1_f, 2);
face_areas = 0.5 * vecnorm(face_normals, 2, 2);
valid_faces = face_areas > 1e-9;
faces = faces(valid_faces,:);
face_normals = face_normals(valid_faces,:);
face_areas = face_areas(valid_faces,:);
face_normals = face_normals ./ (2 * face_areas);

edges = [faces(:,[1,2]); faces(:,[2,3]); faces(:,[3,1])];
edges = sort(edges, 2);
[unique_edges, ~, ~] = unique(edges, 'rows');
fprintf('   ...Setup complete in %.4f seconds.\n', toc(tic_setup));

%% =================== 3. MULTI-BAND RCS SWEEP ========================
fprintf('2. Performing Multi-Band RCS Sweep...\n');
tic_sweep_total = tic;
rcs_matrix = zeros(length(bands), length(azimuth_angles_deg));

azimuth_angles_rad = deg2rad(azimuth_angles_deg);
el = deg2rad(elevation_angle_deg);
k_inc_dirs = [cos(el)*cos(azimuth_angles_rad); cos(el)*sin(azimuth_angles_rad); sin(el)*ones(size(azimuth_angles_rad))];
k_scat_dirs = -k_inc_dirs;

switch upper(polarization)
    case 'VV'
        e_inc_pols = repmat([0; 0; 1], 1, length(azimuth_angles_deg));
    case 'HH'
        e_inc_pols = [-sin(azimuth_angles_rad); cos(azimuth_angles_rad); zeros(size(azimuth_angles_rad))];
    otherwise
        error("Invalid polarization specified. Use 'VV' or 'HH'.");
end

h_sim_waitbar = waitbar(0, 'Starting simulation...', 'Name', 'RCS Sweep Progress');

for i_band = 1:length(bands)
    waitbar_msg = sprintf('Simulating Band %d/%d: %s', i_band, length(bands), bands(i_band).name);
    waitbar(i_band / length(bands), h_sim_waitbar, waitbar_msg);
    fprintf('\n--- Simulating Band %d/%d: %s (%.1f GHz) ---\n', i_band, length(bands), bands(i_band).name, bands(i_band).freq/1e9);
    tic_band = tic;
    
    lambda = physconst('LightSpeed') / bands(i_band).freq;
    k = 2 * pi / lambda;
    
    E_scat_po = calculate_po_sweep_vec(k, k_inc_dirs, k_scat_dirs, e_inc_pols, faces, vertices, face_normals, face_areas, eta_ram);
    
    if enable_MEC
        E_scat_mec = calculate_mec_sweep_vec(k, k_inc_dirs, k_scat_dirs, vertices, unique_edges);
        E_total = E_scat_po + (mec_scaling_factor * E_scat_mec);
    else
        E_total = E_scat_po;
    end

    rcs_linear = 4 * pi * sum(abs(E_total).^2, 1);
    rcs_db = 10 * log10(rcs_linear + 1e-12);
    
    % --- FAILSAFE DATA SANITIZATION (GUARANTEES NO NaN) ---
    % Replace any NaN or Inf values with a very low number to prevent
    % errors in averaging calculations.
    rcs_db(isnan(rcs_db) | isinf(rcs_db)) = -120;
    
    bands(i_band).rcs_db = rcs_db;
    rcs_matrix(i_band, :) = rcs_db;
    
    bands(i_band).max_rcs = max(rcs_db);
    bands(i_band).avg_rcs = 10*log10(mean(10.^(rcs_db/10)));
    frontal_mask = (azimuth_angles_deg <= 45) | (azimuth_angles_deg >= 315);
    rear_mask = (azimuth_angles_deg >= 135) & (azimuth_angles_deg <= 225);
    side_mask = (azimuth_angles_deg > 45 & azimuth_angles_deg < 135) | (azimuth_angles_deg > 225 & azimuth_angles_deg < 315);
    bands(i_band).frontal_rcs = 10*log10(mean(10.^(rcs_db(frontal_mask)/10)));
    bands(i_band).rear_rcs = 10*log10(mean(10.^(rcs_db(rear_mask)/10)));
    bands(i_band).side_rcs = 10*log10(mean(10.^(rcs_db(side_mask)/10)));
    
    fprintf('   > Band simulation time: %.2f seconds.\n', toc(tic_band));
    fprintf('   > VALUE TABLE FOR CURRENT BAND:\n');
    fprintf('     - Average RCS: %8.2f dBsm, Peak RCS: %8.2f dBsm, Frontal RCS: %8.2f dBsm\n', bands(i_band).avg_rcs, bands(i_band).max_rcs, bands(i_band).frontal_rcs);
    
    elapsed_time = toc(tic_sweep_total);
    remaining_time = (elapsed_time / i_band) * (length(bands) - i_band);
    fprintf('   > Progress: %d%% complete. Est. time remaining: %.1f minutes.\n', round(i_band/length(bands)*100), remaining_time/60);
end

close(h_sim_waitbar);
fprintf('\n   ...Total sweep complete in %.2f minutes.\n', toc(tic_sweep_total)/60);

%% =================== 4. GENERATE TEXT REPORT & VISUALIZATIONS ===================
fprintf('4. Generating final analysis report and visualizations...\n');
tic_report_viz = tic;

fprintf('\n========================================================================================\n');
fprintf('                             RCS PERFORMANCE ANALYSIS REPORT\n');
fprintf('========================================================================================\n');
fprintf('Model File: %s,      Polarization: %s,      MEC Enabled: %s\n\n', stl_file, upper(polarization), string(enable_MEC));
for i = 1:length(bands)
    fprintf('--- Analysis for %s (%.1f GHz) ---\n', bands(i).name, bands(i).freq/1e9);
    fprintf('   > Average RCS        : %8.2f dBsm\n', bands(i).avg_rcs);
    fprintf('   > Peak RCS (Hotspot) : %8.2f dBsm\n', bands(i).max_rcs);
    fprintf('   > Frontal Sector RCS : %8.2f dBsm\n', bands(i).frontal_rcs);
    fprintf('   > Side Sector RCS    : %8.2f dBsm\n', bands(i).side_rcs);
    fprintf('   > Rear Sector RCS    : %8.2f dBsm\n\n', bands(i).rear_rcs);
end
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('                               CROSS-BAND PERFORMANCE SUMMARY\n');
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('%-10s | %-12s | %-12s | %-12s | %-12s | %-12s\n', 'Band', 'Avg RCS', 'Peak RCS', 'Frontal RCS', 'Side RCS', 'Rear RCS');
fprintf('%-10s | %-12s | %-12s | %-12s | %-12s | %-12s\n', '----------', '------------', '------------', '-------------', '------------', '-----------');
for i = 1:length(bands)
    fprintf('%-10s | %11.2f | %11.2f | %12.2f | %11.2f | %12.2f\n', ...
        bands(i).name, bands(i).avg_rcs, bands(i).max_rcs, bands(i).frontal_rcs, bands(i).side_rcs, bands(i).rear_rcs);
end
fprintf('========================================================================================\n\n');

num_viz_tasks = 2;
h_viz_waitbar = waitbar(0, 'Initializing...', 'Name', 'Generating Visualizations');
waitbar(0/num_viz_tasks, h_viz_waitbar, 'Generating Simplified Heatmap...');
create_simplified_heatmap(azimuth_angles_deg, bands, rcs_matrix);
waitbar(1/num_viz_tasks, h_viz_waitbar, 'Generating Polar & Cartesian Plots...');
create_2d_plots(deg2rad(azimuth_angles_deg), bands);
waitbar(1, h_viz_waitbar, 'Complete.');
close(h_viz_waitbar);
fprintf('   ...Report and visualizations generated in %.4f seconds.\n', toc(tic_report_viz));
fprintf('Analysis complete.\n');


%% =================== 5. CORE PHYSICS & PLOTTING FUNCTIONS ===================

function create_simplified_heatmap(angles_deg, bands, rcs_matrix)
    figure('Name', 'Simplified RCS Heatmap', 'Position', [50, 50, 800, 450]);
    imagesc(angles_deg, 1:length(bands), rcs_matrix); colormap(jet); h = colorbar;
    ylabel(h, 'RCS (dBsm)'); xlabel('Azimuth Angle (degrees)'); ylabel('Radar Band');
    title('RCS Heatmap: Azimuth vs. Frequency Band');
    band_labels = arrayfun(@(b) sprintf('%s (%.1f GHz)', b.name, b.freq/1e9), bands, 'UniformOutput', false);
    yticks(1:length(bands)); yticklabels(band_labels);
    xlim([0 360]); xticks(0:45:360); grid on; ax = gca; ax.Layer = 'top';
end

function create_2d_plots(angles_rad, bands)
    figure('Name', 'Multi-Band RCS Polar Signature', 'Position', [100, 100, 600, 500]);
    ax_polar = polaraxes; hold(ax_polar, 'on'); colors = jet(length(bands));
    for i = 1:length(bands); polarplot(ax_polar, angles_rad, bands(i).rcs_db, 'LineWidth', 1.5, 'Color', colors(i,:), 'DisplayName', sprintf('%.1f GHz', bands(i).freq/1e9)); end
    ax_polar.ThetaZeroLocation = 'top'; ax_polar.ThetaDir = 'clockwise';
    title('360° RCS Signature (dBsm)'); legend('show', 'Location', 'southoutside', 'NumColumns', 3);
    
    figure('Name', 'Multi-Band RCS vs. Azimuth', 'Position', [750, 100, 700, 500]);
    ax_cartesian = axes; hold(ax_cartesian, 'on');
    for i = 1:length(bands); plot(ax_cartesian, rad2deg(angles_rad), bands(i).rcs_db, 'LineWidth', 1.5, 'Color', colors(i,:), 'DisplayName', sprintf('%.1f GHz', bands(i).freq/1e9)); end
    grid(ax_cartesian, 'on'); xlim(ax_cartesian, [0 360]);
    title('RCS vs. Azimuth Angle'); xlabel(ax_cartesian, 'Azimuth Angle (degrees)'); ylabel(ax_cartesian, 'RCS (dBsm)');
    legend('show', 'Location', 'southoutside', 'NumColumns', 3);
end

% --- ROBUST PHYSICAL OPTICS FUNCTION (v3.1) ---
function E_scat = calculate_po_sweep_vec(k, k_inc_dirs, k_scat_dirs, e_inc_pols, faces, vertices, normals, face_areas, eta)
    num_angles = size(k_inc_dirs, 2);
    E_scat_total_per_angle = complex(zeros(3, num_angles));

    p1_all = vertices(faces(:,1),:); p2_all = vertices(faces(:,2),:); p3_all = vertices(faces(:,3),:);
    face_centroids = (p1_all + p2_all + p3_all) / 3;

    for i = 1:size(faces, 1)
        dot_prod_illumination = normals(i, :) * k_inc_dirs;
        is_illuminated_mask = dot_prod_illumination < 0;
        if ~any(is_illuminated_mask); continue; end
        
        k_inc_subset = k_inc_dirs(:, is_illuminated_mask);
        k_scat_subset = k_scat_dirs(:, is_illuminated_mask);
        e_inc_subset = e_inc_pols(:, is_illuminated_mask);
        num_illuminated_angles = size(k_inc_subset, 2);
        
        n_rep = repmat(normals(i,:)', 1, num_illuminated_angles);
        
        % --- ROBUSTNESS FIX for Normal Incidence (prevents NaN) ---
        e_perp = cross(k_inc_subset, n_rep);
        norm_e_perp = vecnorm(e_perp, 2, 1);
        is_normal_incidence_mask = norm_e_perp < 1e-9;
        is_oblique_incidence_mask = ~is_normal_incidence_mask;
        
        E_refl = complex(zeros(3, num_illuminated_angles));
        
        % Case 1: Oblique Incidence (stable)
        if any(is_oblique_incidence_mask)
            e_perp_oblique = e_perp(:, is_oblique_incidence_mask);
            e_perp_oblique = e_perp_oblique ./ norm_e_perp(is_oblique_incidence_mask);
            e_para_oblique = cross(n_rep(:, is_oblique_incidence_mask), e_perp_oblique);

            E_para_comp = dot(e_inc_subset(:, is_oblique_incidence_mask), e_para_oblique, 1);
            E_perp_comp = dot(e_inc_subset(:, is_oblique_incidence_mask), e_perp_oblique, 1);
            
            cos_theta_i = -dot_prod_illumination(is_illuminated_mask);
            cos_theta_oblique = cos_theta_i(is_oblique_incidence_mask);
            
            if eta == 0; Rv = -1; Rh = 1;
            else; Rv = (eta * cos_theta_oblique - 1) ./ (eta * cos_theta_oblique + 1); Rh = (cos_theta_oblique - eta) ./ (cos_theta_oblique + eta); end
            
            E_refl(:, is_oblique_incidence_mask) = E_para_comp .* Rv .* e_para_oblique + E_perp_comp .* Rh .* e_perp_oblique;
        end
        
        % Case 2: Normal Incidence (special handling)
        if any(is_normal_incidence_mask)
            if eta == 0; R_normal = -1; % PEC reflects with phase inversion
            else; R_normal = (eta - 1) / (eta + 1); end
            E_refl(:, is_normal_incidence_mask) = R_normal .* e_inc_subset(:, is_normal_incidence_mask);
        end
        
        % --- Continue with PO Integral Calculation ---
        v_vec = k_scat_subset - k_inc_subset;
        r_c = face_centroids(i, :)';
        phase_term = exp(1i * k * sum(v_vec .* r_c, 1));
        
        v_dot_n = sum(v_vec .* n_rep, 1);
        v_proj = v_vec - v_dot_n .* n_rep;
        norm_v_proj = vecnorm(v_proj, 2, 1);
        
        sinc_arg = k * norm_v_proj .* sqrt(face_areas(i)) / pi;
        sinc_term = sinc(sinc_arg / 2); % Gordon's method uses L/lambda = kL/2pi
        sinc_term(norm_v_proj < 1e-9) = 1;
        
        I = face_areas(i) .* sinc_term .* phase_term;
        pre_factor = (1j*k / (2*pi));
        E_contrib = pre_factor .* I .* cross(k_scat_subset, cross(E_refl, k_scat_subset));
        
        E_scat_total_per_angle(:, is_illuminated_mask) = E_scat_total_per_angle(:, is_illuminated_mask) + E_contrib;
    end
    E_scat = E_scat_total_per_angle;
end

% --- METHOD OF EDGE CURRENTS FUNCTION ---
function E_scat = calculate_mec_sweep_vec(k, k_inc_dirs, k_scat_dirs, vertices, edges)
    edge_vectors = vertices(edges(:,2),:) - vertices(edges(:,1),:);
    edge_len = vecnorm(edge_vectors, 2, 2);
    t_hat = edge_vectors ./ edge_len;
    edge_midpoints = (vertices(edges(:,1),:) + vertices(edges(:,2),:)) / 2;
    D_s = -1 / (2 * sqrt(2 * pi * k * 1i)); D_h = D_s;
    phase_term = exp(-1i * k * (k_scat_dirs - k_inc_dirs)' * edge_midpoints');
    sinc_arg = (k / (2*pi)) * (t_hat * (k_scat_dirs - k_inc_dirs));
    integral_term = edge_len .* sinc(sinc_arg);
    combined_term = (D_s + D_h) .* integral_term .* phase_term';
    E_scat = t_hat' * combined_term;
end

if ~exist('sinc', 'file'); sinc = @(x) sin(pi*x)./(pi*x); sinc(x(x==0))=1; end