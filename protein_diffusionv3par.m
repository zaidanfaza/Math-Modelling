
function protein_diffusion()
    % Parameters
    n = 83;
    alpha = 1;
    dt = 0.001;
    steps = 30000;
    C0 = zeros(n, 1);
    C0(83,1) = 0.5;

    % Load master matrices from CSV files (if not already loaded)
    if ~exist('masterw1normal', 'var')
        masterw1normal = csvread('master-w1 (normal).csv');
    end
    if ~exist('masterw2normal', 'var')
        masterw2normal = csvread('master-w2 (normal).csv');
    end
    if ~exist('masternnormal', 'var')
        masternnormal = csvread('master-n (normal).csv');
    end

    % Define parts mapping and regions
    parts_mapping = {83, ... %brainstem
                     [1:10, 42:51], ... %frontal
                     [11, 16:20, 52, 57:61], ... %parietal
                     [12:15, 26:27, 53:56, 67:68], ... %limbic
                     [21:24, 62:65], ... %occipital
                     [25, 28:34, 40, 66, 69:75, 81], ... %temporal
                     [35:39, 76:80, 41, 82]}; %basal ganglia
    regions = {'brainstem', 'frontal', 'parietal', 'limbic', 'occipital', 'temporal', 'basal ganglia'};

    conv_indices_all = struct();

    % Loop over different rho values
    for rho = 1.19:0.01:1.2
        fprintf('Results for rho = %.2f:\n', rho);
        
        % Solve differential equations for each rho
        concentrations_over_time_master_w1 = solve_differential_equation(masterw1normal, alpha, rho, C0, dt, steps);
        concentrations_over_time_master_w2 = solve_differential_equation(masterw2normal, alpha, rho, C0, dt, steps);
        concentrations_over_time_master_n = solve_differential_equation(masternnormal, alpha, rho, C0, dt, steps);
        
        % Rank parts by concentration time for each master matrix
        [part_indices_n, part_indices_w1, part_indices_w2] = rank_parts_by_concentration_time(concentrations_over_time_master_n, concentrations_over_time_master_w1, concentrations_over_time_master_w2, dt, parts_mapping, regions);

        % Print part indices for debugging
        fprintf('part_indices_n: %s\n', mat2str(part_indices_n));
        fprintf('part_indices_w1: %s\n', mat2str(part_indices_w1));
        fprintf('part_indices_w2: %s\n', mat2str(part_indices_w2));

        % Define literature data
        literature_data = [1, 2, 3, 4, 4, 5, 5];

        % Define conversion rules for simulation data
        conversion_rules = [1, 4, 5, 4, 5, 3, 2];

        % Convert simulation data
        conv_indices_n = convert_simulation_data(part_indices_n, conversion_rules);
        conv_indices_w1 = convert_simulation_data(part_indices_w1, conversion_rules);
        conv_indices_w2 = convert_simulation_data(part_indices_w2, conversion_rules);

        % Save converted indices in struct
        field_name_n = sprintf('rho_%.2f_master_n', rho);
        field_name_w1 = sprintf('rho_%.2f_master_w1', rho);
        field_name_w2 = sprintf('rho_%.2f_master_w2', rho);
        field_name_n = strrep(field_name_n, '.', '_');
        field_name_w1 = strrep(field_name_w1, '.', '_');
        field_name_w2 = strrep(field_name_w2, '.', '_');
        
        conv_indices_all.(field_name_n) = conv_indices_n;
        conv_indices_all.(field_name_w1) = conv_indices_w1;
        conv_indices_all.(field_name_w2) = conv_indices_w2;

        % Print converted indices for debugging
        fprintf('conv_indices_n: %s\n', mat2str(conv_indices_n));
        fprintf('conv_indices_w1: %s\n', mat2str(conv_indices_w1));
        fprintf('conv_indices_w2: %s\n', mat2str(conv_indices_w2));
    end

    % Save the struct to a file
    save('conv_indices_all.mat', 'conv_indices_all');
end

function dCdt_result = dCdt(C, master, alpha, rho, n)
    dCdt_result = zeros(size(C));
    for i = 1:n
        sum1 = alpha * C(i) * (1 - C(i));
        sum2 = 0;
        for j = 1:n
            sum2 = sum2 + master(i,j) * C(i) - master(i,j) * C(j);
        end
        dCdt_result(i) = sum1 - rho * sum2;
    end
end

function concentrations_over_time = solve_differential_equation(master, alpha, rho, C0, dt, steps)
    C = C0;
    concentrations_over_time = zeros(steps + 1, length(C0));
    concentrations_over_time(1, :) = C';
    for step = 1:steps
        dC = dCdt(C, master, alpha, rho, length(C0));
        C = C + dt * dC;
        concentrations_over_time(step + 1, :) = C';
    end
end

function [part_indices_n, part_indices_w1, part_indices_w2] = rank_parts_by_concentration_time(concentrations_n, concentrations_w1, concentrations_w2, dt, parts_mapping, regions)
    % Initialize arrays to store the time when each part reaches 0.5 concentration
    time_to_0_5_n = zeros(length(parts_mapping), 1);
    time_to_0_5_w1 = zeros(length(parts_mapping), 1);
    time_to_0_5_w2 = zeros(length(parts_mapping), 1);

    % Find the time when each part reaches 0.5 concentration for master_n
    for part_idx = 1:length(parts_mapping)
        part_nodes = parts_mapping{part_idx};
        first_node_time_n = inf;
        for node_idx = part_nodes
            time_idx = find(concentrations_n(:, node_idx) >= 0.5, 1);
            if ~isempty(time_idx)
                first_node_time_n = min(first_node_time_n, time_idx * dt);
            end
        end
        time_to_0_5_n(part_idx) = first_node_time_n;
    end

    % Find the time when each part reaches 0.5 concentration for master_w1
    for part_idx = 1:length(parts_mapping)
        part_nodes = parts_mapping{part_idx};
        first_node_time_w1 = inf;
        for node_idx = part_nodes
            time_idx = find(concentrations_w1(:, node_idx) >= 0.5, 1);
            if ~isempty(time_idx)
                first_node_time_w1 = min(first_node_time_w1, time_idx * dt);
            end
        end
        time_to_0_5_w1(part_idx) = first_node_time_w1;
    end

    % Find the time when each part reaches 0.5 concentration for master_w2
    for part_idx = 1:length(parts_mapping)
        part_nodes = parts_mapping{part_idx};
        first_node_time_w2 = inf;
        for node_idx = part_nodes
            time_idx = find(concentrations_w2(:, node_idx) >= 0.5, 1);
            if ~isempty(time_idx)
                first_node_time_w2 = min(first_node_time_w2, time_idx * dt);
            end
        end
        time_to_0_5_w2(part_idx) = first_node_time_w2;
    end

    % Sort the times for each master matrix and get the corresponding part numbers
    [sorted_times_n, part_indices_n] = sort(time_to_0_5_n);
    [sorted_times_w1, part_indices_w1] = sort(time_to_0_5_w1);
    [sorted_times_w2, part_indices_w2] = sort(time_to_0_5_w2);

    % Display the ranked parts for master_n
    fprintf('Ranked parts based on time to reach concentration of 0.5 for master_n:\n');
    for rank = 1:length(part_indices_n)
        fprintf('%s: %.3f\n', regions{part_indices_n(rank)}, sorted_times_n(rank));
    end

    % Display the ranked parts for master_w1
    fprintf('\nRanked parts based on time to reach concentration of 0.5 for master_w1:\n');
    for rank = 1:length(part_indices_w1)
        fprintf('%s: %.3f\n', regions{part_indices_w1(rank)}, sorted_times_w1(rank));
    end

    % Display the ranked parts for master_w2
    fprintf('\nRanked parts based on time to reach concentration of 0.5 for master_w2:\n');
    for rank = 1:length(part_indices_w2)
        fprintf('%s: %.3f\n', regions{part_indices_w2(rank)}, sorted_times_w2(rank));
    end
end

function converted_sim_data = convert_simulation_data(sim_data, conversion_rules)
    % Initialize the converted simulation data
    converted_sim_data = zeros(size(sim_data));
    
    % Apply the conversion rules
    for i = 1:length(sim_data)
        converted_sim_data(i) = conversion_rules(sim_data(i));
    end
end
