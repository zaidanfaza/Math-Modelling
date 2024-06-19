function calculate_abs_diff()
    % Load the saved conv_indices_all struct
    load('conv_indices_all.mat', 'conv_indices_all');
    
    % Define literature data
    %literature_data = [1, 2, 3, 4, 4, 5, 5];
    literature_data = [1, 2, 3, 4, 4, 4, 4];
    
    % Initialize arrays to store total differences
    total_diffs_n = [];
    total_diffs_w1 = [];
    total_diffs_w2 = [];
    rho_values = 0.01:0.01:10;  % Adjust rho_values as needed
    
    % Loop over the saved data
    fields = fieldnames(conv_indices_all);
    for i = 1:length(fields)
        field_name = fields{i};
        conv_indices = conv_indices_all.(field_name);
        
        % Calculate the total sum of absolute differences
        diff = conv_indices - literature_data';
        abs_diff = abs(diff);
        total_diff = sum(abs_diff);
        
        % Store the results in corresponding arrays
        if contains(field_name, 'master_n')
            total_diffs_n = [total_diffs_n, total_diff];
        elseif contains(field_name, 'master_w1')
            total_diffs_w1 = [total_diffs_w1, total_diff];
        elseif contains(field_name, 'master_w2')
            total_diffs_w2 = [total_diffs_w2, total_diff];
        end
    end
    
    % Display the results
    fprintf('Total sum of absolute differences for master_n:\n');
    disp(total_diffs_n);
    fprintf('Total sum of absolute differences for master_w1:\n');
    disp(total_diffs_w1);
    fprintf('Total sum of absolute differences for master_w2:\n');
    disp(total_diffs_w2);
    
    % Define a custom discrete colormap with 6 colors
    discrete_colors = [0, 1, 0;    % Green
                       1, 1, 0;    % Yellow
                       1, 0.5, 0;  % Orange
                       1, 0, 0;    % Red
                       0.65, 0.16, 0.16; % Brown
                       0, 0, 0];   % Black

    % Plot separate heatmaps
    figure;
    
    subplot(3, 1, 1);
    imagesc(rho_values, 1, total_diffs_n);
    colormap(gca, discrete_colors); % Use the custom discrete colormap
    caxis([0 10]); % Set color limits based on the max value of 10
    xlabel('k');
    ylabel('master-n');
    title('Error score for master-n');
    set(gca, 'YTick', []); % Remove Y ticks for a simpler look
    
    subplot(3, 1, 2);
    imagesc(rho_values, 1, total_diffs_w1);
    colormap(gca, discrete_colors); % Use the custom discrete colormap
    caxis([0 10]); % Set color limits based on the max value of 10
    xlabel('k');
    ylabel('master-w1');
    title('Error score for master-w1');
    set(gca, 'YTick', []); % Remove Y ticks for a simpler look
    
    subplot(3, 1, 3);
    imagesc(rho_values, 1, total_diffs_w2);
    colormap(gca, discrete_colors); % Use the custom discrete colormap
    caxis([0 10]); % Set color limits based on the max value of 10
    xlabel('k');
    ylabel('master-w2');
    title('Error score for master-w2');
    set(gca, 'YTick', []); % Remove Y ticks for a simpler look

    % Add a color bar at the bottom
    h = colorbar('Ticks', [0, 2, 4, 6, 8, 10], 'TickLabels', {'0', '2', '4', '6', '8', '10'}, 'Position', [0.92, 0.11, 0.02, 0.815]); 
    colormap(h, discrete_colors); % Apply the custom discrete colormap
    caxis([0 10]); % Ensure color limits are consistent
    ylabel(h, 'Error Score');
end
