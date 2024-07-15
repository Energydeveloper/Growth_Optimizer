function run_SA
    % Constants and system parameters
    rng(123); % Set seed for reproducibility
    T = 48;  % Time horizon
    params.rho = 1;
    params.pi = 1;
    params.epsilon = 1;
    params.gamma = 1;
    params.sigma = 1;
    params.alpha = 1;
    params.penalty_factor = 100;  % Penalty factor for power imbalance
    params.battery_capacity = 10;  % Battery capacity in Ah

    % Load profile and grid tariff prices
    P_d = [8047, 7943, 7749, 7498, 7240, 7031, 6864, 6780, 6776, 6885, 7227, 7631, ...
           8356, 9105, 9653, 9749, 9593, 9149, 8983, 8724, 8552, 8451, 8306, 8089, ...
           7776, 7602, 7504, 7531, 7581, 7765, 8082, 8556, 9061, 9663, 10319, 10797, ...
           10807, 10777, 10474, 10280, 10030, 9749, 9462, 9131, 8982, 8707, 8577, 8440];
    Cg_t = [99.34, 95.51, 88.41, 87.9, 82.77, 80.29, 78.36, 75.29, 81.26, 87, 92.83, ...
            91.46, 122.65, 172.54, 153.32, 128.48, 92.38, 76.01, 79.65, 72.44, 63.66, ...
            60.05, 56.5, 59.83, 56.57, 62.71, 74.45, 74.71, 76.89, 76.76, 76.67, 102.19, ...
            151.03, 175.66, 600, 600, 377.38, 196.43, 154.77, 162.63, 154.79, 138.49, ...
            126.93, 103.93, 109.45, 108.73, 107.64, 101.32];

    % Normalize P_d and Cg_t
    P_d_max = max(P_d);
    Cg_t_max = max(Cg_t);

    P_d = P_d / P_d_max;
    Cg_t = Cg_t / Cg_t_max;

    % Print normalized values
    disp('Normalized P_d:');
    disp(P_d);
    disp('Normalized Cg_t:');
    disp(Cg_t);

    % Define problem bounds
    xmax = 0.9 * ones(1, 2 * T);
    xmin = [-0.9 * ones(1, T), zeros(1, T)];

    % Run Simulated Annealing
    [best_solution, best_fitness, fitness_history, metrics_history] = runSA(2 * T, xmin, xmax, P_d, Cg_t, params);

    % Calculate the battery state of charge (SoC) in percentage for the first solution
    P_batt = best_solution(1:2:end);
    P_grid = best_solution(2:2:end);
    battery_state = params.battery_capacity;
    soc_history = zeros(1, length(P_batt));
    for i = 1:length(P_batt)
        battery_state = battery_state - P_batt(i);
        soc_history(i) = (battery_state / params.battery_capacity) * 100;
    end

    % Print final objective values for the first solution
    final_objectives = combined_objectives(best_solution, P_d, Cg_t, params);
    disp('Final Objectives:');
    disp(final_objectives);

    % Define time axis in hh:mm format for 48 time steps each representing half an hour
    timeAxis = datetime(2021, 1, 1, 0, 0, 0) + minutes((0:30:30*(48-1)));
    
    % Create large plots
    figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]); % Make the figure large
    
    % Plot Optimal Battery Power Configuration
    subplot(4,1,1);
    plot(timeAxis, P_batt, 'b-', 'LineWidth', 2); % Use timeAxis here
    title('Optimal Battery Power Configuration', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    xlabel('Time (hh:mm)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('Power (kW)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    datetick('x', 'HH:MM', 'keepticks'); % Adjust the x-axis to show time
    grid on;
    
    % Plot Optimal Grid Power Configuration
    subplot(4,1,2);
    plot(timeAxis, P_grid, 'r-', 'LineWidth', 2); % Use timeAxis here
    title('Optimal Grid Power Configuration', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    xlabel('Time (hh:mm)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('Power (kW)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    datetick('x', 'HH:MM', 'keepticks'); % Adjust the x-axis to show time
    grid on;
    
    % Plot Power Balance
    subplot(4,1,3);
    plot(timeAxis, P_d, 'k--', 'LineWidth', 2); % Plot normalized values
    hold on;
    plot(timeAxis, P_batt + P_grid, 'g-', 'LineWidth', 2); % Plot normalized values
    title('Power Balance', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    xlabel('Time (hh:mm)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('Normalized Power', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    legend('Load Profile', 'Total Power', 'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 20);
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    datetick('x', 'HH:MM', 'keepticks'); % Adjust the x-axis to show time
    ylim([min([P_d, P_batt + P_grid]), max([P_d, P_batt + P_grid])]); % Fit the y-axis limits to the data
    grid on;

    % Plot Battery State of Charge (SoC)
    subplot(4,1,4);
    plot(timeAxis, soc_history, 'm-', 'LineWidth', 2);
    title('Battery State of Charge (SoC)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    xlabel('Time (hh:mm)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('SoC (%)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    datetick('x', 'HH:MM', 'keepticks'); % Adjust the x-axis to show time
    grid on;
end

% Simulated Annealing Optimizer
function [best_solution, best_fitness, fitness_history, metrics_history] = runSA(dimension, xmin, xmax, P_d, Cg_t, params)
    global fitness_history_global metrics_history_global;
    fitness_history_global = [];
    metrics_history_global = struct();

    options = optimoptions('simulannealbnd', ...
        'MaxIterations', 10000, ...
        'OutputFcn', @(optimValues, flag) saOutputFcn(optimValues, flag, P_d, Cg_t, params), ...
        'Display', 'iter');
    
    [best_solution, best_fitness] = simulannealbnd(@(x) combined_objective_wrapper(x, P_d, Cg_t, params), ...
        xmin + (xmax - xmin) .* rand(1, dimension), xmin, xmax, options);

    fitness_history = fitness_history_global;
    metrics_history = metrics_history_global;
end

function [cost, metrics] = combined_objective_wrapper(x, P_d, Cg_t, params)
    [cost, metrics] = combined_objective(x, P_d, Cg_t, params);
end

function stop = saOutputFcn(optimValues, flag, P_d, Cg_t, params)
    global fitness_history_global metrics_history_global;
    stop = false;
    switch flag
        case 'iter'
            [cost, metrics] = combined_objective_wrapper(optimValues.x, P_d, Cg_t, params);
            fitness_history_global = [fitness_history_global; cost];
            metrics_history_global = update_metrics_history(metrics_history_global, metrics, length(fitness_history_global));
    end
end

function metrics_history = update_metrics_history(metrics_history, metrics, FEs)
    metric_names = fieldnames(metrics);
    for i = 1:numel(metric_names)
        if ~isfield(metrics_history, metric_names{i})
            metrics_history.(metric_names{i}) = [];
        end
        metrics_history.(metric_names{i})(FEs) = metrics.(metric_names{i});
    end
end

function [objectives, constraints] = combined_objectives(x, P_d, Cg_t, params)
    T = numel(P_d);
    P_batt = x(1:2:end);
    P_grid = x(2:2:end);

    % Initialize the battery state based on its maximum capacity
    battery_state = params.battery_capacity; % Assuming this is the initial full capacity
    
    % Array to hold the updated battery state over time and penalties for capacity issues
    updated_battery_states = zeros(1, length(P_batt));
    battery_capacity_penalties = zeros(1, length(P_batt));
    
    % Calculate battery state update and enforce constraints
    for i = 1:length(P_batt)
        % Update battery state
        battery_state = battery_state - P_batt(i); % Discharge if P_batt is positive, charge if negative
        
        % Store the updated state
        updated_battery_states(i) = battery_state;
        
        % Apply penalty for deviations from the desired state
        if battery_state < 0
            % Apply penalty for under-discharge (falling below zero)
            battery_capacity_penalties(i) = 1 * abs(battery_state - params.battery_capacity);
        elseif battery_state > params.battery_capacity
            % Apply penalty for overcharge (exceeding the battery capacity)
            battery_capacity_penalties(i) = 1 * (battery_state - params.battery_capacity);
        end
    end
    
    % Sum up all penalties to get the total battery capacity penalty
    battery_capacity_penalty = sum(battery_capacity_penalties);
    
    power_imbalance_penalty = params.penalty_factor * sum(abs(P_batt + P_grid - P_d));

    % Cyber layer
    reliability_cost = params.rho * safe_exp(-params.rho * (P_batt + P_grid - P_d));
    resiliency_cost = params.pi * ((P_batt + P_grid) ./ P_d);

    % Physical layer
    peak_to_average = params.epsilon * (max(P_grid + P_batt) / mean(P_grid + P_batt));
    Battery_degradation = params.gamma * (P_batt / 1);

    % Social layer
    Battery_automation = params.sigma * (P_batt ./ (P_grid + 1));
    energy_cost = params.alpha .* (P_grid .* Cg_t + abs(P_batt) .* 0.1);

    % Sum all costs
    objectives = [sum(energy_cost), sum(resiliency_cost), peak_to_average, sum(Battery_degradation), ...
                  sum(reliability_cost), sum(Battery_automation), power_imbalance_penalty + battery_capacity_penalty];

    % Store individual metrics
    metrics.energy_cost = sum(energy_cost);
    metrics.peak_to_average = peak_to_average;
    metrics.reliability_cost = sum(reliability_cost);
    metrics.resiliency_cost = sum(resiliency_cost);
    metrics.Battery_degradation = sum(Battery_degradation);
    metrics.Battery_automation = sum(Battery_automation);
    metrics.power_imbalance_penalty = power_imbalance_penalty;
    metrics.battery_capacity_penalty = battery_capacity_penalty;

    % Define nonlinear constraints for power balance
    constraints = power_balance_constraint(x, P_d);
end

function [c, ceq] = power_balance_constraint(x, P_d)
    T = numel(P_d);
    P_batt = x(1:2:end);
    P_grid = x(2:2:end);

    % Power balance constraint
    ceq = P_batt + P_grid - P_d; % Equality constraint (should be zero)
    c = []; % No inequality constraints
end

function y = safe_exp(x)
    y = exp(min(x, 20));  % Prevent overflow
end

function [cost, metrics] = combined_objective(x, P_d, Cg_t, params)
    [objectives, ~] = combined_objectives(x, P_d, Cg_t, params);
    cost = sum(objectives); % Combine the objectives into a single cost for SA
    metrics = struct('energy_cost', objectives(1), 'resiliency_cost', objectives(2), 'peak_to_average', objectives(3), ...
                     'Battery_degradation', objectives(4), 'reliability_cost', objectives(5), 'Battery_automation', objectives(6), ...
                     'power_imbalance_penalty', 0, 'battery_capacity_penalty', objectives(7));
end
