function run_nsga2
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
    params.battery_capacity =10;  % Battery capacity in Ah

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

    % NSGA-II parameters
    popsize = 500;
    max_gen = 200;
    options = optimoptions('gamultiobj', ...
        'PopulationSize', popsize, ...
        'MaxGenerations', max_gen, ...
        'ParetoFraction', 0.35, ...
        'PlotFcn', {@gaplotpareto}, ...
        'Display', 'iter');

    % Run NSGA-II
    [solutions, fitness] = gamultiobj(@(x) combined_objectives(x, P_d, Cg_t, params), ...
        2 * T, [], [], [], [], xmin, xmax, options);

    % Display results
    disp('Solutions:');
    disp(solutions);
    disp('Fitness:');
    disp(fitness);

    % Calculate the battery state of charge (SoC) in percentage for the first solution
    P_batt = solutions(1, 1:2:end);
    P_grid = solutions(1, 2:2:end);
    battery_state = params.battery_capacity;
    soc_history = zeros(1, length(P_batt));
    for i = 1:length(P_batt)
        battery_state = battery_state - P_batt(i);
        soc_history(i) = (battery_state / params.battery_capacity) * 100;
    end

    % Print final objective values for the first solution
    final_objectives = combined_objectives(solutions(1, :), P_d, Cg_t, params);
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

function objectives = combined_objectives(x, P_d, Cg_t, params)
    T = numel(P_d);
    P_batt = x(1:2:end);
    P_grid = x(2:2:end);

    % Calculate metrics
    battery_state = params.battery_capacity;
    battery_capacity_penalties = zeros(1, length(P_batt));
    
    for i = 1:length(P_batt)
        battery_state = battery_state - P_batt(i);
        if battery_state < 0
            battery_capacity_penalties(i) = 5 * abs(battery_state -  params.battery_capacity);
        elseif battery_state > params.battery_capacity
            battery_capacity_penalties(i) = 5 * (battery_state - params.battery_capacity);
        end
    end
    
    battery_capacity_penalty = sum(battery_capacity_penalties);
    power_imbalance_penalty = params.penalty_factor * sum(abs(P_batt + P_grid - P_d).^2);
    reliability_cost = params.rho * safe_exp(-params.rho * (P_batt + P_grid - P_d));
    resiliency_cost = params.pi * ((P_batt + P_grid) ./ P_d);
    peak_to_average = params.epsilon * (max(P_grid + P_batt) / mean(P_grid + P_batt));
    Battery_degradation = params.gamma * (P_batt / 1);
    Battery_automation = params.sigma * (P_batt ./ (P_grid + 1));
    energy_cost = params.alpha .* (P_grid .* Cg_t + abs(P_batt) .* 0.1);

    % Define objectives
    objectives = [
        sum(energy_cost);
        sum(reliability_cost);
        sum(resiliency_cost);
        peak_to_average;
        sum(Battery_degradation);
        sum(Battery_automation);
        power_imbalance_penalty;
        battery_capacity_penalty;
    ];
end

function y = safe_exp(x)
    y = exp(min(x, 20));  % Prevent overflow
end
