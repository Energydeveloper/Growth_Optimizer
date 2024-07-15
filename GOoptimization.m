function main
    % Constants and system parameters
      % Set seed for reproducibility
    rng(123); % You can use any integer here as the seed
    T = 48;  % Time horizon
    params.rho = 1;
    params.pi = 1;
    params.epsilon = 1;
    params.gamma = 1;
    params.sigma = 1;
    params.alpha = 1;

    params.penalty_factor = 100000;  % Penalty factor for power imbalance
    params.battery_capacity = 150000;  % Battery capacity in Ah

    % Load profile and grid tariff prices
    P_d = [8047, 7943, 7749, 7498, 7240, 7031, 6864, 6780, 6776, 6885, 7227, 7631,...
           8356, 9105, 9653, 9749, 9593, 9149, 8983, 8724, 8552, 8451, 8306, 8089,...
           7776, 7602, 7504, 7531, 7581, 7765, 8082, 8556, 9061, 9663, 10319, 10797,...
           10807, 10777, 10474, 10280, 10030, 9749, 9462, 9131, 8982, 8707, 8577, 8440];
    Cg_t = [99.34, 95.51, 88.41, 87.9, 82.77, 80.29, 78.36, 75.29, 81.26, 87, 92.83,...
            91.46, 122.65, 172.54, 153.32, 128.48, 92.38, 76.01, 79.65, 72.44, 63.66,...
            60.05, 56.5, 59.83, 56.57, 62.71, 74.45, 74.71, 76.89, 76.76, 76.67, 102.19,...
            151.03, 175.66, 600, 600, 377.38, 196.43, 154.77, 162.63, 154.79, 138.49,...
            126.93, 103.93, 109.45, 108.73, 107.64, 101.32];

    xmax = 10000 * ones(1, 2 * T);
    xmin = [-10000 * ones(1, T), zeros(1, T)];
    popsize = 5000;
    dimension = 2 * T;
    MaxFEs = 4000000;

    % Execute the Growth Optimizer
    [best_solution, best_fitness, fitness_history, metrics_history] = GO(popsize, dimension, xmax, xmin, MaxFEs, P_d, Cg_t, params);

    % Calculate the battery state of charge (SoC) in percentage
    P_batt = best_solution(1:2:end);
    battery_state = params.battery_capacity;
    soc_history = zeros(1, length(P_batt));
    for i = 1:length(P_batt)
        battery_state = battery_state - P_batt(i);
        soc_history(i) = (battery_state / params.battery_capacity) * 100;
    end

    % Display results
    disp(['Best fitness: ', num2str(best_fitness)]);
    figure;
    plot(fitness_history);
    xlabel('Function Evaluations', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('Fitness', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    title('Fitness Evolution', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman',  'FontWeight', 'bold');
    
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
    plot(timeAxis, best_solution(2:2:end), 'r-', 'LineWidth', 2); % Use timeAxis here
    title('Optimal Grid Power Configuration', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    xlabel('Time (hh:mm)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('Power (kW)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    datetick('x', 'HH:MM', 'keepticks'); % Adjust the x-axis to show time
    grid on;
    
    % Plot Power Balance
    subplot(4,1,3);
    plot(timeAxis, P_d, 'k--', 'LineWidth', 2);
    hold on;
    plot(timeAxis, P_batt + best_solution(2:2:end), 'g-', 'LineWidth', 2);
    title('Power Balance', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    xlabel('Time (hh:mm)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('Power (kW)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    legend('Load Profile', 'Total Power', 'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 20);
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    datetick('x', 'HH:MM', 'keepticks'); % Adjust the x-axis to show time
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

    % Define larger figure size
    figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
    
    % Define the plot configuration function
    function configPlot(titleText, xLabel, yLabel)
        title(titleText, 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel(xLabel, 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(yLabel, 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
        grid on;
    end
    
    % Cyber Layer
    subplot(8, 1, 1); % Adjusted for 8 plots now
    plot(metrics_history.reliability_cost, 'r-', 'LineWidth', 2);
    configPlot('Reliability Cost Evolution', 'Function Evaluations', 'Cost');
    
    subplot(8, 1, 2);
    plot(metrics_history.resiliency_cost, 'g-', 'LineWidth', 2);
    configPlot('Resiliency Cost Evolution', 'Function Evaluations', 'Cost');
    
    % Physical Layer
    subplot(8, 1, 3);
    plot(metrics_history.peak_to_average, 'b-', 'LineWidth', 2);
    configPlot('Peak to Average Evolution', 'Function Evaluations', 'Ratio');
    
    subplot(8, 1, 4);
    plot(metrics_history.Battery_degradation, 'm-', 'LineWidth', 2);
    configPlot('Battery Degradation Evolution', 'Function Evaluations', 'Degradation');
    
    % Social Layer
    subplot(8, 1, 5);
    plot(metrics_history.Battery_automation, 'c-', 'LineWidth', 2);
    configPlot('Battery Automation Evolution', 'Function Evaluations', 'Automation');
    
    subplot(8, 1, 6);
    plot(metrics_history.energy_cost, 'k-', 'LineWidth', 2);
    configPlot('Energy Cost Evolution', 'Function Evaluations', 'Cost');
    
    % Penalty Analysis
    subplot(8, 1, 7);
    plot(metrics_history.power_imbalance_penalty, 'y-', 'LineWidth', 2);
    configPlot('Power Imbalance Penalty Evolution', 'Function Evaluations', 'Penalty');
    
    subplot(8, 1, 8);
    plot(metrics_history.battery_capacity_penalty, 'or-', 'LineWidth', 2);
    configPlot('Battery Capacity Penalty Evolution', 'Function Evaluations', 'Penalty');
end

function [cost, metrics] = combined_objective(x, P_d, Cg_t, params)
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
            battery_capacity_penalties(i) = 100000 * abs(battery_state -  params.battery_capacity);
        elseif battery_state > params.battery_capacity
            % Apply penalty for overcharge (exceeding the battery capacity)
            battery_capacity_penalties(i) = 100000 * (battery_state - params.battery_capacity);
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
    Battery_degradation = params.gamma * (P_batt / 5000);

    % Social layer
    Battery_automation = params.sigma * (P_batt ./ (P_grid + 1000));
    energy_cost = params.alpha .* (P_grid .* Cg_t + abs(P_batt) .* 50);

    % Sum all costs
    cost = sum(energy_cost + resiliency_cost + peak_to_average + Battery_degradation + reliability_cost + Battery_automation) + ...
           power_imbalance_penalty + battery_capacity_penalty;
    
    % Store individual metrics
    metrics.energy_cost = sum(energy_cost);
    metrics.peak_to_average = peak_to_average;
    metrics.reliability_cost = sum(reliability_cost);
    metrics.resiliency_cost = sum(resiliency_cost);
    metrics.Battery_degradation = sum(Battery_degradation);
    metrics.Battery_automation = sum(Battery_automation);
    metrics.power_imbalance_penalty = power_imbalance_penalty;
    metrics.battery_capacity_penalty = battery_capacity_penalty;
        % Display the final values of each metric
    disp(['Energy Cost: ', num2str(metrics.energy_cost)]);
    disp(['Peak to Average Ratio: ', num2str(metrics.peak_to_average)]);
    disp(['Reliability Cost: ', num2str(metrics.reliability_cost)]);
    disp(['Resiliency Cost: ', num2str(metrics.resiliency_cost)]);
    disp(['Battery Degradation: ', num2str(metrics.Battery_degradation)]);
    disp(['Battery Automation: ', num2str(metrics.Battery_automation)]);
    disp(['Power Imbalance Penalty: ', num2str(metrics.power_imbalance_penalty)]);
    disp(['Battery Capacity Penalty: ', num2str(metrics.battery_capacity_penalty)]); 
end

function y = safe_exp(x)
    y = exp(min(x, 20));  % Prevent overflow
end

% Growth Optimizer: A powerful metaheuristic algorithm for solving different optimization problems
function [gbestX,gbestfitness,gbesthistory, metrics_history] = GO(popsize,dimension,xmax,xmin,MaxFEs, P_d, Cg_t, params)
    FEs=0;
    fitness = zeros(1, popsize);
    x = xmin + (xmax - xmin) .* rand(popsize, dimension);
    gbestfitness = inf;
    gbesthistory = [];
    metrics_history = struct();

    for i = 1:popsize
        [fitness(i), metrics] = combined_objective(x(i, :), P_d, Cg_t, params);
        FEs = FEs + 1;
        if gbestfitness > fitness(i)
            gbestfitness = fitness(i);
            gbestX = x(i, :);
        end
        gbesthistory(FEs) = gbestfitness;
        metrics_history = update_metrics_history(metrics_history, metrics, FEs);
        fprintf("FEs: %d, fitness error: %e\n", FEs, gbestfitness);
    end

    while FEs < MaxFEs
        [~, ind] = sort(fitness);
        Best_X = x(ind(1), :);

        %% Learning phase
        for i = 1:popsize
            Worst_X = x(ind(randi([popsize-4,popsize])), :);
            Better_X = x(ind(randi([2,5])), :);
            random = selectID(popsize, i, 2);
            L1 = random(1);
            L2 = random(2);
            D_value1 = Best_X - Better_X;
            D_value2 = Best_X - Worst_X;
            D_value3 = Better_X - Worst_X;
            D_value4 = x(L1, :) - x(L2, :);
            Distance1 = norm(D_value1);
            Distance2 = norm(D_value2);
            Distance3 = norm(D_value3);
            Distance4 = norm(D_value4);
            rate = Distance1 + Distance2 + Distance3 + Distance4;
            LF1 = Distance1 / rate;
            LF2 = Distance2 / rate;
            LF3 = Distance3 / rate;
            LF4 = Distance4 / rate;
            SF = fitness(i) / max(fitness);
            Gap1 = LF1 * SF * D_value1;
            Gap2 = LF2 * SF * D_value2;
            Gap3 = LF3 * SF * D_value3;
            Gap4 = LF4 * SF * D_value4;
            newx = x(i, :) + Gap1 + Gap2 + Gap3 + Gap4;
            newx = max(newx, xmin);
            newx = min(newx, xmax);
            [newfitness, metrics] = combined_objective(newx, P_d, Cg_t, params);
            FEs = FEs + 1;

            % Update
            if fitness(i) > newfitness
                fitness(i) = newfitness;
                x(i, :) = newx;
            else
                if rand < 0.001 && ind(i) ~= 1
                    fitness(i) = newfitness;
                    x(i, :) = newx;
                end
            end
            if gbestfitness > fitness(i)
                gbestfitness = fitness(i);
                gbestX = x(i, :);
            end
            gbesthistory(FEs) = gbestfitness;
            metrics_history = update_metrics_history(metrics_history, metrics, FEs);
            fprintf("FEs: %d, fitness error: %e\n", FEs, gbestfitness);
        end
        if FEs >= MaxFEs
            break;
        end

        %% Reflection phase
        for i = 1:popsize
            newx = x(i, :);
            j = 1;
            while j <= dimension
                if rand < 0.3
                    R = x(ind(randi(5)), :);
                    newx(j) = x(i, j) + (R(j) - x(i, j)) * rand;
                    if rand < (0.01 + (0.1 - 0.01) * (1 - FEs / MaxFEs))
                        newx(j) = xmin(j) + (xmax(j) - xmin(j)) * rand;
                    end
                end
                j = j + 1;
            end
            newx = max(newx, xmin);
            newx = min(newx, xmax);
            [newfitness, metrics] = combined_objective(newx, P_d, Cg_t, params);
            FEs = FEs + 1;

            % Update
            if fitness(i) > newfitness
                fitness(i) = newfitness;
                x(i, :) = newx;
            else
                if rand < 0.001 && ind(i) ~= 1
                    fitness(i) = newfitness;
                    x(i, :) = newx;
                end
            end
            if gbestfitness > fitness(i)
                gbestfitness = fitness(i);
                gbestX = x(i, :);
            end
            gbesthistory(FEs) = gbestfitness;
            metrics_history = update_metrics_history(metrics_history, metrics, FEs);
            fprintf("FEs: %d, fitness error: %e\n", FEs, gbestfitness);
        end
        if FEs >= MaxFEs
            break;
        end
    end

    % Deal with the situation of too little or too much evaluation
    if FEs < MaxFEs
        gbesthistory(FEs + 1:MaxFEs) = gbestfitness;
    else
        if FEs > MaxFEs
            gbesthistory(MaxFEs + 1:end) = [];
        end
    end
end

function random = selectID(popsize, i, num)
    random = randperm(popsize, num);
    while any(random == i)
        random = randperm(popsize, num);
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
