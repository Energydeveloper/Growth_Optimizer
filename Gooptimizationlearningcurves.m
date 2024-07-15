function main
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

    % Load profiles and grid tariff prices for different locations
    load_profiles = {
        [8047, 7943, 7749, 7498, 7240, 7031, 6864, 6780, 6776, 6885, 7227, 7631, ...
        8356, 9105, 9653, 9749, 9593, 9149, 8983, 8724, 8552, 8451, 8306, 8089, ...
        7776, 7602, 7504, 7531, 7581, 7765, 8082, 8556, 9061, 9663, 10319, 10797, ...
        10807, 10777, 10474, 10280, 10030, 9749, 9462, 9131, 8982, 8707, 8577, 8440], ...
        [5840, 5655, 5550, 5499, 5433, 5376, 5424, 5492, 5579, 5675, 5899, 6101, ...
        6513, 6818, 6835, 6629, 6254, 5825, 5450, 5156, 4981, 4959, 4893, 4871, ...
        4929, 5065, 5100, 5162, 5340, 5603, 5904, 6294, 6714, 7060, 7340, 7591, ...
        7602, 7558, 7394, 7299, 7149, 6920, 6787, 6554, 6363, 6174, 6115, 5909], ...
        [1525, 1541, 1497, 1475, 1419, 1385, 1352, 1334, 1328, 1330, 1335, 1378, ...
        1433, 1540, 1616, 1731, 1810, 1796, 1732, 1600, 1484, 1432, 1381, 1328, ...
        1263, 1263, 1277, 1256, 1261, 1312, 1442, 1461, 1585, 1698, 1890, 2030, ...
        2139, 2146, 2132, 2074, 2036, 1963, 1899, 1815, 1726, 1669, 1603, 1665], ...
        [1027, 1014, 1014, 1017, 1006, 1010, 1008, 1018, 1048, 1043, 1108, 1195, ...
        1286, 1361, 1451, 1499, 1494, 1479, 1441, 1451, 1444, 1411, 1376, 1353, ...
        1336, 1308, 1310, 1291, 1371, 1254, 1295, 1367, 1394, 1437, 1465, 1493, ...
        1479, 1458, 1453, 1411, 1408, 1371, 1335, 1289, 1242, 1211, 1181, 1160], ...
        [4586, 4421, 4280, 4117, 3999, 3936, 3908, 3902, 3913, 4034, 4211, 4400, ...
        4900, 5383, 5872, 6246, 6340, 6398, 6410, 6262, 6092, 5911, 5731, 5638, ...
        5481, 5395, 5274, 5059, 5011, 5184, 5444, 5775, 6208, 6492, 6797, 7000, ...
        7002, 6896, 6771, 6648, 6434, 6222, 6017, 5693, 5356, 5157, 5119, 5046]
    };

    tariff_prices = {
        [99.34, 95.51, 88.41, 87.9, 82.77, 80.29, 78.36, 75.29, 81.26, 87, 92.83, ...
        91.46, 122.65, 172.54, 153.32, 128.48, 92.38, 76.01, 79.65, 72.44, 63.66, ...
        60.05, 56.5, 59.83, 56.57, 62.71, 74.45, 74.71, 76.89, 76.76, 76.67, 102.19, ...
        151.03, 175.66, 600, 600, 377.38, 196.43, 154.77, 162.63, 154.79, 138.49, ...
        126.93, 103.93, 109.45, 108.73, 107.64, 101.32], ...
        [94.53, 90.69, 83.43, 83.05, 78.61, 76.9, 75.27, 73.87, 79.05, 85.34, 95.44, ...
        93.62, 127.56, 310.08, 153.59, 112.18, 74.29, 52.81, 50.76, 40.57, 29.59, ...
        34.26, 26.64, 31.58, 32.52, 40.22, 43.79, 48.13, 57.44, 64.17, 66.33, 91.2, ...
        135.77, 163.53, 600, 600, 372.43, 194.91, 151.28, 165.32, 154.06, 131.56, ...
        121.15, 99.05, 108.73, 103.86, 105.46, 95.61], ...
        [66.86, 68.47, 56.32, 49.55, 44.23, 34.71, 26.7, 28.44, 30.28, 40.64, 46.87, ...
        55.92, 81.63, 107, 114.65, 139.21, 132.13, 118.9, 118.56, 100.42, 83.68, ...
        66.95, 62.82, 49.65, 48.9, 53.92, 56.11, 48.89, 63.45, 67.74, 74.74, 96.26, ...
        126.48, 142.12, 166.15, 200.62, 159.67, 156.59, 134.51, 130.49, 122.16, ...
        115.26, 103.48, 93, 85.25, 83.99, 77.98, 206.72], ...
        [22.61, 20.53, 15.97, 13.03, 13.09, 13.15, 19.57, 18.02, 23.94, 23.83, 27.29, ...
        33.96, 54.96, 72.37, 65.25, 69.1, 58.78, 54.69, 64.42, 60.94, 50.5, 48.35, ...
        47.19, 33.31, 34.7, 29.31, 32.57, 36.13, 48.62, 53.47, 177.81, 58.13, ...
        70.26, 69.25, 65.75, 53.51, 48.78, 56.82, 53.92, 57.97, 55.19, 52.02, 48.22, ...
        38.79, 33.74, 33.15, 30.72, 29.31], ...
        [98.12, 93.17, 85.48, 84.82, 80.07, 76.62, 74.93, 73.62, 78.18, 83.64, ...
        88.69, 89.19, 119.49, 160.82, 143.17, 121.11, 86.37, 70.31, 74.08, 67.7, ...
        58.03, 53.71, 50.24, 52.76, 48.72, 54.35, 66.5, 67.68, 70.86, 72.13, 71.96, ...
        95.4, 142.51, 166.38, 600, 600, 377.04, 196.59, 153.62, 161.65, 151.72, ...
        135.17, 124.55, 101.33, 106.77, 106.22, 105.48, 98.56]
    };

    locations = {'NSW', 'QLD', 'SA', 'TAS', 'VIC'};
    colors = {'b', 'r', 'g', 'm', 'c'};

    all_P_batt = [];
    all_P_grid = [];
    all_soc_history = [];
    all_metrics = struct('energy_cost', [], 'peak_to_average', [], 'reliability_cost', [], 'resiliency_cost', [], 'Battery_degradation', [], 'Battery_automation', [], 'power_imbalance_penalty', [], 'battery_capacity_penalty', []);
    final_objectives = zeros(1, length(locations));
    all_fitness_history = cell(1, length(locations));  % To store fitness history of all locations

    for i = 1:length(locations)
        P_d = load_profiles{i};
        Cg_t = tariff_prices{i};

        % Normalize P_d and Cg_t
        P_d_max = max(P_d);
        Cg_t_max = max(Cg_t);

        P_d = P_d / P_d_max;
        Cg_t = Cg_t / Cg_t_max;

        xmax = 0.9 * ones(1, 2 * T);
        xmin = [-0.9 * ones(1, T), zeros(1, T)];
        popsize = 500;
        dimension = 2 * T;
        MaxFEs = 1000000;

        % Execute the Growth Optimizer
        [best_solution, best_fitness, fitness_history, metrics_history] = GO(popsize, dimension, xmax, xmin, MaxFEs, P_d, Cg_t, params);

        all_fitness_history{i} = fitness_history;  % Store fitness history for current location

        % Calculate the battery state of charge (SoC) in percentage
        P_batt = best_solution(1:2:end);
        P_grid = best_solution(2:2:end);
        battery_state = params.battery_capacity;
        soc_history = zeros(1, length(P_batt));
        for j = 1:length(P_batt)
            battery_state = battery_state - P_batt(j);
            soc_history(j) = (battery_state / params.battery_capacity) * 100;
        end

        all_P_batt = [all_P_batt; P_batt];
        all_P_grid = [all_P_grid; P_grid];
        all_soc_history = [all_soc_history; soc_history];

        % Store the final metrics
        final_metrics = metrics_history(end);
        all_metrics.energy_cost = [all_metrics.energy_cost; final_metrics.energy_cost];
        all_metrics.peak_to_average = [all_metrics.peak_to_average; final_metrics.peak_to_average];
        all_metrics.reliability_cost = [all_metrics.reliability_cost; final_metrics.reliability_cost];
        all_metrics.resiliency_cost = [all_metrics.resiliency_cost; final_metrics.resiliency_cost];
        all_metrics.Battery_degradation = [all_metrics.Battery_degradation; final_metrics.Battery_degradation];
        all_metrics.Battery_automation = [all_metrics.Battery_automation; final_metrics.Battery_automation];
        all_metrics.power_imbalance_penalty = [all_metrics.power_imbalance_penalty; final_metrics.power_imbalance_penalty];
        all_metrics.battery_capacity_penalty = [all_metrics.battery_capacity_penalty; final_metrics.battery_capacity_penalty];

        % Store the final objective value
        final_objectives(i) = best_fitness;
    end

    % Plot the learning curves of each location in one plot with thicker lines
    figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
    hold on;
    for i = 1:length(locations)
        plot(all_fitness_history{i}, 'Color', colors{i}, 'LineWidth', 2);
    end
    title('Learning Curves for Each Location', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    xlabel('Function Evaluations', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('Fitness', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
    legend(locations, 'FontName', 'Times New Roman', 'FontSize', 20, 'Location', 'best');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
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
            battery_capacity_penalties(i) = 5 * abs(battery_state -  params.battery_capacity);
        elseif battery_state > params.battery_capacity
            % Apply penalty for overcharge (exceeding the battery capacity)
            battery_capacity_penalties(i) = 5 * (battery_state - params.battery_capacity);
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
end

function y = safe_exp(x)
    y = exp(min(x, 20));  % Prevent overflow
end

% Growth Optimizer: A powerful metaheuristic algorithm for solving different optimization problems
function [gbestX, gbestfitness, gbesthistory, metrics_history] = GO(popsize, dimension, xmax, xmin, MaxFEs, P_d, Cg_t, params)
    FEs = 0;
    fitness = zeros(1, popsize);
    x = xmin + (xmax - xmin) .* rand(popsize, dimension);
    gbestfitness = inf;
    gbesthistory = [];
    metrics_history = struct('energy_cost', [], 'peak_to_average', [], 'reliability_cost', [], 'resiliency_cost', [], 'Battery_degradation', [], 'Battery_automation', [], 'power_imbalance_penalty', [], 'battery_capacity_penalty', []);

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
            Worst_X = x(ind(randi([popsize - 4, popsize])), :);
            Better_X = x(ind(randi([2, 5])), :);
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
