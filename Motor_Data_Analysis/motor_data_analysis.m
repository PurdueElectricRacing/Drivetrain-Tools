% Set up the Import Options and import the data
num_sheets = 29;
opts = spreadsheetImportOptions("NumVariables", 7);
motor_data = zeros(69, 7, 29);

% Specify sheet and range
opts.DataRange = "A2:G70";

% Specify column names and types
opts.VariableNames = ["VoltageV", "CurrentA", "SpeedRPM", "InputPowerW", "OutputPowerW", "TorqueNcm", "Efficiency"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% Import the data
for i = 1:1:29
    voltage = (i*10) + 50;
    if voltage ~= 200
        opts.Sheet = num2str(voltage) + "V";
        motor_data(:, :, i) = table2array(readtable("all_motor_data.xlsx", opts, "UseExcel", false));
    end
end

% Clear temporary variables
clear opts


%% Motor Data Analysis setup constants
zero_torque_efficiency = 0.4;
zero_rpm_efficiency = 0.4;
rpm2radps = 0.104719755;
continuous_power_loss_limit = 1800; 
num_datasets = 29;

max_torque = 0;
max_torques = zeros(1, num_datasets);
max_rpm = zeros(1, num_datasets);
voltages = (60:10:340);
rpm_no_load = [1901.4 2225.1 2545.5 2858.5 3168.8 3481.0 3776.5 4103.2 4416.2 4735.6 5053.0 5358.8 5648.8 5968.3 0 6597.0 6895.6 7204.3 7508.2 7816.1 8132.9 8445.9 8747.4 9060.4 9288.4 9606.0 9907.0 10204.8 10531.2];
current_no_load = [2.8 2.7 2.9 2.8 2.9 2.9 3 2.9 3 3.1 3.1 3.2 3.3 3.3 0 3.2 3.3 3.4 3.4 3.4 3.3 3.5 3.5 3.5 3.4 3.4 3.4 3.4 3.5];
all_torque = [];
all_rpm = [];
all_efficiency = [];
all_voltage = [];

% Iterate through each data set to construct dataset of torque, rpm, efficiency
for i = 1:1:num_datasets
    if voltages(i) ~= 200
        rpm_data = motor_data(:, 3, i);
        torque_data = motor_data(:, 6, i);
        max_rpm(i) = pchip(torque_data, rpm_data, 2500);
        valid_region = max_rpm(i) < rpm_data;
        valid_rpm = [rpm_no_load(i)+voltages(i); rpm_data(valid_region); max_rpm(i)];
        valid_torque = [0; torque_data(valid_region); 2500];
        num_elem = length(valid_torque);
        efficiency_data = [zero_torque_efficiency*100; motor_data(2:num_elem, 7, i)] ./ 100;
    else
        current_range = [0.1 (5:1:30) (32:1:70) 72]';
        valid_rpm = valid_rpm + (10*31.85);
        p_in = current_range .* 200;
        p_out = valid_rpm .* valid_torque .* 0.104719755 ./ 100;
        efficiency_data = p_out ./ p_in;

       if efficiency_data(1) == 0
            efficiency_data(1) = zero_torque_efficiency;
        end
    end

    all_torque = [all_torque; valid_torque];
    all_rpm = [all_rpm; valid_rpm];
    all_efficiency = [all_efficiency; efficiency_data];
    all_voltage = [all_voltage; voltages(i) .* ones(length(valid_torque), 1)];

    max_torque = max_torque + max(torque_data);
    max_torques(i) = max(torque_data);
end

% Final contruction of torque, rpm, efficiency dataset
all_rpm = [all_rpm; 0; max(all_rpm); max(all_rpm); 0];
all_torque = [all_torque; 0; 2500; 0; 2500];
all_efficiency = [all_efficiency; zero_rpm_efficiency; 0.9; zero_torque_efficiency; zero_rpm_efficiency];
all_voltage = [all_voltage; 0; 500; 340; 0];


%% Convert Dataset into a nice, evenly spaced values
Tx_resolution = 161;
RPM_resolution = 107;

rpm_sweep = linspace(0, 10600, RPM_resolution);
torque_sweep = linspace(0, 2500, Tx_resolution);

rpm_extend = ones(Tx_resolution, 1);
torque_extend = ones(1, RPM_resolution);

rpm_grid = rpm_extend * rpm_sweep;
torque_grid = torque_sweep' * torque_extend;

%rpm_grid = rpm_grid(:);
%torque_grid = torque_grid(:);

efficiency_grid = (griddata(all_rpm,all_torque,all_efficiency,rpm_grid,torque_grid)') .* 100;
voltage_grid = (griddata(all_rpm,all_torque,all_voltage,rpm_grid,torque_grid)')';
voltage_grid(:,1) = (torque_sweep ./ 25);

max_torque = flip([all_torque(1874:end-4); 2500]');
max_torque_double_sided = [-flip(max_torque(1:end-1)) max_torque];

max_rpm = flip([all_rpm(1874:end-4); 0]');
max_rpm_double_sided = [flip(max_rpm(1:end-1)) max_rpm];


%% Unit Conversions
rpm_sweep = rpm_sweep .* rpm2radps;
torque_sweep = torque_sweep ./ 100;

max_torque = max_torque ./ 100;
max_rpm = max_rpm .* rpm2radps;

rpm_grid = rpm_grid .* rpm2radps;
torque_grid = torque_grid ./ 100;

A = unique(all_rpm);
efficiency_grid_temp  = efficiency_grid' ./ 100;

figure(10)
scatter3(rpm_grid, torque_grid, efficiency_grid_temp)
figure(20)
scatter3(rpm_grid, torque_grid, voltage_grid)

x1 = rpm_grid(:);
y1 = torque_grid(:);
z1 = efficiency_grid_temp(:);


%% Fitting
[xData, yData, zData] = prepareSurfaceData( x1, y1, z1 );

% Set up fittype and options.
ft = fittype( 'loess' );

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'untitled fit 1', 'z1 vs. x1, y1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x1', 'Interpreter', 'none' );
ylabel( 'y1', 'Interpreter', 'none' );
zlabel( 'z1', 'Interpreter', 'none' );
grid on
view( 9.5, 23.6 );

z2 = feval(fitresult,[x1,y1]);

efficiency_grid = reshape(z2, [Tx_resolution, RPM_resolution]);

%% Find Torque-RPM curve for particular power loss limit
power_loss_grid = abs(torque_grid .* rpm_grid .* ((1 ./efficiency_grid) - 1));
power_loss_grid(:,1) = (torque_sweep ./ 25);
power_input_grid = abs(torque_grid .* rpm_grid ./ efficiency_grid);
power_input_grid(:,1) = (torque_sweep ./ 25);
valid_power_limit = power_loss_grid <= continuous_power_loss_limit;

power_loss_continuous = power_loss_grid(valid_power_limit);
rpm_continuous = rpm_grid(valid_power_limit);
torque_continuous = torque_grid(valid_power_limit);

figure(1)
scatter(max_rpm, max_torque)

figure(2)
scatter3(all_rpm, all_torque, all_efficiency)

figure(3)
scatter3(rpm_grid(:), torque_grid(:), efficiency_grid_temp(:))

figure(6)
scatter3(rpm_grid(:), torque_grid(:), power_loss_grid(:));

figure(7)
scatter(rpm_continuous, torque_continuous)
hold on
scatter(max_rpm, max_torque)

figure(20)
scatter3(x1, y1, z2)