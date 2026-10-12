clear;
close all;
clc;
rng(7);

init;
coarse = simulate_mode(params, 'coarse');
fine = simulate_mode(params, 'fine');

coarse_q_estimated = zeros(4, numel(coarse.time));
coarse_rpy_estimated = zeros(3, numel(coarse.time));
coarse_covariance = deg2rad(5)^2 * eye(3);
first_valid_sun = find(coarse.sun_valid, 1, 'first');
sun_eci = sun_vector(coarse.time(first_valid_sun), params);
coarse_q_estimated(:, first_valid_sun) = coarse_triad(...
    coarse.sun_measured_body(:, first_valid_sun), ...
    coarse.magnetic_measured_body(:, first_valid_sun), ...
    coarse.position_eci(:, first_valid_sun), sun_eci, params);
[~, coarse_initial_dcm] = kinematics(...
    coarse_q_estimated(:, first_valid_sun), zeros(3, 1), 0);
coarse_rpy_estimated(:, first_valid_sun) = dcm_to_rpy(coarse_initial_dcm);
for sample = first_valid_sun + 1:numel(coarse.time)
    sun_eci = sun_vector(coarse.time(sample), params);
    [coarse_q_estimated(:, sample), coarse_covariance] = coarse_mekf(...
        coarse.sun_measured_body(:, sample), ...
        coarse.magnetic_measured_body(:, sample), ...
        coarse.position_eci(:, sample), sun_eci, ...
        coarse_q_estimated(:, sample - 1), ...
        coarse_covariance, params.simulation.dt, params);
    [~, coarse_dcm] = kinematics(...
        coarse_q_estimated(:, sample), zeros(3, 1), 0);
    coarse_rpy_estimated(:, sample) = dcm_to_rpy(coarse_dcm);
end

% Eclipse mode: initialize once, then use only gyro and magnetometer.
eclipse_q_estimated = zeros(4, numel(coarse.time));
eclipse_bias_estimated = zeros(3, numel(coarse.time));
eclipse_rpy_estimated = zeros(3, numel(coarse.time));
eclipse_initial_attitude_sigma = deg2rad(5);
eclipse_initial_bias_sigma = deg2rad(0.01) / 3600;
eclipse_covariance = diag([eclipse_initial_attitude_sigma^2 * ones(1, 3), ...
    eclipse_initial_bias_sigma^2 * ones(1, 3)]);
eclipse_start = first_valid_sun;
sun_eci = sun_vector(coarse.time(eclipse_start), params);
eclipse_q_estimated(:, eclipse_start) = coarse_triad(...
    coarse.sun_measured_body(:, eclipse_start), ...
    coarse.magnetic_measured_body(:, eclipse_start), ...
    coarse.position_eci(:, eclipse_start), sun_eci, params);
[~, eclipse_initial_dcm] = kinematics(...
    eclipse_q_estimated(:, eclipse_start), zeros(3, 1), 0);
eclipse_rpy_estimated(:, eclipse_start) = dcm_to_rpy(eclipse_initial_dcm);

for sample = eclipse_start + 1:numel(coarse.time)
    [eclipse_q_estimated(:, sample), eclipse_bias_estimated(:, sample), ...
        eclipse_covariance] = eclipse_mekf(...
        coarse.gyro_measurement(:, sample), ...
        coarse.magnetic_measured_body(:, sample), ...
        coarse.position_eci(:, sample), eclipse_q_estimated(:, sample - 1), ...
        eclipse_bias_estimated(:, sample - 1), eclipse_covariance, ...
        params.simulation.dt, params);
    [~, eclipse_dcm] = kinematics(...
        eclipse_q_estimated(:, sample), zeros(3, 1), 0);
    eclipse_rpy_estimated(:, sample) = dcm_to_rpy(eclipse_dcm);
end

sample_count = numel(fine.time);
true_quaternion = zeros(4, sample_count);
for sample = 1:sample_count
    current_dcm = fine.true_attitude_body_eci(:, :, sample);
    true_quaternion(:, sample) = dcm_to_quaternion(current_dcm);
end

q_estimated = zeros(4, sample_count);
bias_estimated = zeros(3, sample_count);
estimated_rpy = zeros(3, sample_count);
attitude_initial_sigma = deg2rad(1);
bias_initial_sigma = deg2rad(0.01) / 3600;
P_estimated = diag([attitude_initial_sigma^2 * ones(1, 3), ...
    bias_initial_sigma^2 * ones(1, 3)]);

q_estimated(:, 1) = dcm_to_quaternion(fine.star_tracker_measurement(:, :, 1));
[~, initial_estimated_dcm] = kinematics(q_estimated(:, 1), zeros(3, 1), 0);
estimated_rpy(:, 1) = dcm_to_rpy(initial_estimated_dcm);
bias_estimated(:, 1) = zeros(3, 1);

for sample = 2:sample_count
    [q_estimated(:, sample), bias_estimated(:, sample), P_estimated] = mekf( ...
        fine.gyro_measurement(:, sample), ...
        fine.star_tracker_measurement(:, :, sample), ...
        q_estimated(:, sample - 1), ...
        bias_estimated(:, sample - 1), ...
        P_estimated, params.simulation.dt, params);

    [~, estimated_dcm] = kinematics(q_estimated(:, sample), zeros(3, 1), 0);
    estimated_rpy(:, sample) = dcm_to_rpy(estimated_dcm);
end

aligned_fine_quaternion = align_quaternion_signs(...
    true_quaternion, q_estimated);
fine_quaternion_rmse = sqrt(mean((aligned_fine_quaternion - ...
    true_quaternion).^2, 2));
fine_attitude_error = zeros(1, sample_count);
for sample = 1:sample_count
    [~, estimated_dcm] = kinematics(...
        q_estimated(:, sample), zeros(3, 1), 0);
    trace_error = (trace(estimated_dcm * ...
        fine.true_attitude_body_eci(:, :, sample)') - 1) / 2;
    fine_attitude_error(sample) = acos(max(-1, min(1, trace_error)));
end
fprintf('Fine quaternion component RMSE: %.6g %.6g %.6g %.6g\n', ...
    fine_quaternion_rmse);
fprintf('Fine attitude error RMS/max (deg): %.6g %.6g\n', ...
    rad2deg(rms(fine_attitude_error)), rad2deg(max(fine_attitude_error)));



radius = vecnorm(coarse.position_eci);
altitude = radius - params.orbit.earth_radius;
speed = vecnorm(coarse.velocity_eci);
fprintf('Coarse Mode and Fine Mode simulated separately.\n');
fprintf('Orbit period: %.1f s\n', params.simulation.orbit_period);
fprintf('Radius variation: %.3f m\n', max(radius) - min(radius));
fprintf('Altitude: %.3f to %.3f km\n', min(altitude) / 1e3, max(altitude) / 1e3);
fprintf('Speed: %.3f to %.3f km/s\n', min(speed) / 1e3, max(speed) / 1e3);
fprintf('Coarse eclipse samples: %d / %d\n', ...
    nnz(coarse.eclipse), numel(coarse.time));
fprintf('Fine eclipse samples: %d / %d\n', ...
    nnz(fine.eclipse), numel(fine.time));

figure('Name', 'Coarse Mode Attitude');
subplot(3, 1, 1);
plot(coarse.time, rad2deg(coarse.true_rpy(1, :)), 'LineWidth', 1.1);
hold on;
plot(coarse.time, rad2deg(coarse_rpy_estimated(1, :)), '--', 'LineWidth', 1.1);
grid on;
ylabel('Roll (deg)');
title('Coarse MEKF Roll');
legend('True', 'Coarse MEKF');

subplot(3, 1, 2);
plot(coarse.time, rad2deg(coarse.true_rpy(2, :)), 'LineWidth', 1.1);
hold on;
plot(coarse.time, rad2deg(coarse_rpy_estimated(2, :)), '--', 'LineWidth', 1.1);
grid on;
ylabel('Pitch (deg)');
title('Coarse MEKF Pitch');
legend('True', 'Coarse MEKF');

subplot(3, 1, 3);
plot(coarse.time, rad2deg(coarse.true_rpy(3, :)), 'LineWidth', 1.1);
hold on;
plot(coarse.time, rad2deg(coarse_rpy_estimated(3, :)), '--', 'LineWidth', 1.1);
grid on;
xlabel('Time (s)');
ylabel('Yaw (deg)');
title('Coarse MEKF Yaw');
legend('True', 'Coarse MEKF');

figure('Name', 'Coarse Eclipse Attitude');
subplot(3, 1, 1);
plot(coarse.time, rad2deg(coarse.true_rpy(1, :)), 'LineWidth', 1.1);
hold on;
plot(coarse.time, rad2deg(eclipse_rpy_estimated(1, :)), '--', 'LineWidth', 1.1);
grid on;
ylabel('Roll (deg)');
title('Eclipse Mode: Gyro + Magnetometer Roll');
legend('True', 'Estimate');

subplot(3, 1, 2);
plot(coarse.time, rad2deg(coarse.true_rpy(2, :)), 'LineWidth', 1.1);
hold on;
plot(coarse.time, rad2deg(eclipse_rpy_estimated(2, :)), '--', 'LineWidth', 1.1);
grid on;
ylabel('Pitch (deg)');
title('Eclipse Mode: Gyro + Magnetometer Pitch');
legend('True', 'Estimate');

subplot(3, 1, 3);
plot(coarse.time, rad2deg(coarse.true_rpy(3, :)), 'LineWidth', 1.1);
hold on;
plot(coarse.time, rad2deg(eclipse_rpy_estimated(3, :)), '--', 'LineWidth', 1.1);
grid on;
xlabel('Time (s)');
ylabel('Yaw (deg)');
title('Eclipse Mode: Gyro + Magnetometer Yaw');
legend('True', 'Estimate');

figure('Name', 'Propagated Orbit');
plot3(coarse.position_eci(1, :) / 1e3, coarse.position_eci(2, :) / 1e3, ...
    coarse.position_eci(3, :) / 1e3, 'LineWidth', 1.2);
grid on;
axis equal;
xlabel('ECI X (km)'); ylabel('ECI Y (km)'); zlabel('ECI Z (km)');
title('Propagated Circular Orbit');

figure('Name', 'Coarse Mode Sun Sensor Measurements');
sun_axis_labels = {'X', 'Y', 'Z'};
for axis_index = 1:3
    subplot(3, 1, axis_index);
    plot(coarse.time, coarse.sun_true_body(axis_index, :), 'k', ...
        'LineWidth', 1.2);
    hold on;
    plot(coarse.time, coarse.sun_measured_body(axis_index, :), 'r', ...
        'LineWidth', 1.0);
    grid on;
    ylabel(['Sun vector ' sun_axis_labels{axis_index}]);
    legend('True', 'Measured');
end
xlabel('Time (s)');
sgtitle('Coarse Mode: Sun Sensor, Body Frame');

figure('Name', 'Coarse Mode Magnetometer Measurements');
plot(coarse.time, coarse.magnetic_true_body * 1e9, 'LineWidth', 1.2);
hold on;
plot(coarse.time, coarse.magnetic_measured_body * 1e9, '.');
grid on;
xlabel('Time (s)');
ylabel('Magnetic field (nT)');
title('Coarse Mode: Magnetometer, Body Frame');
legend('True X', 'True Y', 'True Z', 'Measured X', 'Measured Y', 'Measured Z');

figure('Name', 'Fine Mode Sensor Measurements');
subplot(2, 1, 1);
plot(fine.time, rad2deg(fine.star_tracker_error) * 3600, 'k');
grid on;
ylabel('Attitude error (arcsec)');
title('Fine Mode: Star Tracker Attitude Error');

figure('Name', 'Fine Mode Euler Attitude');
subplot(3, 1, 1);
plot(fine.time, rad2deg(fine.true_rpy(1, :)), 'LineWidth', 1.1);
hold on;
plot(fine.time, rad2deg(estimated_rpy(1, :)), '--', 'LineWidth', 1.1);
grid on;
xlabel('Time (s)');
ylabel('Roll (deg)');
title('Fine Mode: True vs MEKF Estimated Roll');
legend('True', 'MEKF');

subplot(3, 1, 2);
plot(fine.time, rad2deg(fine.true_rpy(2, :)), 'LineWidth', 1.1);
hold on;
plot(fine.time, rad2deg(estimated_rpy(2, :)), '--', 'LineWidth', 1.1);
grid on;
xlabel('Time (s)');
ylabel('Pitch (deg)');
title('Fine Mode: True vs MEKF Estimated Pitch');
legend('True', 'MEKF');

subplot(3, 1, 3);
plot(fine.time, rad2deg(fine.true_rpy(3, :)), 'LineWidth', 1.1);
hold on;
plot(fine.time, rad2deg(estimated_rpy(3, :)), '--', 'LineWidth', 1.1);
grid on;
xlabel('Time (s)');
ylabel('Yaw (deg)');
title('Fine Mode: True vs MEKF Estimated Yaw');
legend('True', 'MEKF');

figure('Name', 'Fine Mode Gyro Measurements');
subplot(2, 1, 1);
plot(fine.time, fine.gyro_measurement * 180 / pi * 3600, 'LineWidth', 1);
hold on;
plot(fine.time, fine.gyro_bias * 180 / pi * 3600, '--');
plot(fine.time, fine.true_rate_body * 180 / pi * 3600, ':', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Angular rate (arcsec/s)');
title('Fine Mode: Gyro Measurement and Bias');
legend('Measured X', 'Measured Y', 'Measured Z', ...
    'Bias X', 'Bias Y', 'Bias Z', 'True X', 'True Y', 'True Z');

gyro_bias_deg_per_hour = fine.gyro_bias * 180 / pi * 3600;
fprintf('Gyro bias max abs (deg/h): %.6g %.6g %.6g\n', ...
    max(abs(gyro_bias_deg_per_hour), [], 2));

figure('Name', 'Fine Mode Gyro Bias');
plot(fine.time, gyro_bias_deg_per_hour, 'LineWidth', 1.1);
grid on;
xlabel('Time (s)');
ylabel('Bias (deg/h)');
title('Fine Mode: Estimated Sensor Bias Used in Simulation');
legend('Bias X', 'Bias Y', 'Bias Z');
    figure('Name', 'True Attitude Quaternion');

subplot(4, 1, 1);
plot(fine.time, true_quaternion(1, :), 'LineWidth', 1.1);
grid on;
ylabel('q_0');
title('True Attitude Quaternion');

subplot(4, 1, 2);
plot(fine.time, true_quaternion(2, :), 'LineWidth', 1.1);
grid on;
ylabel('q_1');

subplot(4, 1, 3);
plot(fine.time, true_quaternion(3, :), 'LineWidth', 1.1);
grid on;
ylabel('q_2');

subplot(4, 1, 4);
plot(fine.time, true_quaternion(4, :), 'LineWidth', 1.1);
grid on;
xlabel('Time (s)');
ylabel('q_3');

plot_quaternion_comparison(...
    coarse.time, true_quaternion, coarse_q_estimated, ...
    'Coarse Mode Quaternion: True vs Estimate', 'Coarse MEKF');
plot_quaternion_comparison(...
    coarse.time, true_quaternion, eclipse_q_estimated, ...
    'Eclipse Mode Quaternion: True vs Estimate', ...
    'Gyro + Magnetometer MEKF');
plot_quaternion_comparison(...
    fine.time, true_quaternion, q_estimated, ...
    'Fine Mode Quaternion: True vs Estimate', 'Fine MEKF');

function plot_quaternion_comparison(time, true_quaternion, ...
        estimated_quaternion, figure_title, estimate_label)
    estimated_quaternion = align_quaternion_signs(...
        true_quaternion, estimated_quaternion);
    component_names = {'q_0', 'q_1', 'q_2', 'q_3'};

    figure('Name', figure_title);
    for component = 1:4
        subplot(4, 1, component);
        plot(time, true_quaternion(component, :), 'LineWidth', 1.1);
        hold on;
        plot(time, estimated_quaternion(component, :), '--', ...
            'LineWidth', 1.1);
        grid on;
        ylabel(component_names{component});
        if component == 1
            title(figure_title);
        end
        if component == 4
            xlabel('Time (s)');
        end
        legend('True', estimate_label);
    end
end

function aligned_quaternion = align_quaternion_signs(...
        true_quaternion, estimated_quaternion)
    aligned_quaternion = estimated_quaternion;
    for sample = 1:size(estimated_quaternion, 2)
        if norm(estimated_quaternion(:, sample)) > 0 && ...
                dot(true_quaternion(:, sample), ...
                estimated_quaternion(:, sample)) < 0
            aligned_quaternion(:, sample) = -estimated_quaternion(:, sample);
        end
    end
end