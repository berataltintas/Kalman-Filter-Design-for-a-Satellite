% Simulation parameters. Sensor samples are generated in main.m.
params.orbit.earth_radius = 6371e3;
params.orbit.a = params.orbit.earth_radius + 505e3;
params.orbit.e = 0;
params.orbit.i = deg2rad(97.43);
% At the epoch, the ascending node is at 18:00 local solar time.
params.orbit.raan = 0;
params.orbit.arg_perigee = 0;
params.orbit.mean_anomaly = 0;

params.simulation.dt = 1;
mu = 3.986004418e14;
params.simulation.orbit_period = 2 * pi * sqrt(params.orbit.a^3 / mu);
params.simulation.duration = params.simulation.orbit_period;
params.simulation.mode = 'fine';

% The DCM maps body-frame vectors into each sensor frame.
params.sunsensor.fov = deg2rad(90);
params.sunsensor.accuracy = deg2rad(1.5);
params.sunsensor.accuracy_sigma = params.sunsensor.accuracy / 3;
params.sunsensor.mounting_dcm = eye(3);
params.sunsensor.misalignment = small_angle_dcm(deg2rad(0.16));
params.sunsensor.covariance = params.sunsensor.accuracy_sigma^2 * eye(3);

params.magnetometer.accuracy = 8.5e-9;
params.magnetometer.bias = 350e-9;
params.magnetometer.mounting_dcm = eye(3);
params.magnetometer.misalignment = small_angle_dcm(deg2rad(1));
params.magnetometer.bias_vector = params.magnetometer.bias * ones(3, 1);
params.magnetometer.covariance = params.magnetometer.accuracy^2 * eye(3);

params.magnetic_field.surface_strength = 3.12e-5;
params.magnetic_field.axis = [0; 0; 1];
params.sun.epoch_eci_direction = [0; -1; 0];
params.sun.annual_rate = 2 * pi / (365.256363004 * 86400);
params.sun.eci_direction = sun_vector(0, params);

params.startracker.accuracy_xy = deg2rad(4.1 / 3 / 3600);
params.startracker.accuracy_z = deg2rad(32 / 3 / 3600);
params.startracker.mounting_dcm = eye(3);
params.startracker.covariance = diag([params.startracker.accuracy_xy^2, ...
    params.startracker.accuracy_xy^2, params.startracker.accuracy_z^2]);

params.gyro.arw = deg2rad(0.005) / sqrt(3600) *100 ;
params.gyro.bias_sigma = deg2rad(0.01) / 3600;
params.gyro.bias_random_walk = params.gyro.bias_sigma / sqrt(3600);
params.gyro.misalignment = small_angle_dcm(500e-6 / sqrt(3));
params.gyro.covariance = params.gyro.arw^2 * eye(3);

% Coarse-mode attitude error and normalized vector measurement covariances.
params.coarse.process_covariance = deg2rad(0.1)^2 * eye(3);
params.coarse.sun_covariance = params.sunsensor.accuracy_sigma^2 * eye(3);
params.coarse.attitude_covariance = deg2rad(1)^2 * eye(3);
magnetic_direction_sigma = params.magnetometer.accuracy / ...
    (params.magnetic_field.surface_strength * ...
    (params.orbit.earth_radius / params.orbit.a)^3);
params.coarse.magnetic_covariance = magnetic_direction_sigma^2 * eye(3);

function dcm = small_angle_dcm(angle)
    skew = [0, -angle, angle; angle, 0, -angle; -angle, angle, 0];
    dcm = eye(3) - skew;
end