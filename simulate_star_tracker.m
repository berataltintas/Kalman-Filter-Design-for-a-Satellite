function [measurement_body_eci, attitude_error] = simulate_star_tracker(...
        attitude_body_eci, params)
    true_sensor_eci = params.startracker.mounting_dcm * attitude_body_eci;
    noise = chol(params.startracker.covariance, 'lower') * randn(3, 1);
    noisy_sensor_eci = small_angle_dcm(noise) * true_sensor_eci;
    measurement_body_eci = params.startracker.mounting_dcm' * noisy_sensor_eci;
    attitude_error = acos(max(-1, min(1, (trace(...
        measurement_body_eci * attitude_body_eci') - 1) / 2)));
end

function dcm = small_angle_dcm(vector)
    angle = norm(vector);
    if angle < eps
        dcm = eye(3);
        return;
    end

    axis = vector / angle;
    skew = [0, -axis(3), axis(2); ...
        axis(3), 0, -axis(1); ...
        -axis(2), axis(1), 0];
    dcm = eye(3) - sin(angle) * skew + ...
        (1 - cos(angle)) * skew^2;
end