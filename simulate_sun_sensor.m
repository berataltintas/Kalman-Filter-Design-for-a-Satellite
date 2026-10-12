function [measurement_body, true_body, valid] = simulate_sun_sensor(...
        sun_eci, attitude_body_eci, params, in_eclipse)
    true_body = attitude_body_eci * sun_eci;
    true_sensor = params.sunsensor.mounting_dcm * true_body;
    valid = ~in_eclipse && ...
        acos(max(-1, min(1, true_sensor(1)))) <= params.sunsensor.fov / 2;

    if ~valid
        measurement_body = nan(3, 1);
        return;
    end

    noisy_sensor = params.sunsensor.misalignment * true_sensor + ...
        chol(params.sunsensor.covariance, 'lower') * randn(3, 1);
    noisy_sensor = noisy_sensor / norm(noisy_sensor);
    measurement_body = params.sunsensor.mounting_dcm' * noisy_sensor;
    measurement_body = measurement_body / norm(measurement_body);
end