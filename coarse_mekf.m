function [q_hat, P_hat] = coarse_mekf(...
        sun_measurement_body, magnetic_measurement_body, position_eci, ...
    sun_eci, q_prev, P_prev, dt, params)
    % Coarse mode uses calibrated sun and magnetic direction measurements.
    sun_body = calibrate_sun_measurement(sun_measurement_body, params);
    magnetic_body = calibrate_magnetic_measurement(...
        magnetic_measurement_body, params);

    P_pred = P_prev + params.coarse.process_covariance * dt;

    if all(isfinite(sun_body))
        q_measurement = triad_measurement(...
            sun_body, magnetic_body, position_eci, sun_eci, params);
        if dot(q_measurement, q_prev) < 0
            q_measurement = -q_measurement;
        end

        q_error = quaternion_multiply(q_measurement, quaternion_conjugate(q_prev));
        if q_error(1) < 0
            q_error = -q_error;
        end
        innovation = 2 * q_error(2:4);
    else
        innovation = zeros(3, 1);
    end

    H = eye(3);
    R = params.coarse.attitude_covariance;
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    correction = K * innovation;

    delta_q = [1; 0.5 * correction];
    delta_q = delta_q / norm(delta_q);
    q_hat = quaternion_multiply(delta_q, q_prev);
    q_hat = q_hat / norm(q_hat);

    identity = eye(3);
    P_hat = (identity - K * H) * P_pred * (identity - K * H)' + ...
        K * R * K';
    P_hat = 0.5 * (P_hat + P_hat');
end

function q_measurement = triad_measurement(...
        sun_body, magnetic_body, position_eci, sun_eci, params)
    sun_eci = sun_eci / norm(sun_eci);
    magnetic_eci = get_magnetic_field(position_eci, params);
    magnetic_eci = magnetic_eci / norm(magnetic_eci);

    body_basis = make_basis(sun_body, magnetic_body);
    reference_basis = make_basis(sun_eci, magnetic_eci);
    q_measurement = dcm_to_quaternion(body_basis * reference_basis');
end

function basis = make_basis(first_vector, second_vector)
    first_vector = first_vector / norm(first_vector);
    second_axis = cross(first_vector, second_vector);
    second_axis = second_axis / norm(second_axis);
    third_axis = cross(first_vector, second_axis);
    basis = [first_vector, second_axis, third_axis];
end

function body_vector = calibrate_sun_measurement(measurement_body, params)
    if any(~isfinite(measurement_body))
        body_vector = nan(3, 1);
        return;
    end

    sensor_vector = params.sunsensor.mounting_dcm * measurement_body;
    body_vector = params.sunsensor.mounting_dcm' * ...
        params.sunsensor.misalignment' * sensor_vector;
    body_vector = body_vector / norm(body_vector);
end

function body_vector = calibrate_magnetic_measurement(measurement_body, params)
    sensor_vector = params.magnetometer.mounting_dcm * measurement_body;
    sensor_vector = sensor_vector - params.magnetometer.bias_vector;
    body_vector = params.magnetometer.mounting_dcm' * ...
        params.magnetometer.misalignment' * sensor_vector;
    body_vector = body_vector / norm(body_vector);
end

function q_product = quaternion_multiply(q_left, q_right)
    product_scalar = q_left(1) * q_right(1) - ...
        dot(q_left(2:4), q_right(2:4));
    product_vector = q_left(1) * q_right(2:4) + ...
        q_right(1) * q_left(2:4) + ...
        cross(q_left(2:4), q_right(2:4));
    q_product = [product_scalar; product_vector(:)];
end

function q_conjugate = quaternion_conjugate(q)
    q_conjugate = [q(1); -q(2:4)];
end
