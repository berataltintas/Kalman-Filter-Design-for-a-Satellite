function [q_hat, bias_hat, P_hat] = eclipse_mekf(...
        current_gyro, current_magnetic_body, position_eci, ...
        q_prev, bias_prev, P_prev, dt, params)
    gyro_body = params.gyro.misalignment' * current_gyro;
    omega_estimated = gyro_body - bias_prev;
    [q_pred, ~] = kinematics(q_prev, omega_estimated, dt);

    omega_skew = [0, -omega_estimated(3), omega_estimated(2); ...
        omega_estimated(3), 0, -omega_estimated(1); ...
        -omega_estimated(2), omega_estimated(1), 0];
    F = [-omega_skew, -eye(3); zeros(3), zeros(3)];
    Phi = eye(6) + F * dt;

    Q = blkdiag(params.gyro.covariance * dt, ...
        params.gyro.bias_random_walk^2 * dt * eye(3));
    P_pred = Phi * P_prev * Phi' + Q;

    magnetic_body = calibrate_magnetic_measurement(...
        current_magnetic_body, params);
    magnetic_eci = get_magnetic_field(position_eci, params);
    magnetic_eci = magnetic_eci / norm(magnetic_eci);
    [~, C_pred] = kinematics(q_pred, zeros(3, 1), 0);
    magnetic_predicted_body = C_pred * magnetic_eci;

    innovation = magnetic_body - magnetic_predicted_body;
    H_attitude = finite_difference_jacobian(...
        q_pred, magnetic_eci, magnetic_predicted_body);
    H = [H_attitude, zeros(3)];
    R = params.coarse.magnetic_covariance;

    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    correction = K * innovation;

    delta_theta = correction(1:3);
    delta_q = [1; 0.5 * delta_theta];
    delta_q = delta_q / norm(delta_q);
    q_hat = quaternion_multiply(delta_q, q_pred);
    q_hat = q_hat / norm(q_hat);

    bias_hat = bias_prev + correction(4:6);
    identity = eye(6);
    P_hat = (identity - K * H) * P_pred * (identity - K * H)' + K * R * K';
    P_hat = 0.5 * (P_hat + P_hat');
end

function body_vector = calibrate_magnetic_measurement(measurement_body, params)
    sensor_vector = params.magnetometer.mounting_dcm * measurement_body;
    sensor_vector = sensor_vector - params.magnetometer.bias_vector;
    body_vector = params.magnetometer.mounting_dcm' * ...
        params.magnetometer.misalignment' * sensor_vector;
    body_vector = body_vector / norm(body_vector);
end

function H = finite_difference_jacobian(q, reference_eci, nominal_body, step)
    if nargin < 4
        step = 1e-7;
    end

    H = zeros(3, 3);
    for axis = 1:3
        delta_theta = zeros(3, 1);
        delta_theta(axis) = step;
        delta_q = [1; 0.5 * delta_theta];
        delta_q = delta_q / norm(delta_q);
        q_perturbed = quaternion_multiply(delta_q, q);
        [~, C_perturbed] = kinematics(q_perturbed, zeros(3, 1), 0);
        perturbed_body = C_perturbed * reference_eci;
        H(:, axis) = (perturbed_body - nominal_body) / step;
    end
end

function q_product = quaternion_multiply(q_left, q_right)
    product_scalar = q_left(1) * q_right(1) - ...
        dot(q_left(2:4), q_right(2:4));
    product_vector = q_left(1) * q_right(2:4) + ...
        q_right(1) * q_left(2:4) + ...
        cross(q_left(2:4), q_right(2:4));
    q_product = [product_scalar; product_vector(:)];
end
