function q_initial = coarse_triad(...
    sun_measurement_body, magnetic_measurement_body, position_eci, ...
    sun_eci, params)
    sun_body = calibrate_sun(sun_measurement_body, params);
    magnetic_body = calibrate_magnetic(magnetic_measurement_body, params);

    sun_eci = sun_eci / norm(sun_eci);
    magnetic_eci = get_magnetic_field(position_eci, params);
    magnetic_eci = magnetic_eci / norm(magnetic_eci);

    body_basis = make_basis(sun_body, magnetic_body);
    reference_basis = make_basis(sun_eci, magnetic_eci);
    attitude_body_eci = body_basis * reference_basis';
    q_initial = dcm_to_quaternion(attitude_body_eci);
end

function basis = make_basis(first_vector, second_vector)
    first_vector = first_vector / norm(first_vector);
    second_axis = cross(first_vector, second_vector);
    second_axis = second_axis / norm(second_axis);
    third_axis = cross(first_vector, second_axis);
    basis = [first_vector, second_axis, third_axis];
end

function body_vector = calibrate_sun(measurement_body, params)
    sensor_vector = params.sunsensor.mounting_dcm * measurement_body;
    body_vector = params.sunsensor.mounting_dcm' * ...
        params.sunsensor.misalignment' * sensor_vector;
    body_vector = body_vector / norm(body_vector);
end

function body_vector = calibrate_magnetic(measurement_body, params)
    sensor_vector = params.magnetometer.mounting_dcm * measurement_body;
    sensor_vector = sensor_vector - params.magnetometer.bias_vector;
    body_vector = params.magnetometer.mounting_dcm' * ...
        params.magnetometer.misalignment' * sensor_vector;
    body_vector = body_vector / norm(body_vector);
end
