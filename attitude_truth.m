function attitude_body_eci = attitude_truth(position_eci, velocity_eci, sun_eci)
    x_body_eci = sun_eci / norm(sun_eci);

    y_reference = velocity_eci - dot(velocity_eci, x_body_eci) * x_body_eci;
    if norm(y_reference) < 1e-10
        y_reference = position_eci - dot(position_eci, x_body_eci) * x_body_eci;
    end
    y_body_eci = y_reference / norm(y_reference);
    z_body_eci = cross(x_body_eci, y_body_eci);
    z_body_eci = z_body_eci / norm(z_body_eci);
    y_body_eci = cross(z_body_eci, x_body_eci);

    attitude_body_eci = [x_body_eci'; y_body_eci'; z_body_eci'];
end