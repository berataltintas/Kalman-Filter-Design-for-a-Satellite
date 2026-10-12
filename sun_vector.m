function sun_vec = sun_vector(t, params)
    if nargin < 2
        sun_vec = [1; 0; 0];
        return;
    end

    epoch_direction = params.sun.epoch_eci_direction;
    epoch_direction = epoch_direction / norm(epoch_direction);
    epoch_angle = atan2(epoch_direction(2), epoch_direction(1));
    angle = epoch_angle + params.sun.annual_rate * t;
    sun_vec = [cos(angle); sin(angle); 0];
end