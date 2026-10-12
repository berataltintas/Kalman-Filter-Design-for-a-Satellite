function [position, velocity] = propagate_orbit(t, orbit_params)
    mu = 3.986004418e14; % Earth's gravitational parameter in m^3/s^2
    n = sqrt(mu / orbit_params.a^3); % Mean motion
    M = orbit_params.mean_anomaly + n * t; % Mean anomaly at time t

    %eccentric anomaly with newton - raphson
    E = M; % Initial guess for eccentric anomaly
    if orbit_params.e > 1e-8
        for iter = 1:10
            E = E - (E - orbit_params.e * sin(E) - M) / (1 - orbit_params.e * cos(E));
        end
    end

    %true anomaly
    nu = 2 * atan2(sqrt(1 + orbit_params.e) * sin(E / 2), sqrt(1 - orbit_params.e) * cos(E / 2));

    %perifocal
    p = orbit_params.a * (1 - orbit_params.e^2);
    r = p / (1 + orbit_params.e * cos(nu));

    position_perifocal = [r * cos(nu); r * sin(nu); 0];
    velocity_perifocal = sqrt(mu / p) * [-sin(nu); orbit_params.e + cos(nu); 0];

    %euler transformation
    R3 = [cos(orbit_params.raan), -sin(orbit_params.raan), 0;
          sin(orbit_params.raan), cos(orbit_params.raan), 0;
          0, 0, 1];

    R1 = [1, 0, 0;
          0, cos(orbit_params.i), -sin(orbit_params.i);
          0, sin(orbit_params.i), cos(orbit_params.i)];

    R3_arg = [cos(orbit_params.arg_perigee), -sin(orbit_params.arg_perigee), 0;
              sin(orbit_params.arg_perigee), cos(orbit_params.arg_perigee), 0;
              0, 0, 1];

    %transformation matrix
    R = R3 * R1 * R3_arg;
    position = R * position_perifocal;
    velocity = R * velocity_perifocal;
end