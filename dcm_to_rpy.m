function rpy = dcm_to_rpy(attitude_body_eci)
    pitch = -asin(max(-1, min(1, attitude_body_eci(1, 3))));
    roll = atan2(attitude_body_eci(2, 3), attitude_body_eci(3, 3));
    yaw = atan2(attitude_body_eci(1, 2), attitude_body_eci(1, 1));
    rpy = [roll; pitch; yaw];
end