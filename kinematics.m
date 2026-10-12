function [q, DCM] = kinematics(q_current,w, dt)

Omega_w = [ 0,         -w(1), -w(2), -w(3);
            w(1),  0,          w(3), -w(2);
            w(2), -w(3),  0,          w(1);
            w(3),  w(2), -w(1),  0       ];

q_dot = 0.5 * Omega_w * q_current;
q = q_current + q_dot * dt;
q = q / norm(q); % Normalizasyon (Şart)
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4); %q = [q0; q1; q2; q3]; şeklinde
DCM = [q0^2+q1^2-q2^2-q3^2, 2*(q1*q2 + q0*q3),   2*(q1*q3 - q0*q2);
            2*(q1*q2 - q0*q3),   q0^2-q1^2+q2^2-q3^2, 2*(q2*q3 + q0*q1);
            2*(q1*q3 + q0*q2),   2*(q2*q3 - q0*q1),   q0^2-q1^2-q2^2+q3^2];
end