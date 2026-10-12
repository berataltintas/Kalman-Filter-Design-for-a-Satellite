function q=dcm_to_quaternion(C)
    %skaler part bulunuyor
    trace_c = trace(C);
    q0 = 0.5*sqrt(1+trace_c);
    %x y ve z bileşenleri bulunuyor
    q1 = (C(2,3)-C(3,2))/(4*q0);
    q2 = (C(3,1)-C(1,3))/(4*q0);
    q3 = (C(1,2)-C(2,1))/(4*q0);
    q=[q0;q1;q2;q3];
    %büyüklüğü bulunuyor
    q = q/norm(q);
end