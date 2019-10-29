function [r1_x,r2_x,rc_1,rc_2]=deter_pos_cylinders(A,Q)
r_x = A.r_x;
unit = A.unit;
N = size(A.r,1);
i = find((r_x)*unit<=Q(1), 1, 'first' );
j = find((r_x)*unit>=Q(1), 1, 'last');

if size(j,2) == 0
    r1_x = r_x(1)*unit;
    r2_x = r_x(2)*unit;
    rc_1 = A.r(1)*unit;
    rc_2 = A.r(2)*unit;
elseif size(i,2) == 0
        r1_x = r_x(N-1)*unit;
        r2_x = r_x(N)*unit;
        rc_1 = A.r(N-1)*unit;
        rc_2 = A.r(N)*unit;
else
    if i ==j
        r1_x = r_x(i)*unit;
        r2_x = r_x(i+1)*unit;
        rc_1 = A.r(i)*unit;
        rc_2 = A.r(i+1)*unit;
    else
        r1_x = r_x(j)*unit;
        r2_x = r_x(i)*unit;
        rc_1 = A.r(j)*unit;
        rc_2 = A.r(i)*unit;
    end
end


end