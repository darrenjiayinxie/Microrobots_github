function [r1_y,r2_y,rc_1,rc_2]=deter_pos_cylinders(A,Q)
r_y = A.r_y;
unit = A.unit;
N = size(A.r,1);
i = find((r_y)*unit<=Q(2), 1, 'first' );
j = find((r_y)*unit>=Q(2), 1, 'last');

if size(j,2) == 0
    r1_y = r_y(1)*unit;
    r2_y = r_y(2)*unit;
    rc_1 = A.r(1)*unit;
    rc_2 = A.r(2)*unit;
elseif size(i,2) == 0
        r1_y = r_y(N-1)*unit;
        r2_y = r_y(N)*unit;
        rc_1 = A.r(N-1)*unit;
        rc_2 = A.r(N)*unit;
else
    if i ==j
        r1_y = r_y(i)*unit;
        r2_y = r_y(i+1)*unit;
        rc_1 = A.r(i)*unit;
        rc_2 = A.r(i+1)*unit;
    else
        r1_y = r_y(j)*unit;
        r2_y = r_y(i)*unit;
        rc_1 = A.r(j)*unit;
        rc_2 = A.r(i)*unit;
    end
end


end