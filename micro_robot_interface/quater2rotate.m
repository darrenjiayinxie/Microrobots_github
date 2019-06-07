% input is the q = q0+q1*i + q2*j +q3*k quaternion 
function output = quater2rotate(input)
output = zeros(4,1);
q0 = input(1);
q1 = input(2);
q2 = input(3);
q3 = input(4);

theta = 2*acos(q0);
if (sin(theta/2) == 0)
    w = [1;0;0];
else
    w = (1/sin(theta/2))*[q1;q2;q3];
end
output(1) = theta;
output(2:4) = w;

end