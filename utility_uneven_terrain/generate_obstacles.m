function [r_y,r_z,r] = generate_obstacles(A, N)
delta_1 = 0.05e-3; %distance between robot and first bump

delta_2 = 0.3e-3; %distance between two neighboured bumps 
r1 = 0.3e-3;
r2 = 0.4e-3;
r = r1:(r2-r1)/(N-1):r2;

for i = 1:N
    
    rc(i) = r(i);
    if i == 1
       r_y(i) = - (delta_1 +rc(i) +A.dim(1)/2 + A.dim(3))+A.initial_q(2);
    else
        r_y(i) = - (delta_2 + rc(i) +rc(i-1)) +r_y(i-1);
    end
    
    r_z(i) = 0;
    
end


end