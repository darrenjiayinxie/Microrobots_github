function [r_x,r_y,r] = generate_obstacles(A, N)
delta_1 = 0.1e-3; %distance between robot and first bump

delta_2 = 1e-3; %distance between two neighboured bumps 
r1 = 0.36e-3;
r2 = 0.4e-3;
%r = r1:(r2-r1)/(N-1):r2;

for i = 1:N
    r(i) = r1;
end

for i = 1:N
    rc(i) = r(i);
    if i == 1
       r_x(i) = - (delta_1 +rc(i) +A.dim(1)/2 + A.dim(3))+A.initial_q(1);
    else
        r_x(i) = - (delta_2 + rc(i) +rc(i-1)) +r_x(i-1);
    end
    
    r_y(i) = 0;
    
end


end