function [r_y,r_z] = generate_obstacles(A, N,r,delta_1,delta_2)
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