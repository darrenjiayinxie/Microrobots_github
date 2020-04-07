function A = planner(A)
N =1000;
A.N = N;

A.Impulses = zeros(N,6);

%% external impulses and angular momentums
for i = 1:N
    P_x = 0;
    P_y = 0;
    P_z = 0;
    P_xt = 0;
    P_yt = 0;
    P_zt = 0;
    
    A.Impulses(i,1) = P_x;
    A.Impulses(i,2) = P_y;
    A.Impulses(i,3) = P_z;
    A.Impulses(i,4) = P_xt;
    A.Impulses(i,5) = P_yt;
    A.Impulses(i,6) = P_zt;
    
    
    
end








end