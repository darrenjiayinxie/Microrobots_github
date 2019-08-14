function A = robot_inertia_parameters(A)


%% approximate dimension of the object 

if A.shape == 'cuboid_shape'
    
    len = A.dim(1);
    wid = A.dim(2);
    heg = A.dim(3);
    
    
    A.V_m = len*heg*wid; % m^3
    
    A.mass = A.density*A.V_m; %(kg)
    m = A.mass;
    A.I_xx = (m/12)*(len^2+heg^2);
    A.I_yy = (m/12)*(heg^2+wid^2);
    A.I_zz = (m/12)*(wid^2+len^2);
elseif A.shape == 'spiked_shape'
    
    wid1 = A.dim(3);
    wid2 = A.dim(4);
    len1 = A.dim(1);
    len2 = A.dim(2);
    heg = A.dim(5);
    
    A.V_m = wid2*len1*heg+(wid1-wid2)*len2*heg/2; % m^3
    A.mass = A.V_m*A.density; %(kg)
    
    %% Chenghao, you can put your code here
    %I_xx = 
    %I_yy = 
    %I_zz = 
elseif A.shape == 'spiked_ended'
    

    wid1 = A.dim(3);
    wid2 = A.dim(4);
    len1 = A.dim(1);
    len2 = A.dim(2);
    heg = A.dim(5);
    A.V_m = wid1*len1*heg+2*wid2*len2*heg; % m^3
    
    A.mass = A.V_m*A.density; %(kg)
    
    %% Chenghao, you can put your code here
    %I_xx = 
    %I_yy = 
    %I_zz = 
end




end