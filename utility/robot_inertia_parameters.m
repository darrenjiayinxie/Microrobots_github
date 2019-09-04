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
    length_r = A.dim(1);
    length_s = A.dim(2);
%     height_r = A.dim(3);
%     height_s = A.dim(4);
    height_r = A.dim(4);
    height_s = (A.dim(3)-A.dim(4))/2;
    width_r = A.dim(5);
    width_s = A.dim(5);
    density = A.density;
    
    A.I_xx = density*((height_r*length_r*width_r^3)/12 + (height_r*length_r^3*width_r)/12)...
            + 2*density*((height_s*length_s*width_s^3)/24 + (height_s*length_s^3*width_s)/48);
    A.I_yy = density*((height_r^3*length_r*width_r)/12 + (height_r*length_r*width_r^3)/12)...
            + 2*density*((height_s^3*length_s*width_s)/36 + (height_s*length_s*width_s^3)/24)...
            + 2*density*0.5*length_s*height_s*width_s * ((0.5*height_r + (1/3)*height_s)^2);
    A.I_zz = density*((height_r^3*length_r*width_r)/12 + (height_r*length_r^3*width_r)/12)...
            + 2*density*((height_s^3*length_s*width_s)/36 + (height_s*length_s^3*width_s)/48)...
            + 2*density*0.5*length_s*height_s*width_s * ((0.5*height_r + (1/3)*height_s)^2);
elseif A.shape == 'spiked_ended'
    

    wid1 = A.dim(3);
    wid2 = A.dim(4);
    len1 = A.dim(1);
    len2 = A.dim(2);
    heg = A.dim(5);
    A.V_m = wid1*len1*heg+2*wid2*len2*heg; % m^3
    
    A.mass = A.V_m*A.density; %(kg)
    
    %% Chenghao, you can put your code here
    length_r = A.dim(1);
    length_s = A.dim(2);
    height_r = A.dim(3);
    height_s = A.dim(4);
    width_r = A.dim(5);
    width_s = A.dim(5);
    density = A.density;
    
    A.I_xx = density*((height_r*length_r*width_r^3)/12 + (height_r*length_r^3*width_r)/12)...
            + 4*density*((height_s*length_s*width_s^3)/24 + (height_s*length_s^3*width_s)/48)...
            + 4*density*0.5*height_s*length_s*width_s * ((0.5*length_r - 0.5*length_s)^2);
    A.I_yy = density*((height_r^3*length_r*width_r)/12 + (height_r*length_r*width_r^3)/12)...
            + 4*density*((height_s^3*length_s*width_s)/36 + (height_s*length_s*width_s^3)/24)...
            + 4*density*0.5*height_s*length_s*width_s * ((0.5*height_r + (1/3)*height_s)^2);
    A.I_zz = density*((height_r^3*length_r*width_r)/12 + (height_r*length_r^3*width_r)/12)...
            + 4*density*((height_s^3*length_s*width_s)/36 + (height_s*length_s^3*width_s)/48)...
            + 4*density*0.5*height_s*length_s*width_s * ((0.5*height_r + (1/3)*height_s)^2 + (0.5*length_r - 0.5*length_s)^2);
end




end