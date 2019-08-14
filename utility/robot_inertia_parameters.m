function A = robot_inertia_parameters(A)


%% approximate dimension of the object 

if A.shape == 'cuboid_shape'
    A.dim=[800e-6 400e-6 100e-6]; %(m) length width height
    A.density = 2.1688e3; %kg/m^3
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
    A.dim=[800e-6 200e-6 650e-6 150e-6 300e-6]; %(m) len1 len2 wid1 wid2 heg
    A.density = 2.1688e3; %kg/m^3
    wid1 = A.dim(3);
    wid2 = A.dim(4);
    len1 = A.dim(1);
    len2 = A.dim(2);
    heg = A.dim(5);
    
    A.V_m = wid2*len1*heg+(wid1-wid2)*len2*heg/2; % m^3
    A.mass = A.V_m*A.density; %(kg)
    m = A.mass;
    
    tri_b =len2;
    tri_h = (wid1-wid2)/2;
    tri_a = tri_b/2;
    
    d_a = (tri_h/3+wid2/2);
    m1 = A.density*len1*wid2*heg;
    m_tri = 0.5*A.density*tri_b*tri_h*heg;
    I_xx_tri = A.density*heg*(tri_b^3*tri_h-tri_b^2*tri_h*tri_a+tri_b*tri_h*tri_a^2+tri_b*tri_h^3)/36;
    A.I_xx = (m1/12)*(len1^2+wid2^2) + (I_xx_tri + m_tri*d_a^2)*2;
    A.I_yy = (m1/12)*(heg^2+wid2^2);
    A.I_zz = (m1/12)*(heg^2+len1^2);
    
elseif A.shape == 'spiked_ended'
    A.dim=[800e-6 100e-6 150e-6 125e-6 400e-6]; %(m) len1 len2 wid1 wid2 heg
    A.density = 2.1688e3; %kg/m^3

    wid1 = A.dim(3);
    wid2 = A.dim(4);
    len1 = A.dim(1);
    len2 = A.dim(2);
    heg = A.dim(5);
    A.V_m = wid1*len1*heg+2*wid2*len2*heg; % m^3
    
    A.mass = A.V_m*A.density; %(kg)
    m = A.mass;
    
    tri_b =len2;
    tri_h = wid2;
    tri_a = tri_b/2;
    
    d_a = sqrt((len1-len2/2)^2+(wid1/2+wid2/3)^2);
    m1 = A.density*len1*wid1*heg;
    m_tri = 0.5*A.density*tri_b*tri_h*heg;
    I_xx_tri = A.density*heg*(tri_b^3*tri_h-tri_b^2*tri_h*tri_a+tri_b*tri_h*tri_a^2+tri_b*tri_h^3)/36;

    
    A.I_xx = (m1/12)*(len1^2+wid1^2) + (I_xx_tri + m_tri*d_a^2)*4;
    A.I_yy = (m/12)*(heg^2+wid1^2);
    A.I_zz = (m/12)*(wid1^2+len1^2);
    
elseif A.shape == 'geckod_shape'
    A.dim=[800e-6 400e-6 270e-6 150e-6]; %(m) length width height
    A.density = 2.1688e3; %kg/m^3
    A.mass = 9.4343e-8; %(kg)
    A.V_m = 3.6e-11; % m^3
elseif A.shape == 'curved_shape'
    A.dim=[1250e-6 1100e-6 400e-6 (25/180)*pi]; %(m) r_1 r_2 wid theta
    r1 = A.dim(1);
    r2 = A.dim(2);
    wid = A.dim(3);
    alpha = A.dim(4);
    D = ((r1^3-r2^3)/(r1^2-r2^2))*(2*sin(alpha)/(3*alpha));
    A.density = 2.1688e3; %kg/m^3
%     A.V_m = (pi*r1^2*wid-pi*r2^2*wid)*(alpha*180/pi)/180; % m^3
%     A.mass = A.V_m*A.density; %(kg)
    A.mass = 6.94e-8; %(kg)
    A.V_m = 3.2e-11; % m^3
end




end