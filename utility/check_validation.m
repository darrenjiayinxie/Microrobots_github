function check_validation(A,Z)

if A.num_patches == 1
    
    A.l(1:15,1) = -infty;
    A.l(16:24,1) = 0;
    A.u(1:24,1) = infty;
     
    ECP = [q_x;q_y;a_z;q_x;q_y;a_z];
    
    v_t = v_x - w_z*(ECP(2) - q_y) + w_y*(ECP(3) - q_z);
    v_o = v_y + w_z*(ECP(1) - q_x) - w_x*(ECP(3) - q_z);
    v_r = w_z;
    
    sig = sqrt(e_t^2*v_t^2 + e_o^2*v_o^2 + e_r^2*v_r^2);
    p_n = m*g*h;
    if sig == 0
        p_t = 0;
        p_o = 0;
        p_r = 0;
    else
        p_t = -e_t^2*mu*p_n*v_t/sig;
        p_o = -e_o^2*mu*p_n*v_o/sig;
        p_r = -e_r^2*mu*p_n*v_r/sig;
    end
    
    
    
    
    Con_wrench = [p_t;p_o;p_r];
    
    La =[0;0;0;0;1;0;0]; %% assuming planar contact
    
    A.Z = [nu;ECP;Con_wrench;sig;La;p_n]; 
    
    A.fun = 'mcp_funjac_single_convex_contact_patch_cuboid';
    elseif A.num_patches == 1.5
    A.l(1:15,1) = -infty;
    A.l(16:26,1) = 0;
    A.u(1:26,1) = infty;
     
    ECP = [q_x;q_y;a_z;q_x;q_y;a_z];
    
    v_t = v_x - w_z*(ECP(2) - q_y) + w_y*(ECP(3) - q_z);
    v_o = v_y + w_z*(ECP(1) - q_x) - w_x*(ECP(3) - q_z);
    v_r = w_z;
    
    sig = sqrt(e_t^2*v_t^2 + e_o^2*v_o^2 + e_r^2*v_r^2);
    p_n = R*m*g*h;
    if sig == 0
        p_t = 0;
        p_o = 0;
        p_r = 0;
    else
        p_t = -e_t^2*mu*p_n*v_t/sig;
        p_o = -e_o^2*mu*p_n*v_o/sig;
        p_r = -e_r^2*mu*p_n*v_r/sig;
    end
    
    
    
    
    Con_wrench = [p_t;p_o;p_r];
    
    La =[1;0;0;0;0;0;0;0;0]; %% assuming planar contact
    
    A.Z = [nu;ECP;Con_wrench;sig;La;p_n]; 
    
    A.fun = 'mcp_funjac_convex_hull_T_bar';
elseif A.num_patches == 2
    A.l(1:24,1) = -infty;
    A.l(25:34,1) = 0;
    A.u(1:34,1) = infty;
    
    a11 = [q_x-r_b;q_y;a_z];
    a12 = a11;
    A1 = [a11;a12];

    a21 = [q_x+r_b;q_y;a_z];
    a22 = a21;
    A2 = [a21;a22];
    ECP = [A1;A2];
    
    v1_t = v_x - w_z*(A1(2) - q_y) + w_y*(A1(3) - q_z);
    v1_o = v_y + w_z*(A1(1) - q_x) - w_x*(A1(3) - q_z);
    v1_r = w_z;
    
    sig1 = sqrt(e_t^2*v1_t^2 + e_o^2*v1_o^2 + e_r^2*v1_r^2);
    

    p1_n = R*m*g*h; % from 0 to mgh
    p2_n = (1-R)*m*g*h;
    
    if sig1 == 0
        p1_t = 0;
        p1_o = 0;
        p1_r = 0;
    else
        p1_t = -e_t^2*mu*p1_n*v1_t/sig1;
        p1_o = -e_o^2*mu*p1_n*v1_o/sig1;
        p1_r = -e_r^2*mu*p1_n*v1_r/sig1;
    end
    
    
    
    v2_t = v_x - w_z*(A2(2) - q_y) + w_y*(A2(3) - q_z);
    v2_o = v_y + w_z*(A2(2) - q_x) - w_x*(A2(3) - q_z);
    v2_r = w_z;
    
    sig2 = sqrt(e_t^2*v2_t^2 + e_o^2*v2_o^2 + e_r^2*v2_r^2);
    
    if sig2 == 0
        p2_t = 0;
        p2_o = 0;
        p2_r = 0;
    else
        p2_t = -e_t^2*mu*p2_n*v2_t/sig2;
        p2_o = -e_o^2*mu*p2_n*v2_o/sig2;
        p2_r = -e_r^2*mu*p2_n*v2_r/sig2;
    end
    
    
    P1 = [p1_t;p1_o;p1_r];
    P2 = [p2_t;p2_o;p2_r];
    
    L1 = [1;0;0];
    L2 = [1;0;0];
    
    La = [L1;L2];
    Con_wrench = [P1; P2];
    sig = [sig1;sig2];
    P_n = [p1_n;p2_n];
    
    A.Z = [nu;ECP;Con_wrench;sig;La;P_n];
    
    A.fun = 'mcp_funjac_multiple_convex_contact_patches_2';
elseif A.num_patches == 3
    A.l(1:33,1) = -infty;
    A.l(34:48,1) = 0;
    A.u(1:48,1) = infty;
    
    Z = (1+R)*ones(48,1);
    
    if A.initial_guess==1
        A.initial_guess = Z;
    end
    A.fun = 'mcp_funjac_multiple_convex_contact_patches_3';
else 
    A.l(1:42,1) = -infty;
    A.l(43:62,1) = 0;
    A.u(1:62,1) = infty;
    
    Z = (1+R)*ones(62,1);
    
    if A.initial_guess==1
        A.initial_guess = Z;
    end
    A.fun = 'mcp_funjac_multiple_convex_contact_patches_4';
end



end