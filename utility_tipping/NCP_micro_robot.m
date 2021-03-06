function A = NCP_micro_robot(A)
unit = A.unit;
unit_mass = A.unit_mass;
N = A.N;
index = A.n_i;

global e_t e_o e_r mu mu1 mu2;

e_t = A.ellipsoid(1);
e_o = A.ellipsoid(2);
e_r = A.ellipsoid(3)*unit;
mu =A.cof;
mu1 = mu;
mu2 = mu;

global M_r B F_elect V_m Van;
M_r = A.M_r*unit_mass/unit;
F_elect = A.F_elect*unit*unit_mass;
V_m = A.V_m*unit^3;
B = A.B;
Van = 0;
A.VAN(1) =Van;

global m g;
m = A.mass*unit_mass;
g = A.gravity*unit;

global len wid len1 len2 wid1 wid2 heg  I_xx I_yy I_zz Height phi theta;
if (A.shape == 'cuboid_shape')
    len = A.dim(1)*unit;  % in fixed body frame's y direction
    wid = A.dim(2)*unit;  % in fixed body frame's x direction
    heg = A.dim(3)*unit;  % in fixed body frame's z direction

    
    Height = A.dim(3)*unit/2;

elseif A.shape == 'spiked_shape'
    len1 = A.dim(1)*unit;
    len2 = A.dim(2)*unit;
    wid1 = A.dim(3)*unit;
    wid2 = A.dim(4)*unit;
    heg = A.dim(5)*unit;
    

elseif A.shape == 'spiked_ended'
    len1 = A.dim(1)*unit;
    len2 = A.dim(2)*unit;
    wid1 = A.dim(3)*unit;
    wid2 = A.dim(4)*unit;
    heg = A.dim(5)*unit;
    
        
end

I_xx = A.I_xx*unit^2*unit_mass;
I_yy = A.I_yy*unit^2*unit_mass;
I_zz = A.I_zz*unit^2*unit_mass;

    
phi = A.phi;

theta = A.theta;

global h;
h = A.h;

% applied wrenches
P_x = A.Impulses(:,1);
P_y = A.Impulses(:,2);
P_z = A.Impulses(:,3);
P_xt = A.Impulses(:,4);
P_yt = A.Impulses(:,5);
P_zt = A.Impulses(:,6);

global p_x p_y p_z p_xt p_yt p_zt;

global fqn time;
fqn = A.fqn;
%% initial configuration and state
% q_old - position and orientation vector at l, q_old=[q_xo;q_yo;q_zo;q0_o;q1_o;q2_o;q3_o]
global q_old;

q_old(1:3,1) = A.initial_q(1:3)*unit;
q_old(4:7,1) = A.initial_q(4:7);
% nu_old - generalized velocity vector at l, nu_old=[v_xo;v_yo;v_zo;w_xo;w_yo;w_zo]
global nu_old;

nu_old(1:3,1) = A.initial_v(1:3)*unit;
nu_old(4:6,1) = A.initial_v(4:6); 


%% define the infinity and initial guess and fun to use
 A = initial_guess(A);
 l = A.l;
 u = A.u;
 Z = A.Z;
 
 fun = A.fun;
 Q = q_old;
for i=1:N
    p_x = P_x(i);
    p_y = P_y(i);
    p_z = P_z(i);
    p_xt = P_xt(i);
    p_yt = P_yt(i);
    p_zt = P_zt(i);
    
    time = h*i;
    tic;
    [A.z(:,i,index),~,~,~,status] = pathmcp(Z,l,u,fun);
    time_NCP = toc;
    
     %% first solve for the solution
    if status == 1 
        A.time_NCP(i) = time_NCP;
    else
        A.time_NCP(i) = 0;
    end
    j = 1;


    if i == 1
        while status == 0
        j = j+1;
        Z_new = change_initial_guess(A,Z,Q);
        tic
        [A.z(:,i,index),~,~,~,status] = pathmcp(Z_new,l,u,fun);
        time_NCP = toc;
        if status == 1 
            A.time_NCP(i) = time_NCP;
        else
            A.time_NCP(i) = 0;
        end
            if j>=60
                error('Path can not found the solution, change your initial guess');
            end
        end

    else
        while status == 0
            j = j+1;
            Z_new = change_guess(A,Z,Q);
            tic
            [A.z(:,i,index),~,~,~,status] = pathmcp(Z_new,l,u,fun);
            time_NCP = toc;
            if status == 1 
                A.time_NCP(i) = time_NCP;
            else
                A.time_NCP(i) = 0;
            end
            if j>=60
                error('Path can not found the solution, change your initial guess');
            end
        end
    end
    
%     %% second check whether the constraints is satisfied or not
%     F_evaluation = A.check(A.z(:,i,index),0);
%     F_evaluation = F_evaluation(l==0);
%     while numel(F_evaluation(F_evaluation < tol))>0
%         status =0;
%         resolve(A,Z,Q,status,time_NCP,i);
%         F_evaluation = A.check(A.z(:,i,index),0);
%         F_evaluation = F_evaluation(l==0);
%     end
   
   
    
   [Q,Nu] = kinematic_map(q_old,A.z(:,i,index),h); % the function which returns the state vectors
   
   Z = A.z(:,i,index); % updating the initial guess for each iteration
   A.q(:,i) = Q; 
   
   A.F_evaluation(:,i) = A.check(Z,0);
   q_old = Q; % updating the beginning value for next time step
   nu_old = Nu; 
   %% determine the adhensive force coarsively
   Van = adhensive_force(A,i,index);
   A.VAN(i) = Van/(A.unit*A.unit_mass*A.h);
   i
end

end