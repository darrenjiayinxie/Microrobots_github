function A = NCP_micro_robot(A)
load('history.mat');
load('history_dis.mat');
size_history = size(history,2);
%size_history = 0;
hist = 0;

unit = A.unit;
unit_mass = A.unit_mass;
N = A.N;
index = A.n_i;
global e_t e_o e_r mu ;

e_t = A.ellipsoid(1);
e_o = A.ellipsoid(2);
e_r = A.ellipsoid(3)*unit;
mu =A.cof;

global M_r B F_elect V_m;
M_r = A.M_r*unit_mass/unit;
F_elect = A.F_elect*unit*unit_mass;
V_m = A.V_m*unit^3;
B = A.B;


global m g;
m = A.mass*unit_mass;
g = A.gravity*unit;

global len wid  heg  I_xx I_yy I_zz Height phi theta ;
if (A.shape == 'cuboid_shape')
    len = A.dim(1)*unit;  % in fixed body frame's x direction
    wid = A.dim(2)*unit;  % in fixed body frame's y direction
    heg = A.dim(3)*unit;  % in fixed body frame's z direction

    I_xx = (m/12)*(len^2+heg^2);
    I_yy = (m/12)*(heg^2+wid^2);
    I_zz = (m/12)*(wid^2+len^2);

    Height = A.dim(3)*unit/2;
end
phi = A.phi;

theta = A.theta;

global h;
h = A.h;

global rc_1 r1_y r1_z;
rc_1 = A.r(1)*unit;
r1_y = A.r_y(1)*unit;
r1_z = A.r_z(1)*unit;

global rc_2 r2_y r2_z;
rc_2 = A.r(2)*unit;
r2_y = A.r_y(2)*unit;
r2_z = A.r_z(2)*unit;


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
    [A.z(:,i,index),f,J,Mu,status] = pathmcp(Z,l,u,fun);
    time_NCP = toc;
    if status == 1 
        A.time_NCP(i) = time_NCP;
    else
        A.time_NCP(i) = 0;
    end
    j = 1;
    if i == 1
        while status == 0
        j = j+1;
        Z_new = change_initial_guess(A,Z);
        tic
        [A.z(:,i,index),f,J,Mu,status] = pathmcp(Z_new,l,u,fun);
        time_NCP = toc;
        if status == 1 
            A.time_NCP(i) = time_NCP;
        else
            A.time_NCP(i) = 0;
        end
            if j>=30
                history = A.z;
                history_dis = A.r2_y;
                save('history.mat','history');
                error('Path can not found the solution, change your initial guess');
            end
        end
        
    else
        while status == 0
            j = j+1;
            if hist == 1
                Z_new = A.z(:,i-1,index);
                hist = 0;
            else
                Z_new = change_guess(A,Z,Q,r1_y,r2_y,rc_1,rc_2);
            end
            tic
            [A.z(:,i,index),f,J,Mu,status] = pathmcp(Z_new,l,u,fun);
            time_NCP = toc;
            if status == 1 
                A.time_NCP(i) = time_NCP;
            else
                A.time_NCP(i) = 0;
            end
            if j>=10
                history = A.z;
                history_dis = A.r_y(2);
                save('history_dis.mat','history_dis');
                save('history.mat','history');
                error('Path can not found the solution, change your initial guess');
            end
        end
    end
    %% determine the adhensive force coarsively
   
   % Van = adhensive_force(A,i,index);
    %A.VAN(i) = Van/(A.unit*A.unit_mass*A.h);
   [Q,Nu] = kinematic_map(q_old,A.z(:,i,index),h); % the function which returns the state vectors
   
   if (i <size_history)&&(history_dis == A.r_y(2))
       Z = history(:,i+1);
       hist = 1;
   else
       Z = A.z(:,i,index); % updating the initial guess for each iteration
   end
   
   A.q(:,i) = Q; 
   
   A.F_evaluation(:,i) = A.check(Z,0);
   q_old = Q; % updating the beginning value for next time step
   nu_old = Nu; 
   [r1_y,r2_y,rc_1,rc_2]=deter_pos_cylinders(A,Q);
   
   i
   
end
history = A.z;
history_dis = A.r_y(2);
save('history_dis.mat','history_dis');
save('history.mat','history');
end