function A = NCP_micro_robot_2D(A)

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

global len1 len2 wid1 wid2  heg  I_xx I_yy I_zz phi theta ;

len1 = A.dim(1)*unit;
len2 = A.dim(2)*unit;
wid1 = A.dim(3)*unit;
wid2 = A.dim(4)*unit;
heg = A.dim(5)*unit;

I_xx = A.I_xx*unit^2*unit_mass;
I_yy = A.I_yy*unit^2*unit_mass;
I_zz = A.I_zz*unit^2*unit_mass;

phi = A.phi;

theta = A.theta;

global h;
h = A.h;

global rc_1 r1_x r1_y;
rc_1 = A.r*unit;
r1_x = A.r1_x*unit;
r1_y = A.r1_y*unit;


% applied wrenches
P_x = A.Impulses(:,1);
P_y = A.Impulses(:,2);
P_zt = A.Impulses(:,3);

global p_x p_y p_zt;

global fqn time;
fqn = A.fqn;
%% initial configuration and state
% q_old - position and orientation vector at l
global q_old;

q_old(1:2,1) = A.initial_q(1:2)*unit;
q_old(3,1) = A.initial_q(3);
% nu_old - generalized velocity vector at l
global nu_old;

nu_old(1:2,1) = A.initial_v(1:2)*unit;
nu_old(3,1) = A.initial_v(3); 


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
          
        end
        
    else
        while status == 0
            j = j+1;
             Z_new = change_guess(A,Z);
              
           
            tic
            [A.z(:,i,index),f,J,Mu,status] = pathmcp(Z_new,l,u,fun);
            time_NCP = toc;
            if status == 1 
                A.time_NCP(i) = time_NCP;
            else
                A.time_NCP(i) = 0;
            end   
        end
    end
    %% determine the adhensive force coarsively
   
   % Van = adhensive_force(A,i,index);
    %A.VAN(i) = Van/(A.unit*A.unit_mass*A.h);
   Nu = A.z(1:3,i,index);
   Q = q_old +h*Nu;
   
   Z = A.z(:,i,index); % updating the initial guess for each iteration
   
   
   A.q(:,i) = Q; 
   
   q_old = Q; % updating the beginning value for next time step
   nu_old = Nu; 

   i
   
end
end