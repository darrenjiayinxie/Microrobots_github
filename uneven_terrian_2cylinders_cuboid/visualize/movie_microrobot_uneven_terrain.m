function movie_microrobot(A)
N = A.N;

unit = A.unit;
if A.shape == 'cuboid_shape'
    L_r = A.dim(1)*unit;
    W_r = A.dim(2)*unit;
    H_r = A.dim(3)*unit;
elseif A.shape == 'spiked_shape'
    L_1 = A.dim(1)*unit;
    L_2 = A.dim(2)*unit;
    W_1 = A.dim(3)*unit;
    W_2 = A.dim(4)*unit;
    H_r = A.dim(5)*unit;
elseif A.shape == 'spiked_ended'
    L_1 = A.dim(1)*unit;
    L_2 = A.dim(2)*unit;
    W_1 = A.dim(3)*unit;
    W_2 = A.dim(4)*unit;
    H_r = A.dim(5)*unit;
elseif A.shape == 'geckod_shape'
    L_r = A.dim(1)*unit;
    W_r = A.dim(2)*unit;
    H_r = A.dim(4)*unit;
    H_L = A.dim(3)*unit;
elseif A.shape == 'curved_shape'
    r1 = A.dim(1)*unit;
    r2 = A.dim(2)*unit;
    wid = A.dim(3)*unit;
    alpha = A.dim(4);
    D = ((r1^3-r2^3)/(r1^2-r2^2))*(2*sin(alpha)/(3*alpha));
    
end

A.corner =0;
Q_y = A.q(2,end);
Q_z = A.q(3,end);
Theta = A.theta;
fqn = A.fqn;
h = A.h;
% ECP
a_x =A.z(7,:);
a_y =A.z(8,:);
a_z = A.z(9,:);

% size of the window
ratio = 1;
x1 = -1.5e-3*unit;
x2 = 1e-3*unit;
y1 = -1e-3*unit + Q_y;
y2 = 1e-3*unit;
z1 = -3e-3*unit;
z2 = 3e-3*unit+Q_z;
if abs(y1- y2)<1e-3*unit
    y1 = -6e-3*unit;
    y2 = 6e-3*unit;
    z1 = -6e-3*unit;
    z2 = 6e-3*unit;
elseif y1>y2
    y1 = -6e-3*unit;
    y2 = 6e-3*unit;
    z1 = -6e-3*unit;
    z2 = 6e-3*unit;    
end


size = [x1 x2 y1 y2 z1 z2];
% plane 

if A.corner == 0
    [X,Y] = meshgrid(x1/2:x2-x1/2:x2,y1:y2-y1:y2);
    Z =-tan(Theta)*Y;
    R_prime = [1,0,0;0,cos(2*pi-Theta),-sin(2*pi-Theta);0,sin(2*pi-Theta),cos(2*pi-Theta)];
else
    y1 = y_c;
    [X1,Y1] = meshgrid(x1/2:x2-x1/2:x2,y1:y2-y1:y2);
    Z1 = X1.*0+Y1.*0;
    
    y3 = -1e-3*unit + Q_y;
    [X2,Y2] = meshgrid(x1/2:x2-x1/2:x2,y3:y1-y3:y1);
    Z2 = -tan(Theta)*(Y2-y_c*(Y2./Y2));
    y1 = -1e-3*unit + Q_y;
    R_prime = eye(3);
end


figure (1);
    set(gca,'nextplot','replacechildren'); 
    v = VideoWriter('test.avi');
    open(v);
    for i = 1:5:N
    
    q_x = A.q(1,i);
    q_y = A.q(2,i);
    q_z = A.q(3,i);
    
    q0 = A.q(4,i);
    q1 = A.q(5,i);
    q2 = A.q(6,i);
    q3 = A.q(7,i);
    
    E = [-q1 q0 -q3 q2;
     -q2 q3 q0 -q1;
     -q3 -q2 q1 q0];

    G = [-q1 q0 q3 -q2;
         -q2 -q3 q0 q1;
         -q3 q2 -q1 q0];
    
    R = E*G';
    H = [R,[q_x;q_y;q_z];zeros(1,3),1];
    
    % Microrobot grid
    if A.shape == 'cuboid_shape'
        [Xr,Yr,Zr] = mashgrid_cuboid(L_r,W_r,H_r,H,[0;0;0]);
    elseif A.shape == 'spiked_shape'
        [Xr,Yr,Zr] = mashgrid_cuboid(L_1,H_r,W_2,H,[0;0;0]);
        v_1(1:4,1) = [H_r/2;0;W_1/2;1];
        v_1(1:4,2) = [H_r/2;-L_2/2;W_2/2;1];
        v_1(1:4,3) = [H_r/2; L_2/2;W_2/2;1];
        v_1(1:4,4) = [-H_r/2;0;W_1/2;1];
        v_1(1:4,5) = [-H_r/2;-L_2/2;W_2/2;1];
        v_1(1:4,6) = [-H_r/2; L_2/2;W_2/2;1];
        
        V_1 = H*v_1;
        V_1 = V_1(1:3,1:6);
        
        v_2(1:4,1) = [H_r/2;0;-W_1/2;1];
        v_2(1:4,2) = [H_r/2;-L_2/2;-W_2/2;1];
        v_2(1:4,3) = [H_r/2; L_2/2;-W_2/2;1];
        v_2(1:4,4) = [-H_r/2;0;-W_1/2;1];
        v_2(1:4,5) = [-H_r/2;-L_2/2;-W_2/2;1];
        v_2(1:4,6) = [-H_r/2; L_2/2;-W_2/2;1];
        
        V_2 = H*v_2;
        V_2 = V_2(1:3,1:6);
        Tri = [ 1 2 3; 4 5 6;1 2 4;4 5 2;1 3 4;4 6 3];
        
        trisurf(Tri,V_1(1,:),V_1(2,:),V_1(3,:),'FaceColor',[1 1 0]); 
        hold on
        trisurf(Tri,V_2(1,:),V_2(2,:),V_2(3,:),'FaceColor',[1 1 0]); 
        hold on
    elseif A.shape == 'curved_shape'
        [Xr,Yr,Zr] = mashgrid_curve(r1,r2,D,wid,alpha,H,[0;0;0]);
        
    elseif A.shape == 'spiked_ended'
        [Xr,Yr,Zr] = mashgrid_cuboid(L_1,H_r,W_1,H,[0;0;0]);
        v_1(1:4,1) = [H_r/2;(L_1-L_2)/2;W_1/2+W_2;1];
        v_1(1:4,2) = [H_r/2;L_1/2-L_2;W_1/2;1];
        v_1(1:4,3) = [H_r/2; L_1/2;W_1/2;1];
        v_1(1:4,4) = [-H_r/2;(L_1-L_2)/2;W_1/2+W_2;1];
        v_1(1:4,5) = [-H_r/2;L_1/2-L_2;W_1/2;1];
        v_1(1:4,6) = [-H_r/2; L_1/2;W_1/2;1];
        
        V_1 = H*v_1;
        V_1 = V_1(1:3,1:6);
        
        v_2(1:4,1) = [H_r/2;-(L_1-L_2)/2;W_1/2+W_2;1];
        v_2(1:4,2) = [H_r/2;-(L_1/2-L_2);W_1/2;1];
        v_2(1:4,3) = [H_r/2; -L_1/2;W_1/2;1];
        v_2(1:4,4) = [-H_r/2;-(L_1-L_2)/2;W_1/2+W_2;1];
        v_2(1:4,5) = [-H_r/2;-(L_1/2-L_2);W_1/2;1];
        v_2(1:4,6) = [-H_r/2; -L_1/2;W_1/2;1];
        
        V_2 = H*v_2;
        V_2 = V_2(1:3,1:6);
        
        v_3(1:4,1) = [H_r/2;-(L_1-L_2)/2;-(W_1/2+W_2);1];
        v_3(1:4,2) = [H_r/2;-(L_1/2-L_2);-W_1/2;1];
        v_3(1:4,3) = [H_r/2; -L_1/2;-W_1/2;1];
        v_3(1:4,4) = [-H_r/2;-(L_1-L_2)/2;-(W_1/2+W_2);1];
        v_3(1:4,5) = [-H_r/2;-(L_1/2-L_2);-W_1/2;1];
        v_3(1:4,6) = [-H_r/2; -L_1/2;-W_1/2;1];
        
        V_3 = H*v_3;
        V_3 = V_3(1:3,1:6);
        
        v_4(1:4,1) = [H_r/2;(L_1-L_2)/2;-(W_1/2+W_2);1];
        v_4(1:4,2) = [H_r/2;L_1/2-L_2;-W_1/2;1];
        v_4(1:4,3) = [H_r/2; L_1/2;-W_1/2;1];
        v_4(1:4,4) = [-H_r/2;(L_1-L_2)/2;-(W_1/2+W_2);1];
        v_4(1:4,5) = [-H_r/2;L_1/2-L_2;-W_1/2;1];
        v_4(1:4,6) = [-H_r/2; L_1/2;-W_1/2;1];
        
        V_4 = H*v_4;
        V_4 = V_4(1:3,1:6);
        
        Tri = [ 1 2 3; 4 5 6;1 2 4;4 5 2;1 3 4;4 6 3];
        
        trisurf(Tri,V_1(1,:),V_1(2,:),V_1(3,:),'FaceColor',[1 1 0]); 
        hold on
        trisurf(Tri,V_2(1,:),V_2(2,:),V_2(3,:),'FaceColor',[1 1 0]); 
        hold on
        trisurf(Tri,V_3(1,:),V_3(2,:),V_3(3,:),'FaceColor',[1 1 0]); 
        hold on
        trisurf(Tri,V_4(1,:),V_4(2,:),V_4(3,:),'FaceColor',[1 1 0]); 
        hold on
        
    elseif A.shape == 'geckod_shape'
        [Xr,Yr,Zr] = mashgrid_cuboid(L_r,W_r,H_r,H,[0;0;0]);
        N_w = 10;
        d_w = W_r/N_w;
        
        for i_w = 1:N_w
            D_w = -(W_r/2 -(i_w-1)*d_w - d_w/4);
            b = d_w/(H_L-H_r);
            C = -(D_w+0.5*H_L*b);
            [Y_v,Z_v] = meshgrid(-L_r/2:L_r:L_r/2,H_r/2:(H_L-H_r)/2:H_L/2);
            X_v = Y_v.*0 + Z_v.*0 +D_w ;
            
            [Y_s,Z_s] = meshgrid(-L_r/2:L_r:L_r/2,H_r/2:(H_L-H_r)/2:H_L/2);
            X_s = -b*Z_s-C;
            
            for i_1 = 1:2
                for j_1 = 1:2
                    A_m = H*[X_v(i_1,j_1);Y_v(i_1,j_1);Z_v(i_1,j_1);1];
                    X_v(i_1,j_1) = A_m(1);
                    Y_v(i_1,j_1) = A_m(2);
                    Z_v(i_1,j_1) = A_m(3);
                    A_m = H*[X_s(i_1,j_1);Y_s(i_1,j_1);Z_s(i_1,j_1);1];
                    X_s(i_1,j_1) = A_m(1);
                    Y_s(i_1,j_1) = A_m(2);
                    Z_s(i_1,j_1) = A_m(3);
                end
            end
            surf(X_v,Y_v,Z_v,'FaceColor',[1 1 0],'Linestyle','-');  
            hold on
            surf(X_s,Y_s,Z_s,'FaceColor',[1 1 0],'Linestyle','-');  
            hold on
        end
        
        for i_w = 1:N_w
            D_w = -(W_r/2 -i_w*d_w + d_w/4);
            b = d_w/(H_L-H_r);
            C =-(D_w-0.5*H_L*b);
            [Y_v,Z_v] = meshgrid(-L_r/2:L_r:L_r/2,-H_L/2:(H_L-H_r)/2:-H_r/2);
            X_v = Y_v.*0 + Z_v.*0 +D_w ;
            
            [Y_s,Z_s] = meshgrid(-L_r/2:L_r:L_r/2,-H_L/2:(H_L-H_r)/2:-H_r/2);
            X_s = -b*Z_s-C;
            
            for i_1 = 1:2
                for j_1 = 1:2
                    A_m = H*[X_v(i_1,j_1);Y_v(i_1,j_1);Z_v(i_1,j_1);1];
                    X_v(i_1,j_1) = A_m(1);
                    Y_v(i_1,j_1) = A_m(2);
                    Z_v(i_1,j_1) = A_m(3);
                    A_m = H*[X_s(i_1,j_1);Y_s(i_1,j_1);Z_s(i_1,j_1);1];
                    X_s(i_1,j_1) = A_m(1);
                    Y_s(i_1,j_1) = A_m(2);
                    Z_s(i_1,j_1) = A_m(3);
                end
            end
            surf(X_v,Y_v,Z_v,'FaceColor',[1 1 0],'Linestyle','-');  
            hold on
            surf(X_s,Y_s,Z_s,'FaceColor',[1 1 0],'Linestyle','-');  
            hold on
        end
        
    end

    
    
 

    % Mircorobot plot
    surf(Xr,Yr,Zr,'FaceColor',[1 1 0],'Linestyle','-');  
    
    
    axis equal;
    axis(ratio*size);
    hold on
    
    % ground
    if A.corner == 0
        surf(X,Y,Z,'FaceColor',[0 1 0],'Linestyle','-'); 
    else
        surf(X1,Y1,Z1,'FaceColor',[0 1 0],'Linestyle','-'); 
        hold on
        surf(X2,Y2,Z2,'FaceColor',[0 1 0],'Linestyle','-'); 
    end
     
    axis equal;
    axis(ratio*size);
    hold on
    
    % sphere
    if A.cylinder == 1
        [xs,ys,zs] = cylinder;
        Ra = [cos(pi/2),0,sin(pi/2);0,1,0;-sin(pi/2),0,cos(pi/2)];
        P1 = Ra*[xs(1,:);ys(1,:);zs(1,:)];
        P2 = Ra*[xs(2,:);ys(2,:);zs(2,:)];
        PP = [A.rc_1*[P1(1,:);P2(1,:)]*unit,A.rc_1*[P1(2,:);P2(2,:)]*unit+A.r1_y*unit,A.rc_1*[P1(3,:);P2(3,:)]*unit+A.r1_z*unit];
        surf(3*PP(:,1:21),PP(:,22:42),PP(:,43:63),'FaceColor',[0 1 0],'Linestyle','-');
        surf(-3*PP(:,1:21),PP(:,22:42),PP(:,43:63),'FaceColor',[0 1 0],'Linestyle','-');
        hold on
        
        P1 = Ra*[xs(1,:);ys(1,:);zs(1,:)];
        P2 = Ra*[xs(2,:);ys(2,:);zs(2,:)];
        PP = [A.rc_2*[P1(1,:);P2(1,:)]*unit,A.rc_2*[P1(2,:);P2(2,:)]*unit+A.r2_y*unit,A.rc_2*[P1(3,:);P2(3,:)]*unit+A.r2_z*unit];
        surf(3*PP(:,1:21),PP(:,22:42),PP(:,43:63),'FaceColor',[0 1 0],'Linestyle','-');
        surf(-3*PP(:,1:21),PP(:,22:42),PP(:,43:63),'FaceColor',[0 1 0],'Linestyle','-');
    end
    
    % frame
    xr = -1e-3*unit;
    yr = 0;
    L = 1e-3*unit;
    
    time = h*i;
    theta = -2*pi*fqn*time+A.theta;
    if y2-y1 <=5e-3*unit
        L = 0.5e-3*unit;
        d = 3e-3*unit;
        w = 3*abs((2*L)/(d));
    else
        w = 5*abs((2*L)/(y2-y1));
    end
    
    
    p1 = [xr,yr,0;xr,yr,0;xr,yr,0];
    p21 = [xr+L,yr,0]';
    p22 = [xr,yr+L,0]';
    p23 = [xr,yr,L]';
    p2 = [p21';p22';p23'];
    arrow3(p1,p2,'k1',w,3*w);
    axis equal;
    axis(ratio*size);
    hold on
    
    % direction of magnetic field
    p1 = [xr,yr,0];
    p2 = [xr,yr-cos(theta)*L,sin(theta)*L]';
    arrow3(p1,p2','r1',w,3*w);
    axis equal;
    axis(ratio*size);
    if unit == 1e6
        xlabel('x (\mu m)');
        ylabel('y (\mu m)');
        zlabel('z (\mu m)');
    elseif unit == 1e3
        xlabel('x (mm)');
        ylabel('y (mm)');
        zlabel('z (mm)');
    end
    hold on
    view([-1,-1,0.5]);
    hold off
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.01);
    
    end

   close(v);   

end