function movie_microrobot_spiked_2D_uneven_terrain(A)
N = A.N;

unit = A.unit;

len1 = A.dim(1)*unit;
len2 = A.dim(2)*unit;
wid1 = A.dim(3)*unit;
wid2 = A.dim(4)*unit;
heg = A.dim(5)*unit;

A.corner =0;
Q_x = A.q(1,end);
Q_y = A.q(2,end);
Theta = A.theta;
fqn = A.fqn;
h = A.h;
% ECP

% size of the window
ratio = 1;
x1 = -1e-3*unit+Q_x;
x2 = 1e-3*unit;
y1 = -1e-3*unit;
y2 = 1e-3*unit+ Q_y;

if abs(x1- x2)<1e-3*unit
    x1 = -6e-3*unit;
    x2 = 6e-3*unit;  
elseif x1>x2
    x1 = -6e-3*unit;
    x2 = 6e-3*unit;  
end


Size = [x1 x2 y1 y2];
% vertex of robot in body frame
p(1,:) = [len1/2 wid2/2 1];
p(2,:) = [len1/2 -wid2/2 1];
p(3,:) = [len2/2 -wid2/2 1];
p(4,:) = [0 -wid1/2 1];
p(5,:) = [-len2/2 -wid2/2 1];
p(6,:) = [-len1/2 -wid2/2 1];
p(7,:) = [-len1/2 wid2/2 1];
p(8,:) = [-len2/2 wid2/2 1];
p(9,:) = [0 wid1/2 1];
p(10,:) = [len2/2 wid2/2 1];

figure (1);
    set(gca,'nextplot','replacechildren'); 
    v = VideoWriter('test.avi');
    open(v);
    for i = 1:5:N
    
    q_x = A.q(1,i);
    q_y = A.q(2,i);
    theta_z = A.q(3,i);
    
   
    R = [cos(theta_z), -sin(theta_z);
         sin(theta_z), cos(theta_z)];
    H = [R,[q_x;q_y];zeros(1,2),1];
    
    % Microrobot 
    
    p_w = H*p';
  
    pgon = polyshape(p_w(1,:),p_w(2,:));
    plot(pgon);
    hold on
    plot(x1:x2-x1:x2,0*(x1:x2-x1:x2),'g');
    
    
    axis equal;
    axis(ratio*Size);
    
    
    for ii = 1:size(A.r,1)
        r1 = A.r(ii)*unit;
        xc(ii) = A.r_x(ii)*unit;
        yc = 0;

        theta = linspace(0,2*pi);
        cx1 = r1*cos(theta) + xc(ii);
        cy1 = r1*sin(theta) + yc;
    
    
    
        plot(cx1,cy1);
        hold on
    end
    
    hold off
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.01);
    
    end

   close(v);   

end