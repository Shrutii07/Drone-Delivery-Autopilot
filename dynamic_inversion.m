close all
clear all 
clc 
% physical constants
% global g m L Ixx Iyy Izz k b         
g = 9.81; 
M = 0.5; 
L = 0.25;
Ixx = 5*10^-3;
Iyy = 5*10^-3;
Izz = 10*10^-3; 
k = 3*10^-6; 
b = 10^-7;
I = [Ixx, 0, 0;
     0, Iyy, 0;
     0, 0, Izz];

%constants
Kv = [20.5 ;20.5; 25];  %for positions
Kp = [155 ; 155; 300];
Kpe = [6; 5.2; 7];      %for angles
Kp2 = 4;
Kq2 = 3;
Kr2 = 4.9;

%Desired position and attitude
ref = [3 4 5]';
ref_dot = [0 0 0]';
ref_ddot = [0 0 0]';
Fstar = zeros(1000); %desired thrust
euler_star = [0 0 10]';  %desired euler angle
euler_dstar = [0 0 0]';
tau_star = zeros(3,1000);  %desired torque
w_star = [0 0 0]';

%Body Frame
w = zeros(3,1000);          % angular velocity in body frame
w_r = zeros(4,1000);        % individual rotor speed
% Inertial frame
pos = zeros(3,1000);        
pos(:,1)=[0 0 2];   % initial location
posdot = zeros(3,1000);
posddot = zeros(3,1000);
euler = pi./180*[-10 10 5]';
euler_dot  = [0 0 0]';

%Simulation

%R = [cos(euler(3,k))*cos(euler(2,k)) cos(euler(3,k))*sin(euler(2,k))*sin(euler(1,k))-sin(euler(3,k))*cos(euler(1,k)) cos(euler(3,k))*sin(euler(2,k))*cos(euler(1,k))+sin(euler(3,k))*sin(euler(1,k));
%     sin(euler(3,k))*cos(euler(2,k)) sin(euler(3,k))*sin(euler(2,k))*sin(euler(1,k))+cos(euler(3,k))*cos(euler(1,k)) sin(euler(3,k))*sin(euler(2,k))*cos(euler(1,k))-cos(euler(3,k))*sin(euler(1,k));
%     -sin(euler(2,k))                cos(euler(2,k))*sin(euler(1,k))                                                 cos(euler(2,k))*cos(euler(1,k))]
dt = 0.001;
for k = 1:1000
    Retw = [1 0 -sin(euler(2,k));
            0 cos(euler(1,k)) cos(euler(2,k))*sin(euler(1,k)); 
            0 -sin(euler(1,k)) cos(euler(2,k))*cos(euler(1,k))];
   
    posddot(:,k) = ref_ddot + Kv.*(ref_dot - posdot(:,k)) + Kp.*(ref - pos(:,k)); %commanded acc.
    Fstar(k) = (M*g + M*posddot(3,k))/(cos(euler(2,k))*cos(euler(1,k)));  % U1
    %desired eular angles
    euler_star(1,k+1) = asin((M*posddot(1,k)*sin(euler_star(3,k)) - M*posddot(2,k)*cos(euler_star(3,k)))/Fstar(k));
    euler_star(2,k+1) = asin((M*posddot(1,k)*cos(euler_star(3,k)) + M*posddot(2,k)*sin(euler_star(3,k)))/Fstar(k)*cos(euler_star(1,k)));
    euler_star(3,k+1) = euler_star(3,k);
    
    posdot(:,k+1) =  posdot(:,k) + dt*posddot(:,k);
    pos(:,k+1) = pos(:,k) + dt*posdot(:,k);
    
    euler_dot(:,k) = euler_dstar + Kpe.*(euler_star(:,k) - euler(:,k));
    euler(:,k+1) = euler(:,k) + dt*euler_dot(:,k);
    
    w(:,k) = Retw * euler_dot(:,k);
    
    tau_star(:,k) = [Ixx*Kp2*(w_star(1)-w(1,k)) + (Izz-Iyy)*w(2,k)*w(3,k); 
                     Iyy*Kq2*(w_star(2)-w(2,k)) + (Ixx-Izz)*w(1,k)*w(3,k); 
                     Izz*Kr2*(w_star(3)-w(3,k)) + (Iyy-Ixx)*w(1,k)*w(2,k)];
    w_r(:,k) = sqrt(inv([k k k k;
                         L*k 0 -L*k 0;
                         0 L*k 0 -L*k;
                         b -b b -b] ) * [Fstar(k);tau_star(:,k)]);
end
t = 1:0.1:101;
figure(1);
plot(t,euler_star(3,1)*ones(1001,1),t,euler(3,:),'b--',t,euler_star(2,1)*ones(1001,1),t,euler(2,:),'r:',t,euler_star(1,1)*ones(1001,1),t,euler(1,:),'g--');
figure(2);
plot(t,ref(1)*ones(1001,1),t,ref(2)*ones(1001,1),t,ref(3)*ones(1001,1),t,pos(1,:),'b--',t,pos(2,:),'r:',t,pos(3,:),'g--');

Z = stepinfo(pos(3,:),t,ref(3))
Y = stepinfo(pos(2,:),t,ref(2))
X = stepinfo(pos(1,:),t, ref(1))
