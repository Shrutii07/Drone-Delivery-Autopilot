close all
clear all
clc
% physical constants      
g = 9.81; 
M = 0.5; m = 0.01;
L = 0.25; LR = 0.2;
Ixx = 5*10^-3; Iyy = 5*10^-3; Izz = 10*10^-3; 
k = 3*10^-6; 
b = 10^-7;
I = diag([Ixx Iyy Izz]);

%constants
Kv = [0.63 ;0.63; 31.5];  %for positions
Kp = [0.1225 ; 0.1225; 6.125];
Kpe = [6; 5.2; 7];      %for angles
Kpw = [4; 3; 4.9];
Kvs = [94.5; 94.5];
Kps = [18.375; 18.375];

%Desired position and attitude
ref = [5 4 8]';
ref_dot = [0 0 0]';
ref_ddot = [0 0 0]';
Fstar = [0];
euler_star = pi./180*[0 0 -10]';
euler_dstar = [0 0 0]';
tau_star = [0 0 0]';
w_star = [0 0 0]';
w_dstar = [0 0 0]';
sway_star = [0 0]';
sway_dstar = [0 0]';


%Body Frame
w = [0 0 0]';
w_dot = [0 0 0]';
w_r = [0 0 0 0]';
% Inertial frame
pos = [0 0 3]';                   
posdot = [0 0 0]';
posddot = [0 0 0]';
euler = pi./180*[-10 5 10]';
euler_dot  = [0 0 0]';

sway = [0.2 0.15]';
sway_dot = [0 0]';
sway_ddot = [0 0]';

%cos(euler(1,1))
%Simulation
dt = 0.01;
for k = 1:3000
    
    posddot(:,k) = ref_ddot + Kv.*(ref_dot - posdot(:,k)) + Kp.*(ref - pos(:,k));
    posdot(:,k+1) = posdot(:,k) + dt*posddot(:,k);
    pos(:,k+1) = pos(:,k) + dt*posdot(:,k);
    
    sway_ddot(:,k) = Kvs.*(sway_dstar -sway_dot(:,k)) + Kps.*(sway_star - sway(:,k));
    sway_dot(:,k+1) = sway_dot(:,k) + dt*sway_ddot(:,k);
    sway(:,k+1) = sway(:,k) + dt*sway_dot(:,k);
    
    A = posddot(1,k) - LR*(-sin(sway(1,k))*sin(sway(2,k))*sway_ddot(2,k) - cos(sway(1,k))*sin(sway(2,k))*sway_dot(1,k)*sway_dot(2,k) - sin(sway(1,k))*cos(sway(2,k))*sway_dot(2,k)^2 + cos(sway(1,k))*cos(sway(2,k))*sway_ddot(1,k) - sin(sway(1,k))*cos(sway(2,k))*sway_dot(1,k)^2 - cos(sway(1,k))*sin(sway(2,k))*sway_dot(1,k)*sway_dot(2,k));
    B = posddot(2,k) + LR*(cos(sway(2,k))*sway_ddot(2,k) - sin(sway(2,k))*sway_dot(2,k)^2);
    C = posddot(3,k) - LR*(-cos(sway(1,k))*sin(sway(2,k))*sway_ddot(2,k) - cos(sway(1,k))*cos(sway(2,k))*sway_dot(2,k)*sway_dot(2,k) - sin(sway(1,k))*sin(sway(2,k))*sway_dot(2,k)*sway_dot(1,k) - sin(sway(1,k))*cos(sway(2,k))*sway_ddot(1,k) - cos(sway(1,k))*cos(sway(2,k))*sway_dot(1,k)^2 + sin(sway(1,k))*sin(sway(2,k))*sway_dot(1,k)*sway_dot(2,k));
    
    Fstar(k) = (M*g + M*posddot(3,k) + m*C)/(cos(euler(2,k))*cos(euler(1,k)));
    
    euler_star(1,k+1) = asin(((M*posddot(1,k) + m*A)*sin(euler_star(3,k)) - (M*posddot(2,k) + m*B)*cos(euler_star(3,k)))/Fstar(k));
    euler_star(2,k+1) = asin(((M*posddot(1,k) + m*A)*cos(euler_star(3,k)) + (M*posddot(2,k) + m*B)*sin(euler_star(3,k)))/Fstar(k)*cos(euler_star(1,k)));
    euler_star(3,k+1) = euler_star(3,k);
    
    euler_dot(:,k) = euler_dstar + Kpe.*(euler_star(:,k) - euler(:,k));
    euler(:,k+1) = euler(:,k) + dt*euler_dot(:,k);
    Retw = [1 0 -sin(euler(2,k));
            0 cos(euler(1,k)) cos(euler(2,k))*sin(euler(1,k)); 
            0 -sin(euler(1,k)) cos(euler(2,k))*cos(euler(1,k))];
    w(:,k) = Retw * euler_dot(:,k);
    
    w_dot(:,k) = w_dstar + Kpw.*[w_star(1) - w(1,k);
                                 w_star(2) - w(2,k);
                                 w_star(3) - w(3,k)];
    
    tau_star(:,k) = [Ixx*w_dot(1,k) + (Izz-Iyy)*w(2,k)*w(3,k); 
                     Iyy*w_dot(2,k) + (Ixx-Izz)*w(1,k)*w(3,k); 
                     Izz*w_dot(3,k) + (Iyy-Ixx)*w(1,k)*w(2,k)];
                 
    w_r(:,k) = sqrt(inv([k k k k;
                         L*k 0 -L*k 0;
                         0 L*k 0 -L*k;
                         b -b b -b] ) * [Fstar(k);tau_star(:,k)]);
end
t = 0:0.01:30;

figure(1);
plot(t,euler_star(3,1)*ones(k+1,1),t,euler(3,:),'b--',t,euler_star(2,1)*ones(k+1,1),t,euler(2,:),'r:',t,euler_star(1,1)*ones(k+1,1),t,euler(1,:),'g--'),legend('\psi^*','\psi','\theta^*','\theta','\phi^*','\phi');
figure(2);
plot(t,ref(1)*ones(k+1,1),t,ref(2)*ones(k+1,1),t,ref(3)*ones(k+1,1),t,pos(1,:),'b--',t,pos(2,:),'r--',t,pos(3,:),'g--'),legend('X^*','Y^*','Z^*','X','Y','Z');
figure(3);
plot(t, sway_star(1)*ones(k+1,1), t,sway(1,:),'b--', t, sway_star(2)*ones(k+1,1), t,sway(2,:),'r--'),legend('\alpha^*', '\alpha','\beta^*', '\beta');

Z = stepinfo(pos(3,:),t,ref(3));
Y = stepinfo(pos(2,:),t,ref(2));
X = stepinfo(pos(1,:),t, ref(1));