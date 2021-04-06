close all
clear all 
clc 

%initial positons velocities
x = [0;2;0];           % x, y, z fixed frame
euler = [10;0;30];   % pitch,roll,yaw fixed frame
v = [0;0;0];           % trans vel body frame
omega = [0;0;0];       % rot vel body frame
F1 = 1000;
F2 = 500;
F3 = 500;
F4 = 500;

%desired
x_des = [1;2;3];  %fixed frame
x_ddes = [0;0;0];  
euler_des = [0;0;60];
euler_ddes = [0;0;0];   
x1_des = [x_des(1);x_des(2)];
x1_ddes = [ x_ddes(1);x_ddes(2)];
x5_des = [euler_des(3);x_des(3)];
x5_ddes = [euler_ddes(3);x_ddes(3)];

%constant
Kt = diag([10^-2 10^-2 10^-2]);
Kr = diag([10^-3 10^-3 10^-3]);
m = 2;
d = 0.2;
c = 0.01;
g = 9.81;
Ix = 1.2416; Iy = 1.2416;
Iz = 2*1.2416;
It = diag([Ix Iy Iz]);
A1 = [3,0;0,3];
A2 = [3,0;0,3];
A3 = [3,0;0,3];
A4 = [3,0;0,3];
A5 = [3,0;0,3];
A6 = [3,0;0,3];
A7 = diag([3 3 3 3]);
zero = [0,0;0,0];
old_v = [0,0,0,0,0,0];

v1_dot = [0;0];
v2_dot = [0;0];
v3_dot = [0;0];
v4_dot = [0;0];
v5_dot = [0;0];
v6_dot = [0;0];
dt = 0.01;

for k =1:2000
    if k==1
        Rr = trans_mat_eul(euler(:,k));
        Rt = trans_mat_vel(euler(:,k)); 
    end

    eul_dot_ff = inv(Rr)*omega(:,k);
    x_dot_ff = Rt*v(:,k);
    
    x1(:,k) = [x(1,k);x(2,k)];
    x2(:,k) = [x_dot_ff(1);x_dot_ff(2)];
    x3(:,k) = [euler(1,k);euler(2,k)];
    x4(:,k) = [eul_dot_ff(1);eul_dot_ff(2)];
    x5(:,k) = [euler(3,k);x(3,k)];
    x6(:,k) = [eul_dot_ff(3);x_dot_ff(3)];

    J_phi = [0, 0, 0;
             0, -sin(x3(1,k)), cos(x3(1,k))*cos(x3(2,k));
             0, -cos(x3(1,k)), -sin(x3(1,k))*cos(x3(2,k))];

    J_theta = [ 0, 0, -cos(x3(2,k));
                0, 0, -sin(x3(2,k))*sin(x3(1,k));
                0, 0, -cos(x3(1,k))*sin(x3(2,k))];

    
    S = [0, -omega(3,k), omega(2,k);
         omega(3,k), 0 , -omega(1,k);
         -omega(2,k), omega(1,k), 0];
    
    f = -(Rt*Kt*Rt.'*x_dot_ff)./m - [0 0 g]';          % G matrix
    f_e = -(It*Rr)\(It*(J_phi*x4(1,k) +J_theta*x4(2,k))*eul_dot_ff - Kr*Rr*eul_dot_ff - cross((Rr*eul_dot_ff),(It*Rr*eul_dot_ff)))+[c*cos(euler(1,k))*tan(euler(2,k))*(F1(k)-F2(k)+F3(k)-F4(k))/Iz;
                                                                                                  -c*sin(euler(1,k))*(F1(k)-F2(k)+F3(k)-F4(k))/Iz;
                                                                                                  d*sin(euler(1,k))*sec(euler(2,k))*(F3(k)-F1(k))/Iy];
f
f_e
f0 = [f(1);f(2)];
    phi_0 = [sin(x3(1,k)) ;
             cos(x3(1,k))*sin(x3(2,k))];
    j0 = [cos(x3(1,k)), 0;
           -sin(x3(1,k))*sin(x3(2,k)), cos(x3(1,k))*cos(x3(2,k))];
    g0 =((F1(k) + F2(k) + F3(k) + F4(k))/m).*[sin(x5(1,k)), cos(x5(1,k));
                                              -cos(x5(1,k)), sin(x5(1,k))];

    f1 = [f_e(1);f_e(2)];
    g1 = [1/Ix, (sin(x3(1,k))*tan(x3(2,k)))/Iy;
          0    , cos(x3(1,k))/Iy];
    phi_1 = [d*(F2(k) - F4(k)); d*(F3(k)- F1(k))];                 
    j1 = [0, d,0,-d;
          -d,0,d, 0];

    f2 = [f_e(3);f(3)];
    g2 = [(cos(x3(1,k))*sec(x3(2,k)))/Iz, 0;
           0, (cos(x3(1,k))*cos(x3(2,k)))/m];
    phi_2 = [c*(F1(k) - F2(k) + F3(k) - F4(k)); (F1(k) + F2(k) + F3(k) + F4(k))];
    j2 = [c,-c,c,-c;
          1, 1, 1, 1];
phi_0
j0
g0
g1
phi_1
j1
g2
phi_2
j2
    v1 = A1*(x1_des-x1(:,k)) + x1_ddes;

    v2 = (g0)\((x1_des - x1(:,k)) + A2*(v1 - x2(:,k)) + v1_dot -f0);

    v3 = (j0)\(g0.'*(v1-x2(:,k)) +A3* (v2 - phi_0) +v2_dot);

    v4 = (g1)\(j0.'*(v2 - phi_0) + A4*(v3 - x4(:,k)) +v3_dot -f1);           %%% change x4 to x4(:,k)

    v5 = A5*(x5_des - x5(:,k)) + x5_ddes;

    v6 = (g2)\((x5_des - x5(:,k)) + A6*(v5 - x6(:,k)) +v5_dot -f2);

    u = ([j1;j2])\([g1, zero; zero, g2].' * [v3-x4(:,k);v5-x6(:,k)] + [v4_dot;v6_dot]+ A7*[v4 - phi_1;v6 - phi_2]);
    
    v1_dot = (v1 - old_v(1))/dt;
    v2_dot = (v2 - old_v(2))/dt;
    v3_dot = (v3 - old_v(3))/dt;
    v4_dot = (v4 - old_v(4))/dt;
    v5_dot = (v5 - old_v(5))/dt;
    v6_dot = (v6 - old_v(6))/dt;
    
    old_v=[v1,v2,v3,v4,v5,v6];
    
    F1(k+1) = u(1)*dt + F1(k);
    F2(k+1) = u(2)*dt + F2(k);
    F3(k+1) = u(3)*dt + F3(k);
    F4(k+1) = u(4)*dt + F4(k);
    
    x2_dot = f0 + g0*phi_0;
    x2(:,k+1) = x2_dot*dt + x2(:,k);
    x1(:,k+1) = x2(:,k)*dt + x1(:,k);
    
    x4_dot = f1 + g1*phi_1;
    x4(:,k+1) = x4_dot*dt + x4(:,k);
    x3(:,k+1) = x4(:,k)*dt + x3(:,k);
    
    x6_dot = f2 + g2*phi_2;
    x6(:,k+1) = x6_dot*dt + x6(:,k);
    x5(:,k+1) = x6(:,k)*dt + x5(:,k);
    
    x(:,k+1) = [x1(:,k+1);x5(2,k+1)];
    euler(:,k+1)= [x3(:,k+1);x5(1,k+1)];
    
    Rr = trans_mat_eul([x3(:,k+1);x5(1,k+1)]);
    Rt = trans_mat_vel([x3(:,k+1);x5(1,k+1)]);   
    
    omega(:,k+1) = Rr*[x4(:,k+1);x6(1,k+1)];
    v(:,k+1) = Rt*[x2(:,k+1);x5(2,k+1)];
end

t = 0:0.01:20;
%numel(t)
figure(1)
plot(t,x(1,:),t,x(2,:),t,x(3,:));
legend("x","y","z");
xlabel("time (s)");
ylabel("location (m)");

function R = trans_mat_eul(theta)
    phi = theta(1,1);
    th = theta(2,1);

%transformation matrix for body to inertial - angular quantities
  R = [ 1,       0,           -sin(th);
        0,  cos(phi), cos(th)*sin(phi);
        0, -sin(phi), cos(phi)*cos(th)];
    
end

function R = trans_mat_vel(theta)
    phi = theta(1,1);
    th = theta(2,1);
    psi = theta(3,1);
    
    R = [cos(phi)*cos(psi), cos(psi)*sin(th)*sin(phi)-sin(psi)*cos(phi), cos(phi)*sin(th)*cos(psi)+sin(phi)*sin(psi);
         sin(psi)*cos(th),  sin(phi)*sin(th)*sin(psi)+cos(psi)*cos(phi), sin(th)*sin(psi)*cos(phi)-cos(psi)*sin(phi);
         -sin(phi),         cos(th)*sin(phi)                           , cos(phi)*cos(th)];
    
end

function H = cross(a,b)
ax = a(1,1);
ay = a(2,1);
az = a(3,1);

R = [0, -az, ay;
     az, 0  -ax;
     -ay, ax, 0];
H = R*b;
end
