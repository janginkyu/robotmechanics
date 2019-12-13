clear;
close all;

%% problem 1

% import robot from urdf file
robot1 = importrobot("robot1.urdf");
robot1.DataFormat = 'row';

% initial UGV position
UGVposition = [0.75; -0.5; 0.];

% desired position of UGV seen in ee frame
p0 = [3.0*sqrt(2); 0.0; 0.0; 1.0];

% maximum iteration for inverse kinematics calc
maxiter = 1000;

q = homeConfiguration(robot1);

%%
lambda = 0.04;
for it = 1:maxiter
    % current position of ray end
    g = getTransform(robot1, q, 'ee', 'world');
    p = g * p0;
    p0dir = p0;
    p0dir(4) = 0;
    p0_hat = hat(g * p0dir);
    
    % classical jacobian of ee w.r.t. world
    jac = geometricJacobian(robot1, q, 'ee');
    
    Jv = jac(4:6, :);
    Jw = jac(1:3, :);
    
    Jp = Jv - p0_hat * Jw;
    q = q - lambda * (p(1:3, 1) - UGVposition)' * Jp;
    
    if norm(p(1:3, 1) - UGVposition) < 1e-5
        break;
    end
    
%     config(1).JointPosition = q(1);
%     config(2).JointPosition = q(2);
%     config(3).JointPosition = q(3);
%     config(4).JointPosition = q(4);
%     config(5).JointPosition = q(5);
%     config(6).JointPosition = q(6);
end

figure(1);
show(robot1, q);
hold on;
ee_position = g(1:3, 4);
ray = [ee_position'; p(1:3)'];
plot3(ray(:, 1), ray(:, 2), ray(:, 3), 'r');
plot3(UGVposition(1), UGVposition(2), UGVposition(3), 'kx');


%%
theta = 0:0.01:(2 * pi);
plot3(0.25 * cos(theta) + 0.5, -0.25 * sin(theta) - 0.5, zeros(size(theta)), 'k');
plot3(cos(theta), sin(theta), zeros(size(theta)), 'b');



%%

r_UGV = 0.25;
dt = 0.001;
it = 0;
for t = 0:dt:(2 * pi)
    %disp(q);
    it = it + 1;
    UGVposition = r_UGV * [cos(t); -sin(t); 0; 0] + [0.5; -0.5; 0; 0];
    UGVvel = r_UGV * [-sin(t); -cos(t); 0];
    
    g = getTransform(robot1, q, 'ee', 'world');
    pray_hat = hat(g * [3.0*sqrt(2); 0; 0; 0]);
    gpray = g * [3.0*sqrt(2); 0; 0; 1];
    
    jac = geometricJacobian(robot1, q, 'ee');
    Jv = jac(4:6, :);
    Jw = jac(1:3, :);
    Jp = Jv - pray_hat * Jw;
    
    % jacobian for 1st task priority
    jac1 = Jp;
    r1dot = UGVvel;
    jac1inv = pinv(jac1);
    
    % jacobian for 2nd task priority
    jac2 = Jw;
    r2dot = zeros(3, 1);
    jac2inv = pinv(jac2);
    
    % jacobian for 3rd task priority
    jac3 = [1, 0, 0, 0, 0, 0];
    r3dot = 0;
    jac3inv = pinv(jac3);
    
    q1dot = jac1inv * r1dot;
    P1 = eye(6) - jac1inv * jac1;
    j2P1 = jac2 * P1;
    j2P1inv = pinv(j2P1);
    q2dot = q1dot + j2P1inv * (r2dot - jac2 * q1dot);
    P2 = P1 - j2P1inv * j2P1;
    
    j3P2 = jac3 * P2;
    if norm(j3P2) < 1e-5
        q3dot = q2dot;
    else
        j3P2inv = pinv(j3P2);
        q3dot = q2dot + P2 * j3P2inv * (r3dot - jac3 * q2dot);
    end
    
    if rem(it, 600) == 1
        show(robot1, q);
        hold on;
        ee_position = g(1:3, 4);
        ray = [ee_position'; gpray(1:3)'];
        plot3(ray(:, 1), ray(:, 2), ray(:, 3), 'r');
    end    
    q = q + dt * q3dot';
%     config(1).JointPosition = q(1);
%     config(2).JointPosition = q(2);
%     config(3).JointPosition = q(3);
%     config(4).JointPosition = q(4);
%     config(5).JointPosition = q(5);
%     config(6).JointPosition = q(6);

end