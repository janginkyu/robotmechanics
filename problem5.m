clear;
close all;
robot3 = importrobot('robot3.urdf');

UAVcenter = [0.5; -0.5; 1.0];

robot3.DataFormat = 'row';
robot3.Gravity = [0., 0., -9.81];
ik = inverseKinematics('RigidBodyTree', robot3);
UAVposition = [0.25; 0.0; 0.0] + UAVcenter;
initialPosition = [-UAVposition(1:2) / norm(UAVposition(1:2)); 1.0];
initialTheta = atan2(-UAVposition(2), -UAVposition(1));
initialR = [
    cos(initialTheta + pi), -sin(initialTheta + pi), 0; 
    sin(initialTheta + pi), cos(initialTheta + pi), 0;
    0, 0, 1
    ];
initialPose = [initialR, initialPosition; zeros(1, 3), 1];
c = [-0.4242   -1.7927    2.0884    2.0204   -1.0704    1.2478]
c = homeConfiguration(robot3);
[qinit, solInfo] = ik('ee', initialPose, ones(6, 1), c);
qinit

figure(2);
hold on;
axis equal;
show(robot3, qinit);
temp = 0:0.001:(2 * pi);
plot3(cos(temp), sin(temp), 1 * (ones(size(temp))));
plot3(0.25 * cos(temp) + UAVcenter(1), 0.25 * sin(temp) + UAVcenter(2), UAVcenter(3) * (ones(size(temp))));

% Dd = diag([1, ones(1, 5)]) * 100;
% Bd = diag(10 * ones(1, 6)) * 20;
% Kd = diag([100, 100 * ones(1, 5)]) * 10;
% Dd = diag([5, 5, 5, 2, 1, 1]) * 1;
% Bd = diag([100, 100, 100, 15.0, 1.50, 1.50]) * 10;
% Kd = diag([50, 50, 50, 0.5, 1.0, 1.0]) * 100;
Dd = diag([5, 5, 5, 1, 0.5, 0.5]) * 17;
Bd = diag([100, 100, 100, 5.50, 1.50, 1.50]) * 100;
Kd = diag([50, 50, 50, 1, 2.0, 2.0]) * 600;
fe = [0; 0; 0; 0; 0; 0];
%%
w = 1;
it = 0;
q = qinit;
qdotmeas = zeros(1, 6);
qmeas = meas(q);
qmeashist = qmeas;
Jmeashist = cylJacobian(robot3, qmeas);
timestamp = 0;
timewindow = 16;
qdot = zeros(1, 6);
qhist = [];
qdothist = [];
uhist = [];
tt = 0:0.001:6.29;
thetadeshist = [];
thetadotdeshist = [];
thetaddotdeshist = [];
qdotmeashist = [];
fhist = [];
xphist = [];
q_hist = [];
xp0 = zeros(6, 1);
xdotr0 = zeros(6, 1);
for t = tt
    it = it + 1;
    % measurement of q and qdot, x and xdot
    qmeas = meas(q);
    Jmeas = cylJacobian(robot3, qmeas);
    Mmeas = massMatrix(robot3, qmeas);
    if size(qmeashist, 1) < timewindow
        qmeashist = unwrap([qmeas; qmeashist]);
        timestamp = [t; timestamp];
        Jmeashist = [Jmeas, Jmeashist];
    else
        qmeashist = unwrap([qmeas; qmeashist(1:(timewindow - 1), :)]);
        timestamp = [t; timestamp(1:(timewindow - 1), :)];
        Jmeashist = [Jmeas, Jmeashist(:, 1:(6 * (timewindow - 1)))];
    end
    if size(qmeashist, 1) < 3
        qdotmeas = zeros(1, 6);
        Jdotmeas = zeros(6, 6);
    else
        qdotmeas = (qmeashist(1, :) - qmeashist(end, :)) / (timestamp(1) - timestamp(end));
        Jdotmeas = (Jmeas(:, 1:6) - Jmeas(:, (end-5):(end))) / (timestamp(1) - timestamp(end));
    end
    gmeas = getTransform(robot3, qmeas, 'ee', 'world');
    pmeas = gmeas(1:3, 4);
    xdotmeas = Jmeas * qdotmeas';
    %x = cylf(robot2, q);
    %xdot = cylJacobian(robot2, q) * qdot;
    
    % desired trajectory
    UAVposition = 0.25 * [cos(w * t); -sin(w * t); 0] + UAVcenter;
    UAVvelocity = w * 0.25 * [-sin(w * t); -cos(w * t); 0];
    theta_UAV = atan2(UAVposition(2), UAVposition(1));
    rrdot_UAV = dot(UAVposition(1:2), UAVvelocity(1:2));
    r2_UAV = dot(UAVposition(1:2), UAVposition(1:2));
    r_UAV = sqrt(r2_UAV);
    rdot_r_UAV = rrdot_UAV / r2_UAV;
    thetadot_UAV = -cos(theta_UAV + w * t) * 0.25 / r_UAV;
    thetaddot_UAV = -thetadot_UAV * rdot_r_UAV + (sin(theta_UAV + w * t) * (thetadot_UAV + w * 1)) * 0.25 / r_UAV;
    
    rdes = 1;
    rdotdes = 0;
    rddotdes = 0;
    thetades = theta_UAV + pi;
    thetadotdes = thetadot_UAV;
    thetaddotdes = thetaddot_UAV;
    hdes = 1;
    hdotdes = 0;
    hddotdes = 0;
    
    xdes = [1; thetades; 1];
    xdotdes = [0; thetadotdes; 0];
    xddotdes = [0; thetaddotdes; 0];
    Rdes = [
        cos(thetades + pi), -sin(thetades + pi), 0;
        sin(thetades + pi), cos(thetades + pi), 0;
        0, 0, 1];
    omegades = [0; 0; thetadotdes];
    alphades = [0; 0; thetaddotdes];
    xdotdes = [omegades; xdotdes];
    xddotdes = [alphades; xddotdes];
    
    if t == 0
        xdotp = xdotr0 - xdotdes;
        xp = xp0;
    end
    
    fe = zeros(6, 1);
    fe(4, 1) = -500 * xp(4, 1);
    xddotp = Dd \ (fe - Bd * xdotp - Kd * xp);
    xdotp = xdotp + 0.001 * xddotp;
    xp = xp + 0.001 * xdotp;
    xddotr = xddotp + xddotdes;
    xdotr = xdotp + xdotdes;
    xr = xp;
    xr(4:6, 1) = xr(4:6, 1) + xdes;

    Rmeas = gmeas(1:3, 1:3);
    dR = Rmeas \ Rdes;
    dw = (dR - dR');
    dw = [dw(3, 2), dw(1, 3), dw(2, 1)];
    if dw(2) > 0
        dtheta = atan2(norm(dw), (trace(dR) - 1));
    else
        dtheta = atan2(-norm(dw), (trace(dR) - 1));
    end
    dw = 2 * dw / norm(dw) * dtheta;
    
    xtild = zeros(6, 1);
    xtild(4:6, 1) = cylf(robot3, qmeas) - xr(4:6, 1);
    xtild(1:3, 1) = dw' - xp(1:3, 1);
    xdottild = xdotmeas - xdotr;
    
    %fe = 0.1 * [0, 0, 0, -10, 0, 0]';
    
%     posdes = [cos(thetades); sin(thetades); hdes];
%     [q, solInfo] = ik('ee', [Rdes, posdes; zeros(1, 3), 1], ones(6, 1), q);
%     
    % control loop
    u = velocityProduct(robot3, qmeas, qdotmeas)' - Jmeas' * Mmeas * (Jmeas \ Jdotmeas) * qdotmeas' + gravityTorque(robot3, qmeas)' - Jmeas' * fe...
        + Mmeas * (Jmeas \ ( ...
        xddotr - Dd \ (Bd * xdottild + Kd * xtild)...
        ));
    
%     if rem(it, 100) == 1
%         t
%         dw
%         dR
%     end
    if rem(it, 600) == 1
        show(robot3, q);
        hold on;
        %ee_position = g(1:3, 4);
        %ray = [ee_position'; gpray(1:3)'];
        %plot3(ray(:, 1), ray(:, 2), ray(:, 3), 'r');
        temp = [UAVposition, gmeas(1:3, 4)];
        plot3(temp(1, :), temp(2, :), temp(3, :), 'r');
        temp = [UAVposition, [cos(thetades); sin(thetades); 1]];
        plot3(temp(1, :), temp(2, :), temp(3, :), 'b');
        it
    end
    thetadeshist = [thetadeshist, thetades];
    thetadotdeshist = [thetadotdeshist, thetadotdes];
    thetaddotdeshist = [thetaddotdeshist, thetaddotdes];
%     
%     M = Rp + 0.001 * hat(omegap) * Rp;
%     Rp = M * (sqrtm(M' * M) \ eye(3));
%     rp = rp + 0.001 * rdotp;
%     

    qddot = forwardDynamics(robot3, q, qdot, u');
    qdot = qdot + 0.001 * qddot;
    qold = q;
    q = wrapToPi(q + 0.001 * qdot);
    g = getTransform(robot3, q, 'ee');
    p = g(1:3, 4);
    r = norm(p(1:2));
    theta = atan2(p(2), p(1));
    % numerical integration
    if r > 1
        %lambda = 0;
    else
        g(1:2, 4) = p(1:2, 1) / norm(p(1:2, 1));
        [qnew, solInfo] = ik('ee', g, ones(6, 1), q);
        qdot = wrapToPi(qnew - qold) / 0.001;
        q = wrapToPi(qnew);
%         qddot1 = forwardDynamics(robot2, q, qdot, u');
%         qddot2 = forwardDynamics(robot2, q, qdot, u', 0.01 * [zeros(7, 6); [0, 0, 0, cos(theta), sin(theta), 0]]);
%         qdot1 = qdot + 0.001 * qddot1;
%         q1 = wrapToPi(q + 0.001 * qdot1);
%         qdot2 = qdot + 0.001 * qddot2;
%         q2 = wrapToPi(q + 0.001 * qdot2);
%         dp = jac * wrapToPi((q2' - q1'));
%         lambda = (1 - r^2) / (2 * dot(p(1:2), dp));
%         qddot = forwardDynamics(robot2, q, qdot, u', lambda * 0.01 * [zeros(7, 6); [0, 0, 0, cos(theta), sin(theta), 0]]);
%         qdot = qdot + 0.001 * qddot;
%         q = q + 0.001 * qdot;
    end
    
    qhist = [qhist; qmeas];
    q_hist = [q_hist; q];
    qdothist = [qdothist; qdot];
    qdotmeashist = [qdotmeashist; qdotmeas];
    uhist = [uhist; u'];
    xphist = [xphist; xp'];
    %fhist = [fhist; lambda * 0.01];
    
%     qdotmeas
%     qmeas
%     q
%     qdot
    
end

