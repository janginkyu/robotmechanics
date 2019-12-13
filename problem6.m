clear;
close all;
robot3 = importrobot('robot3.urdf');

UAVcenter = [0.5; -0.5; 0.8];

robot3.DataFormat = 'row';
robot3.Gravity = [0., 0., -9.81];
ik = inverseKinematics('RigidBodyTree', robot3);
UAVposition = [0.25; 0.0; 0.0] + UAVcenter;
initialPosition = [-UAVposition(1:2) / norm(UAVposition(1:2)); 0.8];
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
plot3(cos(temp), sin(temp), 0.8 * (ones(size(temp))));
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
lambdahist = [];
Jphimeashist = [];
Ameashist = [];
xp0 = zeros(6, 1);
xdotr0 = zeros(6, 1);
eye65 = [
    1, 0, 0, 0, 0;
    0, 1, 0, 0, 0;
    0, 0, 1, 0, 0;
    0, 0, 0, 0, 0;
    0, 0, 0, 1, 0;
    0, 0, 0, 0, 1
    ];
eye56 = eye65';
u = zeros(6, 1);
lambdaintmeas = 0;

Kv = 100;
Kp = 100;
Kf = 1000;
lambdades = -10;

for t = tt
    it = it + 1;
    % measurement of q and qdot, x and xdot
    qmeas = meas(q);
    Jmeas = cylJacobian(robot3, qmeas);
    Jphimeas = Jmeas \ eye65;
    Mmeas = massMatrix(robot3, qmeas);
    Ameas = [0, 0, 0, 1, 0, 0] * Jmeas;
    if size(qmeashist, 1) < timewindow
        qmeashist = unwrap([qmeas; qmeashist]);
        timestamp = [t; timestamp];
        Jmeashist = [Jmeas, Jmeashist];
        Jphimeashist = [Jphimeas, Jphimeashist];
        Ameashist = [Ameas; Ameashist];
    else
        qmeashist = unwrap([qmeas; qmeashist(1:(timewindow - 1), :)]);
        timestamp = [t; timestamp(1:(timewindow - 1), :)];
        Jmeashist = [Jmeas, Jmeashist(:, 1:(6 * (timewindow - 1)))];
        Jphimeashist = [Jphimeas, Jphimeashist(:, 1:(5 * (timewindow - 1)))];
        Ameashist = [Ameas; Ameashist(1:(timewindow - 1), :)];
    end
    if size(qmeashist, 1) < 3
        qdotmeas = zeros(1, 6);
        Jdotmeas = zeros(6, 6);
        Jphidotmeas = zeros(6, 5);
        Adotmeas = zeros(1, 6);
    else
        qdotmeas = (qmeashist(1, :) - qmeashist(end, :)) / (timestamp(1) - timestamp(end));
        Jphidotmeas = (Jphimeas(:, 1:5) - Jphimeas(:, (end-4):(end))) / (timestamp(1) - timestamp(end));
        Adotmeas = (Ameashist(1, :) - Ameashist(end, :)) / (timestamp(1) - timestamp(end));
    end
    gmeas = getTransform(robot3, qmeas, 'ee', 'world');
    pmeas = gmeas(1:3, 4);
    xdotmeas = Jmeas * qdotmeas';
    phidotmeas = eye56 * xdotmeas;
    %x = cylf(robot2, q);
    %xdot = cylJacobian(robot2, q) * qdot;
    lambdameas = (Ameas * (Mmeas \ Ameas')) \ (Adotmeas * qdotmeas' + Ameas * (Mmeas \ ( ...
        u - velocityProduct(robot3, qmeas, qdotmeas)' - gravityTorque(robot3, qmeas)'...
    )));
    lambdaintmeas = lambdaintmeas + 0.001 * (lambdameas - lambdades);
    
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
    
    xdes = [1; thetades; 0.8];
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
    phiddotdes = eye56 * xddotdes;

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
    xtild(4:6, 1) = cylf(robot3, qmeas) - xdes;
    xtild(1:3, 1) = dw';
    xdottild = xdotmeas - xdotdes;
    
    phitild = eye56 * xtild;
    phidottild = eye56 * xdottild;
    
    
    %fe = 0.1 * [0, 0, 0, -10, 0, 0]';
    
%     posdes = [cos(thetades); sin(thetades); hdes];
%     [q, solInfo] = ik('ee', [Rdes, posdes; zeros(1, 3), 1], ones(6, 1), q);
%     
%     % control loop
%     u = velocityProduct(robot3, qmeas, qdotmeas)' - Jmeas' * Mmeas * (Jmeas \ Jdotmeas) * qdotmeas' + gravityTorque(robot3, qmeas)' - Jmeas' * fe...
%         + Mmeas * (Jmeas \ ( ...
%         xddotr - Dd \ (Bd * xdottild + Kd * xtild)...
%         ));
    u = Mmeas * Jphimeas * (phiddotdes - Kv * phidottild - Kp * phitild) ...
        + velocityProduct(robot3, qmeas, qdotmeas)' + gravityTorque(robot3, qmeas)' + Mmeas * Jphidotmeas * phidotmeas...
        + Ameas' * (lambdades - Kf * lambdaintmeas);
    
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
        temp = [UAVposition, [cos(thetades); sin(thetades); 0.8]];
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
    Jold = cylJacobian(robot3, qold);
    Jnew = cylJacobian(robot3, q);
    Jphi = Jnew \ eye65;
    M = massMatrix(robot3, q);
    A = [0, 0, 0, 1, 0, 0] * Jnew;
    Aold = [0, 0, 0, 1, 0, 0] * Jold;
    Adot = (A - Aold) / 0.001;
    lambda = (A * inv(M) * A') \ (...
        Adot * qdot' + A * forwardDynamics(robot3, q, qdot, u')');
    
    qhist = [qhist; qmeas];
    q_hist = [q_hist; q];
    qdothist = [qdothist; qdot];
    qdotmeashist = [qdotmeashist; qdotmeas];
    uhist = [uhist; u'];
    lambdahist = [lambdahist; lambda];
    
    
    
%     xphist = [xphist; xp'];
    %fhist = [fhist; lambda * 0.01];
    
%     qdotmeas
%     qmeas
%     q
%     qdot
    
end
%%
figure(2310);
plot((0:(size(lambdahist, 1) - 1)) / 1000, lambdahist, 'k');
xlabel("time [s]");
ylabel("\lambda [N]");
grid on;

