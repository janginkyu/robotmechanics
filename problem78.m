clear;
close all;
robot4 = importrobot('robot4.urdf');

UAVcenter = [0.5; -0.5; 0];

robot4.DataFormat = 'column';
robot4.Gravity = [0., 0., -9.81];
ik = inverseKinematics('RigidBodyTree', robot4);

RFposition = [-3.15; -3.25; 0.3];
xdirdes = [0; 0; -1];

maxiter = 1000;

q = [pi, 0, 0, 0, 0, 0]';

w = 1;
lambda = 0.4;
for it = 1:maxiter
    g = getTransform(robot4, q, 'ee', 'world');
    p = g(1:3, 4);
    
    jac = geometricJacobian(robot4, q, 'ee');
    
    Jp = jac(4:6, :);
    Jw = jac(1:3, :);
    
    xdir = g(1:3, 1);
    Jx = -hat(xdir) * Jw;
    cost = w * norm(xdir - xdirdes)^2 + norm(p - RFposition)^2;
    Jcost = w * (xdir - xdirdes)' * Jx + (p - RFposition)' * Jp;
    q = q - lambda * Jcost';
    
    if cost < 1e-7
        break;
    end
end
qpickup = q;
figure(43);
show(robot4, qpickup);
hold on;

figure(1);
%show(robot4, q);
%axis equal;
%hold on;

UAVposition = [0.25; 0; 0.5] + UAVcenter;
p0 = [4; 0; 0; 1];

q = [
    1.4282
   -1.1627
    1.9559
    0.2964
   -1.3369
    0.9012];

lambda = 0.04;
for it = 1:maxiter
    % current position of ray end
    g = getTransform(robot4, q, 'ee', 'world');
    p = g * p0;
    p0dir = p0;
    p0dir(4) = 0;
    p0_hat = hat(g * p0dir);
    
    % classical jacobian of ee w.r.t. world
    jac = geometricJacobian(robot4, q, 'ee');
    
    Jv = jac(4:6, :);
    Jw = jac(1:3, :);
    
    Jp = Jv - p0_hat * Jw;
    q = q - lambda * ((p(1:3, 1) - UAVposition)' * Jp)';
    
    if norm(p(1:3, 1) - UAVposition) < 1e-5
        break;
    end
end
%show(robot4, q);
axis equal;
hold on;

qstart = q;
figure(43);
show(robot4, qstart);
hold on;
show(robot4, homeConfiguration(robot4));
hold on;
t = 0:0.001:7;
plot3(cos(t), sin(t), zeros(size(t)), 'b');
%%

figure(3);
delq1 = wrapToPi(qpickup);
delq2 = wrapToPi(qstart - qpickup);
t0 = 0;
qmeashist = [];
B = massMatrix(robot4, qpickup) * 20;
K = massMatrix(robot4, qpickup) * 100;
it = 0;
timewindow = 6;
timestamp = [];
qdot = 0;
q = homeConfiguration(robot4);
for t = 0:0.001:12
    it = it + 1;
    if t >= 6
        t0 = 6;
    end
    posbas = 0.5 * (1 - cos((t - t0) * pi / 6));
    velbas = 0.5 * pi / 6 * sin((t - t0) * pi / 6);
    accbas = 0.5 * (pi / 6)^2 * cos((t - t0) * pi / 6);
    if t < 6
        qdes = posbas * delq1;
        qdotdes = velbas * delq1;
        qddotdes = accbas * delq1;
    else
        qdes = posbas * delq2 + qpickup;
        qdotdes = velbas * delq2;
        qddotdes = accbas * delq2;
    end
    qmeas = meas(q);
    Mmeas = massMatrix(robot4, qmeas);
    if size(qmeashist, 1) < timewindow
        qmeashist = unwrap([qmeas'; qmeashist]);
        timestamp = [t; timestamp];
    else
        qmeashist = unwrap([qmeas'; qmeashist(1:(timewindow - 1), :)]);
        timestamp = [t; timestamp(1:(timewindow - 1), :)];
    end
    if size(qmeashist, 1) < 3
        qdotmeas = zeros(6, 1);
    else
        qdotmeas = (qmeashist(1, :) - qmeashist(end, :))' / (timestamp(1) - timestamp(end));
    end
    e = wrapToPi(qmeas - qdes);
    edot = wrapToPi(qdotmeas - qdotdes);
    tau = Mmeas * qddotdes + velocityProduct(robot4, qmeas, qdotmeas) + gravityTorque(robot4, qmeas) - (B * edot + K * e);
    
    qddot = forwardDynamics(robot4, qmeas, qdotmeas, tau);
    qdot = qdot + 0.001 * qddot;
    q = q + 0.001 * qdot;
    if mod(it, 600) == 1
        if it < 6500
            figure(33);
            show(robot4, q);
            hold on;
        end
        if it > 5500
            figure(44);
            show(robot4, q);
            hold on;
        end
        it
    end
end

%%

eye23 = [0, 1, 0; 0, 0, 1];
qdes = qstart;
q = qstart;
qdotdes = zeros(6, 1);
qdestraj = [];
qddotdes = zeros(6, 1);
qdotdestraj = [];
qddotdestraj = [];
for t = 0:0.001:10
    
    % UAV Position
    
    UAVposition = UAVcenter + [0.25 * cos(t); -0.25 * sin(t); 0.5 + 0.3 * sin(0.2 * pi * t)];
    UAVvelocity = [-0.25 * sin(t); -0.25 * cos(t); 0.06 * pi * cos(0.2 * pi * t)];
    
    qdestraj = [qdestraj; qdes'];
    qdotdestraj = [qdotdestraj; qdotdes'];
    g = getTransform(robot4, qdes, 'ee', 'world');
    pEE = g(1:3, 4);
    REE = g(1:3, 1:3);
    qdotdesold = qdotdes;
    
    J1 = eye23 * REE' * ( ...
        Jp + hat(pEE - UAVposition) * Jw...
        );
    rdot1 = eye23 * REE' * UAVvelocity;
    J2 = pEE' * Jp;
    rdot2 = 0;
    
    P1 = eye(6) - pinv(J1) * J1;
    qdotdes = pinv(J1) * rdot1 + pinv(J2 * P1) * (rdot2 - J2 * pinv(J1) * rdot1);
    qdes = qdes + qdotdes * 0.001;
end
qddotdestraj = (qdotdestraj(3:end, :) - qdotdestraj(2:(end-1), :)) / 0.001;
qddotdestraj = [zeros(1, 6); qddotdestraj; zeros(1, 6)];
%%
umass = importrobot('umass.urdf');
umass.DataFormat = 'column';
umass.Gravity = [0, 0, -9.81];
thetahat = 0;
qmeashist = meas(q)';
timewindow = 6;
timestamp = [];
qdotdes = zeros(6, 1);
qdes = q;
thetaest = 0;
B = 20 * massMatrix(robot4, q) * eye(6);
K = 100 * massMatrix(robot4, q) * eye(6);
lambda = 50;
u = zeros(6, 1);
thetatildint = 0;
mass = 10;
qdot = zeros(6, 1);
it = 0;
thetaesthist = [];
thetatildhist = [];
ehist = [];
gamma = 0.2;
delta = 1.5;
for t = 0:0.001:8.001
    it = it + 1;
    % measurement of q and qdot, x and xdot
    qmeas = meas(q);
    Jmeas = geometricJacobian(robot4, qmeas, 'ee');
    Jpmeas = Jmeas(4:6, :);
    Jwmeas = Jmeas(1:3, :);
    M0meas = massMatrix(robot4, qmeas);
    Mumeas = massMatrix(umass, qmeas);
    if size(qmeashist, 1) < timewindow
        qmeashist = unwrap([qmeas'; qmeashist]);
        timestamp = [t; timestamp];
    else
        qmeashist = unwrap([qmeas'; qmeashist(1:(timewindow - 1), :)]);
        timestamp = [t; timestamp(1:(timewindow - 1), :)];
    end
    if size(qmeashist, 1) < 3
        qdotmeas = zeros(6, 1);
    else
        qdotmeas = (qmeashist(1, :) - qmeashist(end, :))' / (timestamp(1) - timestamp(end));
    end
    
    % UAV Position
    
    UAVposition = UAVcenter + [0.25 * cos(t); -0.25 * sin(t); 0.5 + 0.3 * sin(0.2 * pi * t)];
    UAVvelocity = [-0.25 * sin(t); -0.25 * cos(t); 0.06 * pi * cos(0.2 * pi * t)];
    
    % Track qdes
    
    qdes = qdestraj(it, :)';
    qdotdes = qdotdestraj(it, :)';
    qddotdes = qddotdestraj(it, :)';
    
    edot = qdotmeas - qdotdes;
    e = wrapToPi(qmeas - qdes);
    pdes = getTransform(robot4, qdes, 'ee');
    pdes = pdes(1:3, 4);
    
    Mest = M0meas + thetaest * Mumeas;
    ff0 = inverseDynamics(robot4, qmeas, qdotmeas, qddotdes);
    ffu = inverseDynamics(umass, qmeas, qdotmeas, qddotdes);
    ffest = ff0 + thetaest * ffu;
    fbest = (B * edot + K * e);
    fbu = (B * edot + K * e);
    
%     thetatild = dot(u - (ffest - fbest), ffu - fbu) / dot(ffu - fbu, ffu - fbu);
%     if it > 1
%         thetatildint = thetatildint + thetatild * 0.001;
%     end
    u = ffest - fbest; % - lambda * (ffu - fbu) * thetatildint;
    
    %thetaest = thetaest - lambda^2 * thetatildint * 0.001;
    %thetaest = thetaest + thetatild * 0.01;
    thetatilddot = dot(edot + delta * e, ffu) / gamma;
    thetaest = thetaest - 0.001 * thetatilddot;
    
    qddot = (massMatrix(robot4, q) + mass * massMatrix(umass, q)) \ ( ...
        u - velocityProduct(robot4, q, qdot) - mass * velocityProduct(umass, q, qdot) ...
        - gravityTorque(robot4, q) - mass * gravityTorque(umass, q) ...
    );
    qdot = qdot + qddot * 0.001;
    q = q + qdot * 0.001;
    
    if 0 %mod(it, 1000) == 0
        show(robot4, q);
        hold on;
        g_ = getTransform(robot4, q, 'ee');
        p_ = g_(1:3, 4);
        temp = [UAVposition'; p_'];
        plot3(temp(:, 1), temp(:, 2), temp(:, 3), 'r');
        it
    end
    thetaesthist = [thetaesthist; thetaest];
    ehist = [ehist; norm(e)];
    %thetatildhist = [thetatildhist; thetatild];
end
%%
t = (0:0.001:8);
UAVposition = UAVcenter + [0.25 * cos(t); -0.25 * sin(t); 0.5 + 0.3 * sin(0.2 * pi * t)];
plot3(UAVposition(1, :), UAVposition(2, :), UAVposition(3, :), 'k');
hold on;
plot3(cos(t), sin(t), zeros(size(t)), 'b');
axis equal;
grid on;
