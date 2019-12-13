% in the order of omega, r, theta, h
function wJ = cylJacobian(robot, q)

g = getTransform(robot, q, 'ee', 'world');
jacobian = geometricJacobian(robot, q, 'ee');
position = g(1:3, 4);
x = position(1);
y = position(2);
r2 = dot(position(1:2), position(1:2));
r = sqrt(r2);
wJ = [eye(3), zeros(3, 3); zeros(3, 3), [x / r, y / r, 0; -y / r2, x / r2, 0; 0, 0, 1]] * jacobian;