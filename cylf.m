function x = cylf(robot, q)

g = getTransform(robot, q, 'ee', 'world');
position = g(1:3, 4);
x_ = position(1);
y_ = position(2);
h = position(3);

r = norm([x_, y_]);
theta = atan2(y_, x_);

x = [r; theta; h];
