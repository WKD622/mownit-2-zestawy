function draw_points_1()
figure; hold on;
axis ([0 100 0 3])
x_double_2_f = load('out/x_double_2_f');
plot(x_double_2_f(:,1),x_double_2_f(:,2), 'o-r', 'MarkerSize', 4);
x_float_2_f = load('out/x_float_2_f');
plot(x_float_2_f(:,1),x_float_2_f(:,2), 'o-b', 'MarkerSize', 4);
hold off;