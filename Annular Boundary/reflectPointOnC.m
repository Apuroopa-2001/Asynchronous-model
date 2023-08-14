function reflec_angle = reflectPointOnC(x_0, y_0, x_1, y_1, x_2, y_2, circle_R)
    x_3 = x_2 - ((y_2 * (y_1*x_2 - y_2*x_1))/circle_R^2);
    y_3 = ((x_2*(y_1*x_2 - y_2*x_1))/circle_R^2) + y_2;
    x_4 = (2*x_3) - x_1;
    y_4 = (2*y_3) - y_1;
    reflec_angle = atan2((y_4 - y_2),(x_4 - x_2));
end
