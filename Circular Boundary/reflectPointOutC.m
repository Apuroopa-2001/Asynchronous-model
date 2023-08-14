function [x_2, y_2, x_4, y_4, reflec_angle] = reflectPointOutC(x_0, y_0, x_1, y_1, circle_R)
    a = (x_1 - x_0)^2 + (y_1 - y_0)^2;
    b = 2 * (y_1 - y_0) * ((x_1 - x_0)*y_0 - (y_1 - y_0)*x_0);
    c = (x_0*(y_1 - y_0) - y_0*(x_1 - x_0))^2 - ((x_1 - x_0)^2 * circle_R^2);
    b2_4ac = b^2 - 4*a*c;
    rootof_b2_4ac = sqrt(b2_4ac);
    x_2_a = (-b + rootof_b2_4ac) / (2*a);
    y_2_a = sqrt(circle_R^2 - x_2_a^2);
    y_2_a_1 = y_2_a;
    y_2_a_2 = -y_2_a;
    cos_2a10v = round(((x_1 - x_0)*(x_2_a - x_0) + (y_1 - y_0)*(y_2_a_1 - y_0))/(sqrt((x_1 - x_0)^2 + (y_1 - y_0)^2)*sqrt((x_2_a - x_0)^2 + (y_2_a_1 - y_0)^2)), 8);
    diff_theta_2a10_v = acos(cos_2a10v);
    cos_2a20v = round(((x_1 - x_0)*(x_2_a - x_0) + (y_1 - y_0)*(y_2_a_2 - y_0))/(sqrt((x_1 - x_0)^2 + (y_1 - y_0)^2)*sqrt((x_2_a - x_0)^2 + (y_2_a_2 - y_0)^2)), 8);
    diff_theta_2a20_v = acos(cos_2a20v);
    if diff_theta_2a10_v <= diff_theta_2a20_v
        y_2_a = y_2_a_1;
    else
        y_2_a = y_2_a_2;
    end
    
    x_2_b = (-b - rootof_b2_4ac) / (2*a);
    y_2_b = sqrt(circle_R^2 - x_2_b^2);
    y_2_b_1 = y_2_b;
    y_2_b_2 = -y_2_b;
    cos_2b10v = round(((x_1 - x_0)*(x_2_b - x_0) + (y_1 - y_0)*(y_2_b_1 - y_0))/(sqrt((x_1 - x_0)^2 + (y_1 - y_0)^2)*sqrt((x_2_b - x_0)^2 + (y_2_b_1 - y_0)^2)), 8);
    diff_theta_2b10_v = acos(cos_2b10v);
    cos_2b20v = round(((x_1 - x_0)*(x_2_b - x_0) + (y_1 - y_0)*(y_2_b_2 - y_0))/(sqrt((x_1 - x_0)^2 + (y_1 - y_0)^2)*sqrt((x_2_b - x_0)^2 + (y_2_b_2 - y_0)^2)), 8);
    diff_theta_2b20_v = acos(cos_2b20v);
    if diff_theta_2b10_v <= diff_theta_2b20_v
        y_2_b = y_2_b_1;
    else
        y_2_b = y_2_b_2;
    end
    
    if sqrt((x_1 - x_2_a)^2 + (y_1 - y_2_a)^2) <= sqrt((x_1 - x_2_b)^2 + (y_1 - y_2_b)^2)
        x_2 = x_2_a;
        y_2 = y_2_a;
    else
        x_2 = x_2_b;
        y_2 = y_2_b;
    end
    
    x_3 = x_2 - ((y_2 * (y_1*x_2 - y_2*x_1))/circle_R^2);
    y_3 = ((x_2*(y_1*x_2 - y_2*x_1))/circle_R^2) + y_2;
    x_4 = (2*x_3) - x_1;
    y_4 = (2*y_3) - y_1;
    reflec_angle = atan2((y_4 - y_2), (x_4 - x_2));
end
