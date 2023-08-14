close all
clear
clc

%%
load('try.mat') % Load data

iter = 1; % iteration 

pos_t = pos_t(:,:,:,iter);
theta_t = theta_t(:,:,iter);

for t = 1:10:n_iter

    vel_x = cos(theta_t(:,t));
    vel_y = sin(theta_t(:,t));
    pos_x = pos_t(:,1,t);
    pos_y = pos_t(:,2,t);

    quiver(pos_x, pos_y, vel_x, vel_y, 0.2, 'LineWidth', 2.5, 'ShowArrowHead','off',...
        'Color', '#B0E0E6')

    hold all

    plot(pos_x, pos_y, '.', 'Color', '#00A693', 'MarkerSize', 25)
    angles=linspace(0,2*pi,200);
     
    % Generate x-coordinates.
    x=circle_R*cos(angles);
     
    % Generate y-coordinate.
    y=circle_R*sin(angles);
     
    % plot the circle.
    plot(x,y);
    hold off

    axis('equal')
    axis([-circle_R-1 circle_R+1 -circle_R-1 circle_R+1])

    drawnow('limitrate')

end
