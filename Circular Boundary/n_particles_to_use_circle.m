function [theta_t, pos_t, sum_int] = n_particles_to_use_circle (n, r_spon, r_align, sigma_theta, dt, n_iter, ...
    k_alg, S0, circle_R, int_rad)

% INITIAL CONDITIONS - Putting agents in a lattice

radial_pos(:) = rand(1, n) * circle_R;
radial_angle(:) = rand(1, n) * 2 * pi - pi;
pos(:,1) = radial_pos .* cos(radial_angle);
pos(:,2) = radial_pos .* sin(radial_angle);


theta = 2*pi*rand(n,1); % orientation of agents

pos_t = zeros(n,2,n_iter); % Stores positions of agents over the simulation
theta_t = zeros(n,n_iter); % Stores orientation of agents over the simulation

pos_t(:,:,1) = pos; % position at time t = 1
theta_t(:,1) = theta; % orientation at time t = 1

% position and orientation at time t-1
pos_t_1 = pos;
theta_t_1 = theta;

theta_d = theta; % Desired heading angle

% Reaction rates
te_spon = (1/r_spon)*log(1./rand(n,1)); % time at which spontaneous reaction happens
te_align = (1/r_align)*log(1./rand(n,1)); % time at which alignment reaction happens
sum_int = zeros(n_iter,1);

for t = 2:n_iter

    if rem(t, 1e4) == 0
        disp(t*dt);
    end

    sum_spon = (te_spon - dt) <= dt;
    sum_alg = (te_align - dt) <= dt;
    sum_int(t) = sum(sum_spon | sum_alg);

    for i = 1:n
        
        e_spon = te_spon(i) <= dt; % Check if spontaneous reaction happens at this time step
        e_align = te_align(i) <= dt; % Check if alignment reaction happens at this time step

        te_spon(i) = te_spon(i) - dt;
        te_align(i) = te_align(i) - dt;
        
        % If the interaction happened at this time step, chose a next time
        % point where the interaction will happen.
        if te_spon(i) < 0
            te_spon(i) = (1/r_spon)*log(1/rand());
        end

        if te_align(i) < 0
            te_align(i) = (1/r_align)*log(1/rand());
        end

        % Check if any event is happening in this time interval

        %if no event, continue to move in the desired direction of previous step

        if (e_spon + e_align) == 0

            theta_d(i) = theta_d(i);

        else % in case of an event 

            %dist_x = pos_t_1(:,1) - repmat(pos_t_1(i,1),n,1); % distance between x coordinate of i and x coordinates of other agents 
            %dist_x = dist_x - (round(dist_x/box_length))*box_length; % taking the shortest distance because of periodic boundary conditions
            %dist_y = pos_t_1(:,2) - repmat(pos_t_1(i,2),n,1); % similarly for y
            %dist_y = dist_y - (round(dist_y/box_length))*box_length;
            dist_vec = pos_t_1 - repmat(pos_t_1(i,:),n,1);
            dist_mag = sqrt(dist_vec(:,1).^2 + dist_vec(:,2).^2);
            %dist_mag = sqrt(dist_x.^2 + dist_y.^2); %magnitude of rij
            dist_mag(i) = int_rad + 1;

            % Spontaneous rotation

            % theta_d_spon = theta_t_1(i) + random(trunc_dis_ang);
            theta_d_spon = theta_t_1(i) + (-sigma_theta + 2*sigma_theta*rand());
            % theta_d_spon = (-sigma_theta + 2*sigma_theta*rand());

            % alignmnet

            neigh_align = find(dist_mag < int_rad);
            if isempty(neigh_align) == 0
                neigh_align = neigh_align(randperm(length(neigh_align),min(length(neigh_align), k_alg))); %Picking k_alg random neighbor to align
                theta_d_align = mean([cos(theta_t_1(neigh_align)) sin(theta_t_1(neigh_align))],1);
                theta_d_align = atan2(theta_d_align(1,2), theta_d_align(1,1));
            else
                theta_d_align = theta_d(i);
            end

            % resultant vector

            theta_d(i) = atan2(e_spon*sin(theta_d_spon) + e_align*sin(theta_d_align),...
                e_spon*cos(theta_d_spon) + e_align*cos(theta_d_align));

            if theta_d(i) < 0
                theta_d(i) = theta_d(i) + 2*pi;
            elseif theta_d(i) > 2*pi
                theta_d(i) = theta_d(i) - 2*pi;
            end

        end

        theta(i) = theta_d(i);
        init_pos(:) = pos(i,:); %position pre-update
        pos(i,:) = pos(i,:) + S0 * dt *[cos(theta(i)) sin(theta(i))]; % update position
        
        % ensure that agents stay within the box
        %pos(i,1) = pos(i,1) - floor(pos(i,1)/box_length)*box_length;
        %pos(i,2) = pos(i,2) - floor(pos(i,2)/box_length)*box_length;
        r_x_1_y_1 = sqrt(pos(i,1)^2 + pos(i,2)^2); % Calculating distance of the next position from the center
        
        if r_x_1_y_1 == circle_R % Introducing reflective boundary
            x_0 = init_pos(1);
            y_0 = init_pos(2);
            x_1 = init_pos(1) + 1.5 * dt * S0 * cos(theta(i));
            y_1 = init_pos(2) + 1.5 * dt * S0 * sin(theta(i));
            x_2 = pos(i,1);
            y_2 = pos(i,2);
            theta(i) = reflectPointOnC(x_0, y_0, x_1, y_1, x_2, y_2, circle_R);
            
        elseif r_x_1_y_1 > circle_R % Introducing reflective boundary
            x_0 = init_pos(1);
            y_0 = init_pos(2);
            x_1 = pos(i,1);
            y_1 = pos(i,2);
            x_2 = 0;
            y_2 = 0;
            dist_reflec_point = r_x_1_y_1;
            index = 0;
            
            while dist_reflec_point >= circle_R && index == 0
                if dist_reflec_point > circle_R
                    [x_0, y_0, x_1, y_1, theta(i)] = reflectPointOutC(x_0, y_0, x_1, y_1, circle_R);
                    dist_reflec_point = sqrt(x_1^2 + y_1^2);
                elseif dist_reflec_point == circle_R
                    x_2 = x_1;
                    y_2 = y_1;
                    x_1 = x_1 + cos(theta(i));
                    y_1 = y_1 + sin(theta(i));
                    theta(i) = reflectPointOnC(x_0, y_0, x_1, y_1, x_2, y_2, circle_R);
                    index = 1;
                end
            end
            
            pos(i,1) = x_1;
            pos(i,2) = y_1;
        end
        
        if theta(i) < 0
            theta(i) = theta(i) + 2*pi;
        elseif theta(i) > 2*pi
            theta(i) = theta(i) - 2*pi;
        end
        theta_d(i) = theta(i);

    end
    
    theta_d = theta;
    theta_t_1 = theta;
    pos_t_1 = pos;

    pos_t(:,:,t) = pos; % store position 
    theta_t(:,t) = theta; % store orientation

end