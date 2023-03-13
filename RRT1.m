function [Wp] = RRT1(x,y,x_goal,y_goal,x_obs,y_obs,obstacle_dim,beh)

flag_plots  = 0;
x_max       = 15;
y_max       = 15;
obstacle    = zeros(length(x_obs),4);

for j = 1:length(x_obs)
obstacle(j,:)  = [x_obs(j),y_obs(j),obstacle_dim(j,1),obstacle_dim(j,2)];
end

EPS = 1;%0.5;
if length(str2num(beh))
    numNodes = 1;%1500;
else
    numNodes = 0;
end
q_start.coord = zeros(2,1);
q_start.coord(1) = x;
q_start.coord(2) = y;
q_start.cost = 0;
q_start.parent = 0;
q_goal.coord = zeros(2,1);
q_goal.coord(1) = x_goal;
q_goal.coord(2) = y_goal;
q_goal.cost = 0;

nodes(1) = q_start;
if flag_plots == 1
    for j = 1:length(x_obs)
        rectangle('Position',obstacle(j,:),'FaceColor',[.8 .8 .8],'EdgeColor',[.7 .7 .7])
    end
    axis equal
    hold on
end
h_waitbar = waitbar(0, 'Planning...');
for i = 1:1:numNodes
         q_rand = [(rand(1)*x_max) (rand(1)*y_max)];
    %     q_rand = [2*randunif(1)+x 2*randunif(1)+y];
    
    
    if flag_plots == 1
        plot(q_rand(1), q_rand(2), 'x', 'Color',  [0 0.4470 0.7410])
    end
    
    % Break if goal node is already reached
    for j = 1:1:length(nodes)
        if nodes(j).coord == q_goal.coord
            break
        end
    end
    
    % Pick the closest node from existing list to branch out from
    ndist = [];
    for j = 1:1:length(nodes)
        n = nodes(j);
        tmp = dist1(n.coord, q_rand);
        ndist = [ndist tmp];
    end
    [val, idx] = min(ndist);
    q_near = nodes(idx);
    
    q_new.coord = steer(q_rand, q_near.coord, val, EPS);
    if noCollision(q_rand, q_near.coord, obstacle)
        if flag_plots == 1
            line([q_near.coord(1), q_new.coord(1)], [q_near.coord(2), q_new.coord(2)], 'Color', 'k', 'LineWidth', 2);
            drawnow
            hold on
        end
        q_new.cost = dist1(q_new.coord, q_near.coord) + q_near.cost;
        
        % Within a radius of r, find all existing nodes
        q_nearest = [];
        r = 60;
        neighbor_count = 1;
        for j = 1:1:length(nodes)
            if noCollision(nodes(j).coord, q_new.coord, obstacle) && dist1(nodes(j).coord, q_new.coord) <= r
                q_nearest(neighbor_count).coord = nodes(j).coord;
                q_nearest(neighbor_count).cost = nodes(j).cost;
                neighbor_count = neighbor_count+1;
            end
        end
        
        % Initialize cost to currently known value
        q_min = q_near;
        C_min = q_new.cost;
        
        % Iterate through all nearest neighbors to find alternate lower
        % cost paths
        
        for k = 1:1:length(q_nearest)
            if noCollision(q_nearest(k).coord, q_new.coord, obstacle) && q_nearest(k).cost + dist1(q_nearest(k).coord, q_new.coord) < C_min
                q_min = q_nearest(k);
                C_min = q_nearest(k).cost + dist1(q_nearest(k).coord, q_new.coord);
                if flag_plots == 1
                    line([q_min.coord(1), q_new.coord(1)], [q_min.coord(2), q_new.coord(2)], 'Color', 'g');
                end
                %                 hold on
            end
        end
        
        % Update parent to least cost-from node
        for j = 1:1:length(nodes)
            if nodes(j).coord == q_min.coord
                q_new.parent = j;
            end
        end
        
        % Append to nodes
        nodes = [nodes q_new];
    end
    waitbar(i/(numNodes-1));
    
end

D = [];
for j = 1:1:length(nodes)
    tmpdist = dist1(nodes(j).coord, q_goal.coord);
    D = [D tmpdist];
end

% Search backwards from goal to start to find the optimal least cost path
[val, idx] = min(D);
q_final = nodes(idx);
q_goal.parent = idx;
q_end = q_goal;
nodes = [nodes q_goal];
i = 0 ;
while q_end.parent ~= 0
    start = q_end.parent;
%     if flag_plots == 1
        line([q_end.coord(1), nodes(start).coord(1)], [q_end.coord(2), nodes(start).coord(2)], 'Color', 'r', 'LineWidth', 2);
        hold on
%     end
    q_end = nodes(start);
    i = i+1;
    
end
%Wp = zeros(i,2);
q_end = q_goal;
nodes = [nodes q_goal];
k=0;
while q_end.parent ~= 0
    k = k+1;
    start = q_end.parent;
    %     Q = q_end
    q_end = nodes(start);
    Wp(k,:) = q_end.coord;
end

close(h_waitbar);



end

