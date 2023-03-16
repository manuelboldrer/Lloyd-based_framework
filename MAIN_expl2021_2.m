%/*-----------------------------------------------------
%|      Manuel Boldrer, PhD                            |
%|      Department of Cognitive Robotics               |
%|      Delft University of Technology                 |
%|                                                     |
%|      email: m.boldrer@tudelft.nl                    |
%|      https://manuelboldrer.github.io/               |                                                   
%-----------------------------------------------------*/
%The algorithms implemented in this code were used to generate
%part of the simulation results in the following papers:

%[1] Boldrer, M., Palopoli, L., & Fontanelli, D. (2022).
% A unified Lloyd-based framework for multi-agent collective behaviours. 
% Robotics and Autonomous Systems, 156, 104207.

%[2] Boldrer, M., Palopoli, L., & Fontanelli, D. (2020, October). 
% Lloyd-based Approach for Robots Navigation in Human-shared environments. 
% In 2020 IEEE/RSJ International Conference on Intelligent Robots and
% Systems.

%[3] Boldrer, M., Bevilacqua, P., Palopoli, L., & Fontanelli, D. (2021).
% Graph connectivity control of a mobile robot network with mixed dynamic multi-tasks. 
% IEEE Robotics and Automation Letters, 6(2), 1934-1941.

%[4] Boldrer, M., Pasqualetti, F., Palopoli, L., & Fontanelli, D. (2022). 
% Multiagent persistent monitoring via time-inverted kuramoto dynamics. 
% IEEE Control Systems Letters, 6, 2798-2803.


close all force
clearvars
clc
format compact



%% Parameters
beh             = ['1','1','1','2','2','2','c','c','c','e','e','e'];  %behavior of the agent: '1','2','3' = flocking,
% 'c' = coverage,'r' = rendezvous, %'e' = exploration. 
dx              = 0.05 ;%0.075  %space discretization
tend            = 10 ;          %ending time of the simulation
n_agents        = length(beh);  %number of robots
n_wp            = 1   ;         %number of waypoints
sizeA = 0.35*ones(n_agents,1);  %robot physical dimension (radius)

xi    = sizeA*0.35+0.1;
Rcoh(1,:)            = 1000*ones(n_agents,1);
Rcoht0          = Rcoh;
dt              = 0.05;%0.05;%033;%033   ;  %time step
time            = 0:dt:tend ;                %time vector

% Rmax            = 3*ones(n_agents,1);%
for i = 1:n_agents
    if beh(i) == '1' || beh(i) == '2' || beh(i) == '3' || beh(i) == '4' || beh(i) == 'r'
        Rmax(i) = 0.5;   %Maximum sensing radius
    elseif beh(i) == 'c' || beh(i) == 'e'
        Rmax(i) = 0.5;
    end
    Rmaxt0(i) = Rmax(i);
    Rmint0(i) = Rmax(i);
end
RmaxComm        = 5*ones(n_agents,1);  %Communication Radius  
xlim            = [0 30] ;
ylim            = [0 30] ;
bbox            = [xlim(1),ylim(1);xlim(1),ylim(2);xlim(2),ylim(2);xlim(2),ylim(1)];
obs             = [1] ;
n_obsD          = length(obs)  ;
sizeO           = .5; % Radius of dynamic obstalce
epsilon         = .1;    % safety factor
b               = 0.2;
%% Initialize variables
%Delete not used variables

x               = zeros(length(time),n_agents);
y               = zeros(length(time),n_agents);
theta           = zeros(length(time),n_agents);
v               = cell(length(time),n_agents);
c               = cell(length(time),n_agents);
cx              = zeros(length(time),n_agents);
cy              = zeros(length(time),n_agents);
cx2             = zeros(length(time),n_agents);
cy2             = zeros(length(time),n_agents);
XX1             = cell(length(time)-1,n_agents);
YY1             = cell(length(time)-1,n_agents);
XX2             = cell(length(time)-1,n_agents);
YY2             = cell(length(time)-1,n_agents);
x_obs_dyna      = zeros(length(time),n_obsD);
y_obs_dyna      = zeros(length(time),n_obsD);
VobsD           = zeros(length(time),n_obsD,2);
index_prev_vec  = ones(n_agents,1);
vel             = zeros(length(time)-1,n_agents);      %control v (linear velocity)
omega           = zeros(length(time)-1,n_agents);      %control omega (angular velocity)
x_dot           = zeros(length(time),n_agents);        %derivative of x
y_dot           = zeros(length(time),n_agents);        %derivative of y
theta_dot       = zeros(length(time),n_agents);        %derivative of theta
hh              = zeros(2,length(time),n_agents);      %heading of vehicle
hd              = zeros(2,length(time),n_agents);      %heading of vehicle
UU              = zeros(n_agents,2,length(time));
x_Levy          = zeros(n_agents,1);
y_Levy          = zeros(n_agents,1);
ux              = zeros(length(time),n_agents);
uy              = zeros(length(time),n_agents);
xdot            = zeros(length(time),n_agents);
ydot            = zeros(length(time),n_agents);
xddot           = zeros(length(time),n_agents);
yddot           = zeros(length(time),n_agents);
IN1             = cell(n_agents,1);
IN2             = cell(n_agents,1);
Uglob           = cell(length(time),1);
Xvis1           = cell(length(time),n_agents);
Yvis1           = cell(length(time),n_agents);
goal            = zeros(n_agents,2,n_wp);
Wp              = cell(n_agents,1);
vwp             = zeros(n_agents,2,length(time));      %desired velocity
wpp             = ones(n_agents,1);
wp_input_vec    = zeros(2, length(time), n_agents);
wp_input_vec1   = zeros(2, length(time), n_agents);
wp_path         = cell(2,1);
index_wp_smart  = zeros(length(time), n_agents);
tmp1            = zeros(n_agents,1);
Xvis            = cell(length(time),n_agents);
Yvis            = cell(length(time),n_agents);
Xvisp           = cell(length(time),n_agents);
Yvisp           = cell(length(time),n_agents);
count1          = zeros(length(time),n_agents);
griglia         = zeros(length(time),length(Xvis));
behtime         = zeros(length(time),n_agents);
Add             = zeros(length(time),n_agents,n_agents);
JJJ             = zeros(length(time),n_agents);
Jo              = cell(length(time),1);
G               = zeros(n_agents);
% countL = zeros(n_agents);
xbar            = zeros(n_agents,1);
ybar            = zeros(n_agents,1);
sizeComm        = zeros(n_agents,n_agents);
in1             = cell(n_agents,1);
JJ              = cell(length(time),1);
%% Flags

flag_expl        = 1; 
flag_cm          = 1; %connectivity maintenance
nonholo          = 0; %nonholonomic constraints
task_switch_flag = 0; %allow dynamic task switching
vwmst_flag       = 1; %visible weighted minimum spanning tree 

flag_obs        = 18; %select obstacle configuration
video_flag      =  0; 
manual_ics      =  0; %reset initial conditions
planning_flag   =  0; %global path planner
manual_goal     =  0; %reset the final goal position 

%% Obstacle definition
switch flag_obs
    
    case 19
        
        xlim          = [0 5.5] ;
        ylim          = [0 5.5] ;
        x_obs         = [0,0,5,0,2];%,0,3];%,0.5];
        y_obs         = [0,0,0,5,2.7];%,0,0];%,2];
        obstacle_dim  = [5,0.6;0.6,5;0.6,5;5.6,0.6;1,.5];%;2,2;2,2];%;0.8 0.5] ;
    case 18
        
        xlim          = [0 5.5] ;
        ylim          = [0 5.5] ;
        x_obs         = [0,0,5,0,2];%,0,3];%,0.5];
        y_obs         = [0,0,0,5,2.7];%,0,0];%,2];
        obstacle_dim  = [5,0.6;0.6,5;0.6,5;5.6,0.6;1,.5];%;2,2;2,2];%;0.8 0.5] ;
        
    case 17
        x_obs         = [0,0,15,0,5,10,5,10,2,3,12,13];
        y_obs         = [0,0,0,7.5,5,5,0,0,0,5,2,5.5];
        obstacle_dim  = [15,0.2;0.2,7.5;0.2,7.5;15.2,0.2;0.2,2.5;0.2,2.5;0.2,4;0.2,4;...
            1,3;2,1;1,2;2,0.7];
        xlim          = [0 16] ;
        ylim          = [0 8] ;
    case 16
        x_obs         = [0,0,15,0,5,5,0,5,10,10,10,10,2,2,12,12,12,12,2,2];
        y_obs         = [0,0,0,15,0,6.5,7.5,13,0,6.5,7.5,13,2,5,2,5,10,12,10,11];
        obstacle_dim  = [15,0.6;0.6,15;0.6,15;15,0.6;0.4,4.5;...
            0.4,4.5;5,0.4;0.4,2;0.4,4.5;0.4,4.5;5,0.4;0.4,2;1,1;1,1;1,1;1,1,;1,2;1,1,;2,1;1,1];
        xlim          = [0 30] ;
        ylim          = [0 30] ;
    case 15
        x_obs         = [0,0,15,0,5,5,0,5,10,10,10,10];
        y_obs         = [0,0,0,15,0,6.5,7.5,13,0,6.5,7.5,13];
        obstacle_dim  = [15,0.6;0.6,15;0.6,15;15,0.6;0.4,4.5;...
            0.4,4.5;5,0.4;0.4,2;0.4,4.5;0.4,4.5;5,0.4;0.4,2];
        xlim          = [0 30] ;
        ylim          = [0 30] ;
        
    case 14
        x_obs         = [0,0,30,0,10,10,0,10,20,20,20,20];
        y_obs         = [0,0,0,30,0,13,15,25,0,13,15,25];
        obstacle_dim  = [30,0.6;0.6,30;0.6,30;30,0.6;0.4,10;...
            0.4,10;10,0.4;0.4,6;0.4,10;0.4,10;10,0.4;0.4,6];
        xlim          = [0 30] ;
        ylim          = [0 30] ;
        
    case 13
        xlim           = [0 15];
        ylim           = [0 15];
        %x_obs         = [0,0,15,0,0,10,5,10,5] ;
        %y_obs         = [0,0,0,15,5,0,5,10,5]  ;
        x_obs          = [0,0,15,0,10,5,4,10];%,7.5];%,7,7] ;
        y_obs          = [0,0,0,15,0,5,5,10];%,7.5];%,7,8]  ;
        obstacle_dim   = [15,0.6;0.6,15;0.6,15;15,0.6;0.5,5;...
            0.5,5;2.5,0.5;5,0.5];%;1,1];
    case 12
        xlim         = [0 15];
        ylim         = [0 15];
        x_obs        = [0,0,15,0,7];
        y_obs        = [0,0,0,15,5];
        obstacle_dim = [15,0.6;0.6,15;0.6,15;15,0.6;1,5] ;
        
    case 11
        xlim         = [0 15];
        ylim         = [0 15];
        x_obs        = [0,0,15,0,2.5,12,5,12,7,9.5,9,5.5];
        y_obs        = [0,0,0,15,5,7,10,4,4.5,10.5,7.5,7.25];
        obstacle_dim = [15,0.6;0.6,15;0.6,15;15,0.6;1.5,2;0.6,1.2;1.3,1;1,1;1,1.5;3,1;1,1;1.5,1] ;
        
    case 10
        xlim         = [0 15];
        ylim         = [0 15];
        x_obs        = [0,0,15,0,2.5,10,5,10,7];
        y_obs        = [0,0,0,15,5,7,10,4,7];
        obstacle_dim = [15,0.6;0.6,15;0.6,15;15,0.6;1.5,2;0.6,1.2;1.3,1;1,1;1.3,0.8] ;
    case 9
        xlim         = [0 4];
        ylim         = [0 4];
        x_obs        = [0,0,4,0,2,2];
        y_obs        = [0,0,0,4,0,2];
        obstacle_dim = [5,0.5;0.1,5;0.6,5;5,0.6;2,1.5;1,2] ;
        
    case 8  
        xlim         = [0 15];
        ylim         = [0 15];
        x_obs        = [0,0,15,0,2.5,10,5,10];
        y_obs        = [0,0,0,15,5,7,10,4];
        obstacle_dim = [15,0.6;0.6,15;0.6,15;15,0.6;5,3;2,5;3,3;1,1] ;
        
    case 7
        xlim         = [0 15];
        ylim         = [0 15];
        %x_obs       = [0,0,15,0,0,10,5,10,5] ;
        %y_obs       = [0,0,0,15,5,0,5,10,5]  ;
        x_obs        = [0,0,15,0,0,10,5,5,10];%,7.5];%,7,7] ;
        y_obs        = [0,0,0,15,5,0,5,5,10];%,7.5];%,7,8]  ;
        obstacle_dim = [15,0.6;0.6,15;0.6,15;15,0.6;5,0.5;0.5,5;...
            0.5,5;2.5,0.5;5,0.5];%;1,1];
    case 0
        x_obs         = [0,0,4,0,1.7,4];
        y_obs         = [0,0,0,4,1.7,4];
        obstacle_dim  = [4,0.5;0.5,4;0.5,4;4,0.5;1,1;0.5,0.5];
        xlim          = [0 4] ;
        ylim          = [0 4] ;
    case -1
        x_obs         = [-5,-5,5,-5];
        y_obs         = [-5,-5,-5,5];
        obstacle_dim  = [10,0.6;0.6,10;0.6,10;10,0.6];
        xlim          = [-5 5] ;
        ylim          = [-5 5] ;
    case 1
        x_obs         = [10,20,0,6];
        y_obs         = [10,15,15,15];
        obstacle_dim  = [10,10;10,5;4,1;4,1];
        xlim          = [0 30] ;
        ylim          = [0 30] ;
    case 2
        x_obs         = [0,0,30,0,0.1,10,10];
        y_obs         = [0,0,0,30,20,10,0];
        obstacle_dim  = [30,0.6;0.6,30;0.6,30;30,0.6;20,5;20,5;3,5];
        xlim          = [0 30] ;
        ylim          = [0 30] ;
    case 3
        x_obs         = [7.5]; %#ok<*NBRAK>
        y_obs         = [7.5];
        obstacle_dim  = [15,15];
        xlim          = [0 30] ;
        ylim          = [0 30] ;
    case 4
        xlim          = [5 25] ;
        ylim          = [5 25] ;
        x_obs         = [10,12.5,18.75,10,12.5]; %#ok<*NBRAK>
        y_obs         = [10,10,11.25,17.5,18.75];
        obstacle_dim  = [2.5,2.5;7.5,1.25;1.25,8.75;2.5,2.5;6.25,1.25];
        bbox          = [xlim(1),ylim(1);xlim(1),ylim(2);xlim(2),ylim(2);xlim(2),ylim(1)];
    case 5
        x_obs         = [0,0,30,0,5,5,17,17];
        y_obs         = [0,0,0,30,5,17,5,17];
        obstacle_dim  = [30,0.6;0.6,30;0.6,30;30,0.6;8,8;8,8;8,8;8,8];
        xlim          = [0 30] ;
        ylim          = [0 30] ;
    case 6
        x_obs         = [0,10,0,10];
        y_obs         = [10,10,19.5,16.5];
        obstacle_dim  = [10,0.5;9,3.5;10,0.5;9,3.5];
        xlim          = [-5 5] ;
        ylim          = [-5 5] ;
end

obstacle  = zeros(length(x_obs),4);
obstacleD = nan(length(time),n_obsD,4);

for j = 1:length(x_obs)
    obstacle(j,:)  = [x_obs(j),y_obs(j),obstacle_dim(j,1),obstacle_dim(j,2)];
end
%% Define Initial Conditions
%initial robots' state (x,y,theta)
sec_numb = 1;
if manual_ics == 1
    fprintf('%d) Define initial conditions of the vehicles \n',sec_numb); sec_numb = sec_numb + 1;
    figure('Name', 'Environment','units','normalized','outerposition',[0 0 1 1]);
    hold on;
    title('Map','Interpreter','latex');
    grid on;
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    axis equal;
    grid minor
    axis([xlim,ylim])
    obstacles = zeros(length(x_obs),4);
    for j = 1:length(x_obs)
        obstacles(j,:)  = [x_obs(j),y_obs(j),obstacle_dim(j,1),obstacle_dim(j,2)];
    end
    for j = 1:length(x_obs)
        rectangle('Position',obstacles(j,:),'FaceColor',[.8 .8 .8],'EdgeColor',[.7 .7 .7]);
    end
    
    for i = 1:n_agents
        display(['insert vehicle position of agent ', num2str(i)]);
        tmp    = ginput(1);
        x(1,i) = tmp(1);
        y(1,i) = tmp(2);
        vehicle_ref_point = plot(tmp(1), tmp(2), '.', 'markersize', 8, 'color', 'k');
        display(['insert second point for the heading of vehicle ', num2str(i)]);
        tmp2       = ginput(1);
        theta(1,i) = atan2(tmp2(2) - tmp(2), tmp2(1) - tmp(1));
        delete(vehicle_ref_point);
        plot_unicycle(x(1,i), y(1,i), theta(1,i), 'k',sizeA(i));
    end
    x_1     = x(1,:);
    y_1     = y(1,:);
    theta_1 = theta(1,:);
    save('x_1','x_1');
    save('y_1','y_1');
    save('theta_1','theta_1');
else
    load('x_1');
    load('y_1');
    load('theta_1');
    
    x(1,:)     = x_1;
    y(1,:)     = y_1;
    theta(1,:) = theta_1;
end

%goal position (for the case of coverage, exploration it is not used TO FIX)
if manual_goal == 1
    
    for i = 1:n_agents
        for w = 1:n_wp
            display(['insert goal for agent ', num2str(i)]);
            tmp         = ginput(1);
            goal(i,1,w) = tmp(1);
            goal(i,2,w) = tmp(2);
            plot(goal(i,1,w), goal(i,2,w) ,'x');
            
        end
        
    end
    save('goal','goal');
else
    load('goal');
end

%% Meshgrid

x_grid = xlim(1):dx:xlim(2);
y_grid = ylim(1):dx:ylim(2);

[X,Y]  = meshgrid(x_grid,y_grid);
XX     = reshape(X,[size(X,1)*size(X,2),1]);
YY     = reshape(Y,[size(Y,1)*size(Y,2),1]);
%% Definition of the Probability density function \phi(q) (1) in [1]
%[1] Manuel Boldrer, Luigi Palopoli, Daniele Fontanelli,
%A unified Lloyd-based framework for multi-agent collective behaviours,
%Robotics and Autonomous Systems.

U_0    = 10000;
Uw = cell(n_agents,1);
for j = 1:n_agents
    Uw{j} = U_0*ones(length(XX),1);
end
R      = zeros(length(time),n_agents);
syms r_x r_y wpx wpy RR %hx hy
r     = [r_x; r_y];
wp    = [wpx; wpy];
U     = U_0 * exp(-(norm(r-wp))/RR);
U     = matlabFunction(U, 'File','UUfun');

%% Static obstacle inflation
obstacle2 = cell(n_agents,1);

obstacle1   = [obstacle];
for j =1:n_agents
    epsi2(j)       = sizeA(j)*0.35/2+dx;
    obstacle2{j}   = [obstacle1(:,1)-epsi2(j),obstacle1(:,2)-epsi2(j),obstacle1(:,3)+2*epsi2(j),obstacle1(:,4)+2*epsi2(j)];
end


%% Iteratively apply Lloyd's Algorithm
warning('off')
h_waitbar = waitbar(0, 'Simulating...');
%Dynamic obstacles
x_obs_dyna(1,1)          = 0    ;
y_obs_dyna(1,1)          = 0    ;

        
        for kk = 1:length(time)-1
            for j = 1:n_agents
                if planning_flag == 1 %&& length(str2num(beh(j)))>0 %global path planner it is necessary for behaviors 1,2,3 
                    if kk == 1  % is not necessary to compute every step, just for all agents
                        Wp{j}    = RRT1(x(1,j),y(1,j),goal(j,1,1),goal(j,2,1),x_obs-xi(j),y_obs-xi(j),obstacle_dim+2.2*xi(j),beh(j));
                        Wp{j}    = [goal(j,:,1);Wp{j}];
                        Wp{j}    = flip(Wp{j}); %waypoints computed from the RRTstar
                        save('Wp','Wp');
                    end
                else
                    if kk == 1
                        load('Wp');
                    end
                end
                [wp_input_vec(:,kk,j), index_prev_vec(j), wp_path{j}, index_wp_smart(kk,j),~] = generate_wp_path(Wp{j}', x(kk,j), y(kk,j), index_prev_vec(j),kk,obstacle2{j},Rmaxt0(j));

            end

 %% Dynamic obstacles velocity
            VobsD(kk,1,1)    = -0*(x_obs_dyna(kk,1)-x(kk,j))/norm([x_obs_dyna(kk,1)-x(kk,j),y_obs_dyna(kk,1)-y(kk,j)]);%1  ;%3*rand()     ;
            VobsD(kk,1,2)    = -0*(y_obs_dyna(kk,1)-y(kk,j))/norm([x_obs_dyna(kk,1)-x(kk,j),y_obs_dyna(kk,1)-y(kk,j)]);%-.0 ;

            
            for q = 1:n_obsD
                x_obs_dyna(kk+1,q) =  x_obs_dyna(kk,q) + dt*VobsD(kk,q,1);
                y_obs_dyna(kk+1,q) =  y_obs_dyna(kk,q) + dt*VobsD(kk,q,2);

                %Add uncertainty on the obstacle velocity reading 
%                 VobsD(kk,q,1)      =  VobsD(kk,q,1) + VobsD(kk,q,1)*0.3*randn();
%                 VobsD(kk,q,2)      =  VobsD(kk,q,2) + VobsD(kk,q,2)*0.3*randn();
            end
            
            for jj = 1:length(x_obs)
                obstacle(jj,:)  = [x_obs(jj),y_obs(jj),obstacle_dim(jj,1),obstacle_dim(jj,2)];
            end
            
            for qq = 1:n_obsD
                obstacleD(kk,qq,:) =  [x_obs_dyna(kk,qq),y_obs_dyna(kk,qq),sizeO,sizeO];
            end
            neigh = cell(n_agents,1);
            for j=1:n_agents     %Generation of the Voronoi cells           
                [v{kk,j},c{kk,j}] = VoronoiBounded([x(kk,:)';x_obs_dyna(kk,:)';],[y(kk,:)';y_obs_dyna(kk,:)'], bbox);
            end
            Ad = zeros(n_agents); %Adjacency matrix
            for m = 1 :length(Ad)
                for n = 1:length(Ad)
                    if  norm([x(kk,m),y(kk,m)]-[x(kk,n),y(kk,n)]) > min(RmaxComm(m),RmaxComm(n)) %limits on the distance
                        Ad(m,n) = 0;
                        Ad(n,m) = 0;
                    else
                        Ad(m,n) = 1;
                        Ad(m,n) = 1;
                    end
                    if noCollision([x(kk,m),y(kk,m)]',[x(kk,n),y(kk,n)]',obstacle2{j}) ==0  %limits due to the visibility
                        Ad(m,n) = 0;
                        Ad(n,m) = 0;
                    end
                end
            end
            
            
    

            J = zeros(n_agents);
            for m = 1:n_agents
                for n = 1:n_agents
                    if Ad(m,n) > 0
%                         if beh(n) == "1"
%                             J(m,n) = norm([x(kk,m)-wp_input_vec(1,kk,n),y(kk,m)-wp_input_vec(2,kk,n)]);
%                         elseif beh(n) == "2"
%                             J(m,n) = norm([x(kk,m)-wp_input_vec(1,kk,n),y(kk,m)-wp_input_vec(2,kk,n)]);
%                         elseif beh(n) == "3"
%                             J(m,n) = norm([x(kk,m)-wp_input_vec(1,kk,n),y(kk,m)-wp_input_vec(2,kk,n)]);
%                         elseif beh(n) == "4"
%                             J(m,n) = norm([x(kk,m)-wp_input_vec(1,kk,n),y(kk,m)-wp_input_vec(2,kk,n)]);
%                         elseif kk > 1 && beh(n) == "c"
%                             J(m,n) = 1/length(XX1{kk-1,m});
%                         elseif beh(n) == "r"
%                             J(m,n) = sqrt(var([x(kk,m);nonzeros(x(kk,:).*Adr(m,:))])^2+var([y(kk,m);nonzeros(y(kk,:).*Adr(m,:))])^2);
%                         elseif kk > 1 && beh(n) == "e"
%                             J(m,n) = 1/length(XX1{kk-1,m});
%                         end

                        if kk>1
                            J(m,n) = Jo{kk-1}(m,n);
                        end
                    end

                end
            end
            
            
            jcount = zeros(n_agents);
            jcount1 = zeros(n_agents);
            
            if task_switch_flag == 1
                for m = 1:n_agents
                    for n = 1:n_agents
                        %add visiblity constraint
                        if  length(str2num(beh(m)))>0 && length(str2num(beh(n)))>0  && Ad(m,n)>0 &&...
                                J(m,m)>=J(m,n) && J(n,n)>=J(n,m)...
                                && jcount(m,n)==0
                            tmpbeh = beh(n);
                            beh(n) = beh(m);
                            beh(m) = tmpbeh;
                            tmp0m=Wp{m} ;
                            Wp{m} = Wp{n};
                            Wp{n} = tmp0m;
                            jcount(m,n) = 1;
                            jcount(n,m) = 1;
                            tmpind = index_prev_vec(m) ;
                            index_prev_vec(m) = index_prev_vec(n);
                            index_prev_vec(n) = tmpind ;
                            tmpgoal = goal(m,:,1);
                            goal(m,:,1) = goal(n,:,1);
                            goal(n,:,1) = tmpgoal;
                            tmpwp = wp_input_vec(:,kk,m);
                            wp_input_vec(:,kk,m) = wp_input_vec(:,kk,n);
                            wp_input_vec(:,kk,n) = tmpwp;
                        end
                        
                        
                        %             if kk>1 && beh(m) == "c" && length(str2num(beh(n)))>0 ...
                        %                     && beh(m)~=beh(n) && Ad(m,n)>0 &&...
                        %                     J(m,m) >= J(n,m) && J(n,n) >= J(m,n)...
                        %                     && jcount1(m,n)==0
                        %                 tmpbeh = beh(n);
                        %                 beh(n) = beh(m);
                        %                 beh(m) = tmpbeh;
                        %                 tmp0m=Wp{m} ;
                        %                 Wp{m} = Wp{n};
                        %                 Wp{n} = tmp0m;
                        %                 jcount1(m,n) = 1;
                        %                 jcount1(n,m) = 1;
                        %                 tmpind = index_prev_vec(m) ;
                        %                 index_prev_vec(m) = index_prev_vec(n);
                        %                 index_prev_vec(n) = tmpind ;
                        %                 tmpgoal = goal(m,:,1);
                        %                 goal(m,:,1) = goal(n,:,1);
                        %                 goal(n,:,1) = tmpgoal;
                        %             end
                        
                        %             if kk>1 && (( beh(m) == "e" && beh(n) == "1") || (beh(n) == "e" && beh(m) =="2" ))...
                        %                     && beh(m)~=beh(n) && Ad(m,n)>0 &&...
                        %                     J(m,m) >= J(n,m) && J(n,n) >= J(m,n)...
                        %                     && jcount1(m,n)==0
                        if kk>1 &&  beh(m) == "e" && length(str2num(beh(n)))>0 ...
                                && beh(m)~=beh(n) && Ad(m,n)>0 &&...
                                J(m,m) >= J(n,m)  && J(n,n) >= J(m,n)...
                                && jcount1(m,n)==0
                            tmpbeh = beh(n);
                            beh(n) = beh(m);
                            beh(m) = tmpbeh;
                            tmp0m=Wp{m} ;
                            Wp{m} = Wp{n};
                            Wp{n} = tmp0m;
                            jcount1(m,n) = 1;
                            jcount1(n,m) = 1;
                            tmpind = index_prev_vec(m) ;
                            index_prev_vec(m) = index_prev_vec(n);
                            index_prev_vec(n) = tmpind ;
                            tmpgoal = goal(m,:,1);
                            goal(m,:,1) = goal(n,:,1);
                            goal(n,:,1) = tmpgoal;
                        end
                        
                        %             if kk>1 && (beh(m) == "c" || beh(n) == "c") && beh(m)~=beh(n) && Ad(m,n)>0 &&...
                        %                     J(m,m) >= J(n,m) && J(n,n) >= J(m,n)...
                        %                     && jcount1(m,n)==0
                        %                 tmpbeh = beh(n);
                        %                 beh(n) = beh(m);
                        %                 beh(m) = tmpbeh;
                        %                 tmp0m=Wp{m} ;
                        %                 Wp{m} = Wp{n};
                        %                 Wp{n} = tmp0m;
                        %                 jcount1(m,n) = 1;
                        %                 jcount1(n,m) = 1;
                        %                 tmpind = index_prev_vec(m) ;
                        %                 index_prev_vec(m) = index_prev_vec(n);
                        %                 index_prev_vec(n) = tmpind ;
                        %             end
                    end
                end
            end
            
            %     if kk ==17
            %         keyboard
            %     end
            
            
            %
            Adtmp0 = Ad;
            
            if vwmst_flag ==1
                for m =   1:length(Ad)
                    for n = 1:length(Ad)
                        if beh(m) == beh(n)
                            if m == n
                                Ad(m,n) = 1;
                            else
                                if kk ==1
                                else
                                    Ad(m,n) = Ad(m,n) * (sizeComm(m,n));
                                end
                            end
                        elseif beh(m)== 'c' || beh(n) =='c' ||beh(m)== 'e' || beh(n) =='e'
                            Ad(m,n) = 100*Ad(m,n) * (sizeComm(m,n));
                        else
                            Ad(m,n) = 1000*Ad(m,n) * (sizeComm(m,n));
                            
                        end
                    end
                end
                Ad = Ad+Ad';
                
                
                
                G     = graph(Ad);
                [MST] = minspantree(G);
                sAd   =  adjacency(MST);
                Ad    = full(sAd);
            end
            %    Ad
            %     for m = 1:length(Ad)
            %         for n = 1:length(Ad)
            %             Adtmp = Ad;
            %             Adtmp(m,n) = 0;
            %             Adtmp(n,m) = 0;
            %             eigtmp = eig(diag(sum(Adtmp))-Adtmp);
            %             if kk>1 && m ~= n && eigtmp(2) > 1e-6 &&... %sum(Ad(m,:))==3 &&...
            %                     Ad(m,n) == 1 &&  (edgeselected(m) ~=n && edgeselected(n) ~=m)%isempty(find(edgeselected == n))
            %
            %                     % Ad(m,n) == 1 &&  (edgeselected(m) ~=n && edgeselected(n) ~=m ) %isempty(find(edgeselected == n))
            %
            %                 Ad(m,n) = 0;
            %                 Ad(n,m) = 0;
            %             end
            %         end
            %     end
            
            
            
            
            
            %     for m =   1:length(Ad)
            %         for n = 1:length(Ad)
            %             Adtmp = Ad;
            %             Adtmp(m,n) = 0;
            %             Adtmp(n,m) = 0;
            %             eigtmp = eig(diag(sum(Adtmp))-Adtmp);
            %             if kk>1 && m ~= n && sum(Ad(m,:))>3 && eigtmp(2) > 1e-6 &&...
            %                     Ad(m,n) == 1 &&  edgediscarded(m) ==n %isempty(find(edgeselected == n))
            %                 Ad(m,n) = 0;
            %                 Ad(n,m) = 0;
            %             end
            %         end
            %     end
            %         for m = 1 :length(Ad)
            %             for n = 1:length(Ad)
            %                 Adtmp = Ad;
            %                   eigtmp1 = eig(diag(sum(Ad))-Ad);
            %                 Adtmp(m,n) = 0;
            %                 Adtmp(n,m) = 0;
            %                 eigtmp = eig(diag(sum(Adtmp))-Adtmp);
            %     %
            %                  if   m ~= n && Ad(m,n) >0
            %                     if kk>1 && eigtmp(2) > 0.001 && eigtmp(2) < eigtmp1(2)
            %                         Ad(m,n) = 0;
            %                         Ad(n,m) = 0;
            %     %                     countL(m,n) = 1;
            %                     end
            %
            %                 end
            %             end
            %         end
            
            %Ad = [1 1 0 0 0;1 1 1 0 0;0 1 1 1 0;0 0 1 1 1;0 0 0 1 1];
            % Ad = [1 0 1 0 0;0 1 0 1 1;1 0 1 0 1;0 1 0 1 0;0 1 1 0 1];
            % Ad = zeros(n_agents);
            if kk == 1
                Ad0 = Ad;
            end
            
            
            %Ad = [0 1 0 0 0;1 0 1 0 0; 0 1 0 1 0 ; 0 0 1 0 1; 0 0 0 1 0];
            Add(kk,:,:) = Ad;
            %         Ad = neighbours(c{kk},n_agents);
            %     AAA = rand(n_agents)+eye(n_agents)>0.5;
            %     ASD =AAA*AAA'>=1;
            %     Ad = ASD*eye(n_agents);
            
            %     Ad1 = neighbours(c{kk},n_agents);
            %     for m = 1:n_agents
            %         for n = 1:n_agents
            %             if Ad1(m,n) >0 && beh(m)=='r' && beh(n) == 'r'
            %                 Ad1(m,n)=1;
            %             else
            %                 Ad1(m,n)=0;
            %             end
            %         end
            %     end
            %     De   = diag(sum(Ad1)) ; %-eye(length(Ad));
            %
            %     Lap  = De-Ad1 ;
            
            %
            % De   = diag(sum(Adtmp0)) ; %-eye(length(Ad));
            % Lap  = De-Adtmp0 ;
            %     Lap = diag(sum(Ad0))-Ad0;
            ww   = cell(n_agents,1);
            Xcom = cell(n_agents,1);
            Ycom = cell(n_agents,1);
            
            %% Compute Control
            for j = 1 : n_agents %calculate the centroid of each cell
                
                behtime(kk,j) =beh(j);
                
                if norm([mean(x(kk,:))-x(kk,j),mean(y(kk,:))-y(kk,j)]) > sqrt(Rcoht0(j))-0.25
                    Rcoh(j) = (norm([mean(x(kk,:))-x(kk,j),mean(y(kk,:))-y(kk,j)]))^2+0.25;
                else
                    Rcoh(j) = Rcoht0(j);
                end
                                
                obstacle1                 = [obstacle];
                obstacle2{j}              = [obstacle1(:,1)-epsi2(j),obstacle1(:,2)-epsi2(j),obstacle1(:,3)+2*epsi2(j),obstacle1(:,4)+2*epsi2(j)];
                distj                     = sqrt((x(kk,:)-x(kk,j)).^2+(y(kk,:)-y(kk,j)).^2);

                if kk ==1
                    [Xvis{kk,j},Yvis{kk,j}]   = visibilitypoints2(x(kk,:),y(kk,:),squeeze(atan2(hd(2,kk,:),hd(1,kk,:))),obstacle2{j},obstacleD(kk,:,:), VobsD(kk,:,:),1.5,xlim,ylim,Rmax,epsi2(j),j,sizeA.*0.35,Rcoh(j),theta(kk,:),b,Ad,Rmaxt0);
                else
                    [Xvis{kk,j},Yvis{kk,j}]   = visibilitypoints2(x(kk,:),y(kk,:),squeeze(atan2(hd(2,kk-1,:),hd(1,kk-1,:))),obstacle2{j},obstacleD(kk,:,:), VobsD(kk,:,:),1.5,xlim,ylim,Rmax,epsi2(j),j,sizeA.*0.35,Rcoh(j),theta(kk,:),b,Ad,Rmaxt0);
                end
                sizexvis = size(Xvis{kk,j},2);
                %         [Xvis{kk,j},Yvis{kk,j}]   = visibilitypoints2(x(kk,:),y(kk,:),obstacle2{j},obstacleD(kk,:,:), VobsD(kk,:,:),1.5,xlim,ylim,Rmax,epsi2(j),j,[1 1 1],Rcoh(j),theta(kk,:),b,Ad,Rmaxt0);                
                %         [Xvis{kk,j},Yvis{kk,j}]             = visibilitypointsNOWALL(x(kk,:),y(kk,:),obstacle2{j},obstacleD(kk,:,:), VobsD(kk,:,:),1.5,xlim,ylim,Rmax,epsi2(j),j,sizeA.*0.35,Rcoh(j),theta(kk,:),b,Ad,Rmaxt0);
                [Xvis2,Yvis2,count1(kk,j)] = visibilitypointsNOCOH(x(kk,:),y(kk,:),obstacle2{j},obstacleD(kk,:,:), VobsD(kk,:,:),1.5,xlim,ylim,Rmax,epsi2(j),j,sizeA.*0.35,Rcoh(j),theta(kk,:),b,Ad,Rmaxt0);             
                [Xvis3,Yvis3]              = visibilitypointsNONEIGH(x(kk,:),y(kk,:),obstacle2{j},obstacleD(kk,:,:), VobsD(kk,:,:),1.5,xlim,ylim,Rmax,epsi2(j),j,sizeA.*0.35,Rcoh(j),theta(kk,:),b,Ad,Rmaxt0);
                if flag_cm == 1
                    for n = 1:n_agents
                        
                        if j~=n && Adtmp0(j,n) > 0
                            [Xcom{n},Ycom{n}] = communication2(x(kk,:),y(kk,:),obstacle2{j},obstacleD(kk,:,:), VobsD(kk,:,:),1.5,xlim,ylim,Rmax(n)+Rmax(j),epsi2(j),n,sizeA(j)*0.35,Rcoh(j),theta(kk,:),b,Ad,RmaxComm(j));%Rmaxt0(j)+Rmaxt0(n));
                            [Xvis3,Yvis3]     = visibilitypoints2(x(kk,:),y(kk,:),squeeze(atan2(hd(2,kk,:),hd(1,kk,:))),obstacle2{j},obstacleD(kk,:,:), VobsD(kk,:,:),1.5,xlim,ylim,Rmax,epsi2(j),j,sizeA.*0.35,Rcoh(j),theta(kk,:),b,Ad,Rmaxt0);
                            
                            
                            ww{n}         = boundary(Xcom{n}',Ycom{n}',1);
                            [in1{n},~]    = inpolygon(Xvis{kk,j},Yvis{kk,j},Xcom{n}(ww{n}),Ycom{n}(ww{n}));
                            xvise = Xvis{kk,j}(in1{n});
                            yvise = Yvis{kk,j}(in1{n});
                            
                            if Ad(j,n) >0
                                Xvis{kk,j}    = Xvis{kk,j}(in1{n});
                                Yvis{kk,j}    = Yvis{kk,j}(in1{n});
                            end
                            sizeComm(j,n) = size(in1{n},2)/sum(in1{n});%sizexvis-size(xvise,2);%size(Xvis{kk,j},2)-sum(in1);
                            % sum( UUfun(R(kk,j),XX1{kk,j},YY1{kk,j},xbar(j),ybar(j)))     
                        end
                    end    
                end                
                X1          = Xvis{kk,j};
                Y1          = Yvis{kk,j};              
                if ~isempty(X1)
                    k           = boundary(X1',Y1',1);
                    k1          = boundary(Xvis3',Yvis3',1);
                else
                    k = [];
                    k1 = [];
                    
                end
                if ~isempty(k)
                    
                    [in11,on]    = inpolygon(XX,YY,X1(k),Y1(k));
                    XX1{kk,j}   = XX(in11);
                    YY1{kk,j}   = YY(in11);
                    IN1{j} = in11;
                    [in2,~]     = inpolygon(XX,YY,Xvis3(k1),Yvis3(k1));
                    XX2{kk,j}   = XX(in2);
                    YY2{kk,j}   = YY(in2);
                    IN2{j} = in2;
                else
                    XX1{kk,j} = x(kk,j);
                    YY1{kk,j} = y(kk,j);                  
                    XX2{kk,j} = x(kk,j);
                    YY2{kk,j} = y(kk,j);
                end
                

                %         if  (length(Xvis1)-length(Xviss))>50
                %             Rmax(j) = Rmax(j) + dt*(length(Xvis1)-length(Xviss))/10;
                %         elseif  ((length(Xvis1)-length(Xviss))<50 && (length(Xvis1)-length(Xviss))>10)
                %             Rmax(j) = Rmax(j) - dt*0*(length(Xvis1)-length(Xviss))/10;
                %         elseif (length(Xvis1)-length(Xviss))<10
                %             Rmax(j) = Rmax(j) - dt*100/10;
                %         end
                %         Rmax(j) = max(min(Rmax(j),3),Rmaxt0(j));
                %
                %         if  (length(Xvis2)-length(Xvis{kk,j}))>100 && length(Xvis{kk,j}) < 400
                
                %         Rcoh(j) = max(min(Rcoh(j),3*Rcoh(j)),Rcoht0(j));
                
                
                fx = cell(1);
                fy = cell(1);
                if beh(j) == "r"
                    ccount=0;
                    nei = find(Ad(j,:)>0);
                    for p=1:length(nei)
                        if beh(nei(p)) == "r"
                            ccount = ccount+1;
                            fx{1}(ccount)= x(kk,nei(p));
                            fy{1}(ccount)= y(kk,nei(p));
                        end
                    end
                    if beh(j) == "c"
                        fx{1}=x(kk,j);
                        fy{1}=y(kk,j);
                    end
                    
                    xbar(j) = mean(fx{1});
                    ybar(j) = mean(fy{1});
                elseif length(str2num(beh(j)))>0
                    xbar(j) = wp_input_vec(1,kk,j);
                    ybar(j) = wp_input_vec(2,kk,j);
                end
            end
         
            %     for j =1 :n_agents
            %         if count1(kk,j) < 1 && beh(j) == "c"
            %             Rmax(j) = Rmax(j) - 1*(Rmax(j)-Rmaxt0(j))*dt;
            %         else
            %             Rmax(j) = Rmax(j) - 1*(Rmax(j)-Rmint0(j))*dt;%-Rmaxt0(j))*dt;
            %         end
            %
            %     end
  
            for j = 1:n_agents
                if flag_expl == 1 || beh(j) == "e"
                    k = boundary(XX1{kk,j},YY1{kk,j},1);
                    for q = 1:length(XX)
                        % if inpolygon(XX(q),YY(q),x(kk,j) + Rmax(j)*cos(0:0.1:2*pi),y(kk,j)+Rmax(j)*sin(0:0.1:2*pi))%XX1{kk,j}(k),YY1{kk,j}(k)) %&& ...
                        if inpolygon(XX(q),YY(q),XX1{kk,j}(k),YY1{kk,j}(k)) %&& ...
                            
                            %norm([x(kk,j)-XX(q),y(kk,j)-YY(q)])< Rmax(j)%-Rmax(j)/5
                            %             if kk>1 && norm([XX(q)-x(kk-1,j),YY(q)-y(kk-1,j)])< Rmax(j)
                            Uw{j}(q) = max(Uw{j}(q)-150,0.1);
                        else
                            Uw{j}(q) = Uw{j}(q);
                        end
                    end
                end
            end
            for j = 1:n_agents
                if length(str2num(beh(j)))>0 || beh(j) == "r"
                    R(kk,j)     = 0.35;%0.25
                    %             if kk >1
                    %                 if norm([x(kk,j)-cx2(kk-1,j),y(kk,j)-cy2(kk-1,j)])<0.2 %||  ...
                    %                     %norm([x(kk,j)-cx(kk-1,j),y(kk,j)-cy(kk-1,j)])<0.05
                    %                     R(kk,j) = R(kk-1,j) - dt;
                    %                     R(kk,j) = max(R(kk,j),0.05);
                    %                 else
                    %                     R(kk,j) = R(kk-1,j) -dt*(R(kk-1,j)-R(1,j));
                    %                 end
                    %             end
                elseif beh(j) == "c" || beh(j) == "e"
                    R(1,j)= 1000;
                    R(kk,j) = R(1,j);
                end
            end
            for j = 1:n_agents
                if beh(n)~= "e"
                    for n = 1:n_agents
                        Jo{kk}(j,n)        = 1/sum( UUfun(R(kk,n),XX1{kk,j},YY1{kk,j},xbar(n),ybar(n)));
                    end
                else
                    for n = 1:n_agents
                        Jo{kk}(j,n)        = 1/sum(Uw{j}(IN2{j}));
                    end
                end
                
                %         R(kk,j) = 0.1;
                %         R(kk,j) = min(max(R(kk,j),R(1,j)),1);
                %         U_base =  1;
                %         R_base =  1e5;
                %         U_0    = -1;
                %         syms r_x r_y
                %         r     = [r_x; r_y];
                %         Um(j)  = oldpotsum1(kk,x(:,j),y(:,j),r,R(kk,j),U_base,n_agents);
                %         Uf    = matlabFunction(Uw{j}, 'File','UUfun');
                if flag_expl ==0 || beh(j)~= "e"
                    mass        = sum( UUfun(R(kk,j),XX1{kk,j},YY1{kk,j},xbar(j),ybar(j)));
                    %             mass1       = sum( UUfun(R(kk,j),XX2{kk,j},YY2{kk,j},xbar(j),ybar(j)));
                    mass1       = sum( UUfun(R(kk,j),XX2{kk,j},YY2{kk,j},xbar(j),ybar(j)));
                    %             ratio  = mass1/mass;
                    JJJ(kk,j) = mass1;
                else
                    mass        = sum(Uw{j}(IN1{j}));
                    %             mass1       = sum(Uw
                    %             ratio  = mass1/mass;
                    %             mass        = sum(Uglob{kk}(IN1{j}));
                    
                    ue          = Uw{j}(IN1{j});
                    %             ue          = Uglob{kk}(IN1{j});
                    
                    mass1 = sum(Uw{j}(IN2{j}));
                    JJJ(kk,j) = mass1;
                    
                    %             mass1        = sum(Uglob{kk}(IN2{j}));
                    
                    ue1          = Uw{j}(IN2{j});
                    %             ue1          = Uglob{kk}(IN2{j});
                    
                    
                end
                CX           = 0;
                CY           = 0;
                CX1          = 0;
                CY1          = 0;
                for p = 1:length(XX1{kk,j})
                    if flag_expl ==0 || beh(j) ~= "e"
                        
                        CX = CX + XX1{kk,j}(p)*(UUfun(R(kk,j),XX1{kk,j}(p),YY1{kk,j}(p),xbar(j),ybar(j)));
                        CY = CY + YY1{kk,j}(p)*(UUfun(R(kk,j),XX1{kk,j}(p),YY1{kk,j}(p),xbar(j),ybar(j)));
                        
                        %             CX = CX + XX1{kk,j}(p)*(UUfun(XX1{kk,j}(p),YY1{kk,j}(p)));
                        %             CY = CY + YY1{kk,j}(p)*(UUfun(XX1{kk,j}(p),YY1{kk,j}(p)));
                        
                    else
                        CX = CX + XX1{kk,j}(p)*(ue(p));
                        CY = CY + YY1{kk,j}(p)*(ue(p));
                        
                        
                    end
                end
                for p = 1:length(XX2{kk,j})
                    if flag_expl == 0 || beh(j) ~= "e"
                        
                        %                 CX1 = CX1 + XX2{kk,j}(p)*(UUfun(R(kk,j),XX2{kk,j}(p),YY2{kk,j}(p),xbar(j),ybar(j)));
                        %                 CY1 = CY1 + YY2{kk,j}(p)*(UUfun(R(kk,j),XX2{kk,j}(p),YY2{kk,j}(p),xbar(j),ybar(j)));
                        
                        CX1 = CX1 + XX2{kk,j}(p)*(UUfun(R(kk,j),XX2{kk,j}(p),YY2{kk,j}(p),xbar(j),ybar(j)));
                        CY1 = CY1 + YY2{kk,j}(p)*(UUfun(R(kk,j),XX2{kk,j}(p),YY2{kk,j}(p),xbar(j),ybar(j)));
                        %             CX = CX + XX1{kk,j}(p)*(UUfun(XX1{kk,j}(p),YY1{kk,j}(p)));
                        %             CY = CY + YY1{kk,j}(p)*(UUfun(XX1{kk,j}(p),YY1{kk,j}(p)));
                        
                    else
                        CX1 = CX1 + XX2{kk,j}(p)*(ue1(p));
                        CY1 = CY1 + YY2{kk,j}(p)*(ue1(p));
                        
                        
                    end
                end
                
                
                %                CX1 = CX1 + XX2{kk,j}(p)*(UUfun(R(kk,j),XX2{kk,j}(p),YY2{kk,j}(p),xbar(j),ybar(j)));
                %                 CY1 = CY1 + YY2{kk,j}(p)*(UUfun(R(kk,j),XX2{kk,j}(p),YY2{kk,j}(p),xbar(j),ybar(j)));
                %                   CX1 = CX1 + XX2{kk,j}(p)*(ue1(p));
                %                 CY1 = CY1 + YY2{kk,j}(p)*(ue1(p));
                
                
                
                cx(kk,j)    = CX/(mass);
                cy(kk,j)    = CY/(mass);
                cx2(kk,j)   = CX1/mass1;
                cy2(kk,j)   = CY1/mass1;
                
                %         if norm([cx(kk,j)-x(kk,j),cy(kk,j)-y(kk,j)]) > 0.2
                %
                %             mass        = sum( UUfun(1,XX1{kk,j},YY1{kk,j},x(kk,j)+1*(cx(kk,j)-x(kk,j))/norm(cx(kk,j)-x(kk,j)),y(kk,j)+1*(cy(kk,j)-y(kk,j))/norm(cy(kk,j)-y(kk,j))) );
                %             CX          = 0;
                %             CY          = 0;
                %
                %             for p = 1:length(XX1{kk,j})
                %                 CX = CX + XX1{kk,j}(p)*(UUfun(1,XX1{kk,j}(p),YY1{kk,j}(p),x(kk,j)+1*(cx(kk,j)-x(kk,j))/norm(cx(kk,j)-x(kk,j)),y(kk,j)+1*(cy(kk,j)-y(kk,j))/norm(cy(kk,j)-y(kk,j))));
                %                 CY = CY + YY1{kk,j}(p)*(UUfun(1,XX1{kk,j}(p),YY1{kk,j}(p),x(kk,j)+1*(cx(kk,j)-x(kk,j))/norm(cx(kk,j)-x(kk,j)),y(kk,j)+1*(cy(kk,j)-y(kk,j))/norm(cy(kk,j)-y(kk,j))));
                %             end
                %             cx(kk,j)   = CX/(mass);
                %             cy(kk,j)   = CY/(mass);
                %         else
                %             cx(kk,j)    = CX/(mass);
                %             cy(kk,j)    = CY/(mass);
                %             mass        = sum( UUfun(1,XX1{kk,j},YY1{kk,j},x(kk,j)+1*cos(theta(kk,j)),y(kk,j)+1*sin(theta(kk,j) )));
                %             CX          = 0;
                %             CY          = 0;
                %
                %             for p = 1:length(XX1{kk,j})
                %                 CX = CX + XX1{kk,j}(p)*(UUfun(1,XX1{kk,j}(p),YY1{kk,j}(p),x(kk,j)+1*cos(theta(kk,j)),y(kk,j)+1*sin(theta(kk,j))));
                %                 CY = CY + YY1{kk,j}(p)*(UUfun(1,XX1{kk,j}(p),YY1{kk,j}(p),x(kk,j)+1*cos(theta(kk,j)),y(kk,j)+1*sin(theta(kk,j))));
                %             end
                %             cx(kk,j)   = CX/(mass);
                %             cy(kk,j)   = CY/(mass);
                %
                %         end
                %         cx(kk,j) = x(kk,j)+0.5*(cx(kk,j)-x(kk,j))/norm(cx(kk,j)-x(kk,j));
                %         cy(kk,j) = y(kk,j)+0.5*(cy(kk,j)-y(kk,j))/norm(cy(kk,j)-y(kk,j));
                
                
                if kk>1
                    x_dot(kk,j)     = x_dot(kk-1,j); %* cos(theta(kk-1,j)) ;
                    y_dot(kk,j)     = y_dot(kk-1,j); %* sin(theta(kk-1,j)) ;
                    %             theta_dot(kk,j) = omega(kk-1,j);
                end
            end
            for j = 1:n_agents
                
                %         x_dot(kk,j)      = vel(kk,j)*cos(theta(kk,j));
                %         y_dot(kk,j)      = vel(kk,j)*sin(theta(kk,j));
                %         theta_dot(kk,j)  = omega(kk,j);
                if kk == 1
                    x(kk+1,j)        = x(kk,j)    ;% + 5*(cx(kk,j)-x(kk,j)) * dt;
                    y(kk+1,j)        = y(kk,j)    ;% + 5*(cy(kk,j)-y(kk,j)) * dt;                    
                else
                    Cx = cx(kk,j)-x(kk,j);
                    Cy = cy(kk,j)-y(kk,j);
                    Cxm = cx(kk-1,j) - x(kk-1,j);
                    Cym = cy(kk-1,j) - y(kk-1,j);
                    if norm([(cx(kk,j)-cx(kk-1,j)),cy(kk,j)-cy(kk-1,j)]) == 0
                        xdot(kk,j) = 3*(cx(kk,j)-x(kk,j));%+((cx(kk,j)-cx(kk-1,j)));
                        ydot(kk,j) = 3*(cy(kk,j)-y(kk,j));%+((cy(kk,j)-cy(kk-1,j)));
                    elseif beh(j)~= "e"
                        xdot(kk,j) = 6*(cx(kk,j)-x(kk,j));%+.75*((cx(kk,j)-cx(kk-1,j)))/(norm([cx(kk,j)-cx(kk-1,j),cy(kk,j)-cy(kk-1,j)]));%+.1*randn()*(cy(kk,j)-cy(kk-1,j))/norm([cx(kk,j)-cx(kk-1,j),cy(kk,j)-cy(kk-1,j)]));
                        ydot(kk,j) = 6*(cy(kk,j)-y(kk,j));%+.75*((cy(kk,j)-cy(kk-1,j)))/(norm([cx(kk,j)-cx(kk-1,j),cy(kk,j)-cy(kk-1,j)]));%+.1*randn()*(cx(kk,j)-cx(kk-1,j))/norm([cx(kk,j)-cx(kk-1,j),cy(kk,j)-cy(kk-1,j)]));
%                       xdot(kk,j) = 3*(cx(kk,j)-x(kk,j))/norm([cy(kk,j)-y(kk,j),(cx(kk,j)-x(kk,j))]);
%                       ydot(kk,j) = 3*(cy(kk,j)-y(kk,j))/norm([cy(kk,j)-y(kk,j),(cx(kk,j)-x(kk,j))]);
                        
                    elseif beh(j) == "e"
                        
                        xdot(kk,j) = 3*(cx(kk,j)-x(kk,j))+1.75*((cx(kk,j)-cx(kk-1,j)))/(norm([cx(kk,j)-cx(kk-1,j),cy(kk,j)-cy(kk-1,j)]));%+.1*randn()*(cy(kk,j)-cy(kk-1,j))/norm([cx(kk,j)-cx(kk-1,j),cy(kk,j)-cy(kk-1,j)]));
                        ydot(kk,j) = 3*(cy(kk,j)-y(kk,j))+1.75*((cy(kk,j)-cy(kk-1,j)))/(norm([cx(kk,j)-cx(kk-1,j),cy(kk,j)-cy(kk-1,j)]));%+.1*randn()*(cx(kk,j)-cx(kk-1,j))/norm([cx(kk,j)-cx(kk-1,j),cy(kk,j)-cy(kk-1,j)]));
                        
                        %                   xdot(kk,j) = 6*(Cx)+.75*(Cx-Cxm)/(norm([Cx-Cxm,Cy-Cym]));%+.1*randn()*(cy(kk,j)-cy(kk-1,j))/norm([cx(kk,j)-cx(kk-1,j),cy(kk,j)-cy(kk-1,j)]));
                        %                   ydot(kk,j) = 6*(Cy)+.75*(Cy-Cym)/(norm([Cx-Cxm,Cy-Cym]));%+.1*randn()*(cx(kk,j)-cx(kk-1,j))/norm([cx(kk,j)-cx(kk-1,j),cy(kk,j)-cy(kk-1,j)]));
%                       xdot(kk,j) = 3*(cx(kk,j)-x(kk,j))/norm([cy(kk,j)-y(kk,j),(cx(kk,j)-x(kk,j))]);
%                       ydot(kk,j) = 3*(cy(kk,j)-y(kk,j))/norm([cy(kk,j)-y(kk,j),(cx(kk,j)-x(kk,j))]);
                    end  
                end
            end

            for j = 1:n_agents              
                if nonholo==0                    
                    x(kk+1,j)  = x(kk,j)     + xdot(kk,j)*dt ;%+x(kk,j)-x(kk-1,j)))*dt * dt;
                    y(kk+1,j)  = y(kk,j)     + ydot(kk,j)*dt ;
                    %hd(:,kk,j)     = [cx(kk,j)-x(kk,j)' cy(kk,j)-y(kk,j)]/norm([cx(kk,j)-x(kk,j) cy(kk,j)-y(kk,j)]);
                    hd(:,kk,j)     = [wp_input_vec(1,kk,j)-x(kk,j)' wp_input_vec(2,kk,j)-y(kk,j)]/norm([wp_input_vec(1,kk,j)-x(kk,j) wp_input_vec(2,kk,j)-y(kk,j)]);                    
                else                    
                    hd(:,kk,j)     = [wp_input_vec(1,kk,j)-x(kk,j)' wp_input_vec(2,kk,j)-y(kk,j)]/norm([wp_input_vec(1,kk,j)-x(kk,j) wp_input_vec(2,kk,j)-y(kk,j)]);
                    
                    if kk >1
                        [vel(kk,j),omega(kk,j)] = controller(theta(kk,j),  hd(:,kk,j), vel(kk-1,j), dt,x(kk,j),y(kk,j),cx2(kk,j),cy2(kk,j));
                        if norm([x(kk,j)-goal(1),y(kk,j)-goal(2)]) < 1
                            vel(kk,j) = norm([cx2(kk,j)-x(kk,j),cy2(kk,j)-y(kk,j)]);
                        end
                    end
                    x_dot(kk,j)      = vel(kk,j)*cos(theta(kk,j));
                    y_dot(kk,j)      = vel(kk,j)*sin(theta(kk,j));
                    theta_dot(kk,j)  = omega(kk,j);
                    x(kk+1,j)        = x(kk,j)     + vel(kk,j) * cos(theta(kk,j)) * dt;
                    y(kk+1,j)        = y(kk,j)     + vel(kk,j) * sin(theta(kk,j)) * dt;
                    theta(kk+1,j)    = theta(kk,j) + omega(kk,j) * dt;
                end
            end          
            h_waitbar = waitbar(kk/(length(time)-1));
            kkend = kk;
        end        
        close(h_waitbar);    
  
%% Post-processing (Plots)
obstacles            = zeros(length(x_obs),4);
drones               = gobjects(n_agents,9); % initialize array of plots
plot_obj             = gobjects(n_agents,7); % initialize array of plots
verCellHandle        = gobjects(length(time),1);
rect                 = gobjects(1,length(x_obs));
rectD                = gobjects(n_obsD,13);
txt                  = gobjects(n_obsD,1);
rand1                = gobjects(length(time),n_agents);
lev                  = gobjects(length(time),n_agents);
centr                = gobjects(length(time),1);
circ                 = gobjects(length(time),n_agents);
circ1                = gobjects(length(time),n_agents);
arrow                = gobjects(n_agents,n_obsD);
scat                 = gobjects(n_agents,1);
arrow1               = gobjects(n_agents,1);
delaunay             = gobjects(length(time),1);
apol                 = gobjects(n_agents,n_obsD);
cone_obj             = gobjects(1,4); % initialize array of plots
edges = gobjects(1,1) ;
colore = zeros(n_agents,3);

%close
figure('Name', 'Animation','units','normalized','outerposition',[0 0 1 1],'position',[0.13,0.364351851851852,0.4,0.5]);
set(gcf,'color',[1 1 1])
box on
xlabel('x (m)','interpreter','latex','fontsize',18)
ylabel('y (m)','interpreter','latex','fontsize',18)
axis([0 xlim(2)+0.5 0 ylim(2)+0.5])
grid on
hold on
set(gcf,'color','w')
obstacleX   = [obstacle1(:,1),obstacle1(:,2),obstacle1(:,3),obstacle1(:,4)];
for j = 1:length(x_obs)
    rect(j) = rectangle('Position',obstacleX(j,:),'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);
end

for kk = 1:1:length(time)
    
        for j = 1:n_agents
            if behtime(kk,j) == '1'
                colore(j,:) = DodgerBlue;
            elseif behtime(kk,j) == '2'
                colore(j,:) = DarkOrange;
            elseif behtime(kk,j) == 'r'
                colore(j,:) = [1;0;0];
            elseif behtime(kk,j) == 'c'
                colore(j,:) = [0;0;1];
            elseif behtime(kk,j) == 'e'
                colore(j,:) = [0;1;0];
            elseif behtime(kk,j) == '3'
                colore(j,:) = [1;1;0];
            elseif behtime(kk,j) == '4'
                colore(j,:) = [0;1;1];
            end
        end
    %     dot(hd(:,kk-1,m),([x(kk,n),y(kk,n)]-[x(kk,m),y(kk,m)])/norm([x(kk,m),y(kk,m)]-[x(kk,n),y(kk,n)]))
    %     dot(hd(:,kk-1,n),([x(kk,m),y(kk,m)]-[x(kk,n),y(kk,n)])/norm([x(kk,m),y(kk,m)]-[x(kk,n),y(kk,n)]))
     squeeze(Add(kk,:,:));
    [xo,yo] = gplot(ans,[x(kk,:)',y(kk,:)']);
    edges = plot(xo,yo,'m');
    for w = 1:n_obsD
        rectD(w,:) = human(obstacleD(kk,w,1),obstacleD(kk,w,2),VobsD(kk,w,1),VobsD(kk,w,2),sizeO);
    end
    for j = 1:n_agents
        bound = boundary(XX1{kk,j},YY1{kk,j},0.9);
        verCellHandle(kk,j)  = patch(x(kk,j),y(kk,j),[1 1 1/255],'FaceAlpha',0.1,'EdgeColor','k'); % use color i  -- no robot assigned yet
        set(verCellHandle(kk,j), 'XData',XX1{kk,j}(bound),'YData',YY1{kk,j}(bound));
    end


    for j = 1:n_agents
                                                centr(kk,j)             = plot( cx(kk,j), cy(kk,j),'bx');
%         traj =   plot(x(1:kk,j),y(1:kk,j),'-','linewidth',0.1,'Color','k');
        %         traj =   plot(x(1:kk,j),y(1:kk,j),'-','linewidth',0.01,'Color','k');
        %                                 circ2(j)     = circle(x(kk,j),y(kk,j),Rmaxt0(j),'none');
        if nonholo == 0
%             circ1(j)       = circle(x(kk,j),y(kk,j),1,[1 0.1 0.1]);
                       circ1(j)       = circle(x(kk,j),y(kk,j),.35*sizeA(j),colore(j,:));

        else
            plot_obj(j,:)  = plot_unicycle(x(kk,j), y(kk,j), theta(kk,j), 'k', sizeA(j));
        end
    end    
    pause
    box on    
    if video_flag == 1
        F(kk) = getframe(gcf);
    else
        drawnow
    end    
    if kk < length(time)-1
        delete(edges)
        delete(apol)
        delete(cone_obj)
        delete(drones)
        delete(verCellHandle)
        delete(centr)
        delete(rand1)
        delete(arrow)
        delete(arrow1)
        delete(rectD)
        delete(lev)
        delete(plot_obj)
        delete(circ);
        delete(circ1);
        % delete(wpplot);
        delete(delaunay);
        delete(txt);
%         delete(edges);
        % delete(traj);
        delete(scat);
        % delete(obs_circle);
%         delete(ar)
        %                 delete(circ2);
    else
        break
    end
end


