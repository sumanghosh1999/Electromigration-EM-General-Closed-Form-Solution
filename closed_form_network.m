%% EM Closed-form Solution Evaluation on (m-1-n) Interconnect Network Geometry
    % Author: Suman Ghosh
    % Affiliation: Georgia Institute of Technology
    % Date: July, 2024
    % Description: Generates and saves EM data for several iterations using the closed-form
    %              solution on (m-1-n) network geomtery. System specific
    %              parameters are adjustable inside the main simulation loop.

function main

close all;
clear all;
clc;

% Using default random state
rng('default')
s=rng;

sim_N=100; % number of iterations
data=cell(sim_N+1,20);
data(1,:)={'m','n','L','W','H','j','c_i,1','c_i,2','Injection point','Injection current','Time taken for Parameter Calculation(sec)','EM-stress at t=T(LEFT)','EM-stress at t=T(CENTRAL)','EM-stress at t=T(RIGHT)','Time taken for EM-stress evaluation(sec)','Void Location & Void Time(LEFT)','Void Location & Void Time(CENTRAL)','Void Location & Void Time(RIGHT)','Time taken to calculate voiding(sec)','Total Time(sec)'};

for sim_n=1:sim_N
    clc;
    fprintf ("Iteration number: %d\n",sim_n);
    %% Simulation Parameters
    
    % geometrical parameters
    m=randi([1,10])                    % number of segments on left of central segment
    n=randi([1,10])                    % number of segements on right of central segment
    
    % Spatial parameters based on geometry
    L_list = (randi([100,500], 1, m+n+1))*1e-6;               % lengths of (m + 1 + n) segments
    W_list = 2.366e-8 * ones(1,m+n+1);                        % widths of (m + 1 + n) segments
    H_list = 6.8e-8 * ones(1,m+n+1);                          % heights of (m + 1 + n) segments
    L_divs = 10;
    j_list = (randi([0, 1], 1, m+n+1) * 2 - 1).*(5e7 + (5e9 - 5e7) * rand(1, m+n+1));   % current densities of (m + 1 + n) segments
    
    % Current injection
    I_inj = 1e-4; % injection current magnitude
    
    I_inj = I_inj * (randi([0, 1]) * 2 - 1); % injection current directionality (injection / withdrawl)
    inj_loc = randi([0, 1]); % injection location
    
    if inj_loc==0      % left-central node injection
        inj_L = I_inj;
        inj_R = 0;
    elseif inj_loc==1  % right-central node injection
        inj_L = 0;
        inj_R = I_inj;
    end

    % Updating current densities
    A_list=W_list.*H_list;
    I_list=j_list.*A_list;

    if inj_loc==0      % left-central node injection
        V_L=(inj_L +sum(I_list(1:m+1)))/sum(A_list(1:m+1)./L_list(1:m+1));
        j_list(1:m+1)=V_L./L_list(1:m+1);
    elseif inj_loc==1  % right-central node injection
        V_R=(inj_R +sum(I_list(m+1:m+n+1)))/sum(A_list(m+1:m+n+1)./L_list(m+1:m+n+1));
        j_list(m+1:m+n+1)=V_R./L_list(m+1:m+n+1);
    end

    % Simulation time parameters
    T=86400;
    t_step=8640;
    
    % Partial infinite sum parameters
    N = 3;       % no. of terms in partial sum
    start_n = 1; % n=1 start point of partial sum
    
    %%%%%%%%%%%%
    % T = 300K %
    %%%%%%%%%%%%
    temp=300;
    
    % compact model parameters
    alpha0= 3.92987e-13;
    beta= 544;
    sigma_lim= 3.1478e+08;
    gamma= -0.7736;
    delta=2.07e+08;
    zeta=5.7864e+05;
    
    alpha=alpha0;
    
    tic; % start timer

    %% Stress parameters calculation
    fprintf('Calculating stress parameters...\n');

    % c_i,1 calculation
    
    c1_m_list=beta*j_list(1:m);                % left terminal nodes
    c1_n_list=beta*j_list(m+2:length(j_list)); % right terminal nodes
    c1_central=(1/(W_list(m+1)*H_list(m+1)))*(sum(beta*W_list(1:m+1).*H_list(1:m+1).*j_list(1:m+1))-sum(W_list(1:m).*H_list(1:m).*c1_m_list));  % left-central node
    c1_central=(1/(W_list(m+1)*H_list(m+1)))*(sum(beta*W_list(m+1:length(W_list)).*H_list(m+1:length(H_list)).*j_list(m+1:length(j_list)))-sum(W_list(m+2:length(W_list)).*H_list(m+2:length(W_list)).*c1_n_list));  % right-central node
    
    % c_i,2 calculation
    
    syms c2_ [1 m+n+1]
    syms X
    
    equations=[c2_(m)==X];
    
    for i=1:m
        equations=[equations c1_m_list(i)*L_list(i)+c2_(i) == c2_(m+1)];
    end
    
    for i=1:n
        equations = [equations c2_(m+1+i) == c1_central*L_list(m+1)+c2_(m+1)];
    end
   
    %pretty(equations)
    c2_sol=struct2array(vpasolve(equations,c2_));
    phi_list=c2_sol-X;
    
    sum_left = sum( (W_list(1:m).*H_list(1:m)).*(0.5.*c1_m_list.*(L_list(1:m).^2) + phi_list(1:m).*L_list(1:m)) );
    sum_center = (W_list(m+1)*H_list(m+1)).*(0.5.*c1_central*(L_list(m+1)^2) + phi_list(m+1).*L_list(m+1));
    sum_right = sum( (W_list(m+2:length(W_list)).*H_list(m+2:length(H_list))).*(0.5.*c1_n_list.*(L_list(m+2:length(L_list)).^2) + phi_list(m+2:length(phi_list)).*L_list(m+2:length(L_list))) );
    
    X_sol=-(sum_left+sum_center+sum_right)/sum(W_list.*H_list.*L_list);
    
    c2_sol=phi_list+X_sol;
    
    t_calc=toc; % Time taken for parameter calculation

    fprintf('EM stress parameters calculated.\n');
    fprintf("Time taken for EM stress parameters calculation: %d seconds\n",t_calc);

    %% Stress Evaluation
    fprintf ("Evaluating EM stress at t = %d seconds timestep...\n",T);
    m_seg_stress_list = [];
    for i = 1:m
        m_seg_stress_list = [m_seg_stress_list stresscalc(c1_m_list(i),c2_sol(i),j_list(i), L_list(i), L_divs, alpha, beta, start_n, N, T)];
    end
    
    central_seg_stress = stresscalc(c1_central,c2_sol(m+1),j_list(m+1), L_list(m+1), L_divs, alpha, beta, start_n, N, T);
    
    n_seg_stress_list = [];
    for i = 1:n
        n_seg_stress_list = [n_seg_stress_list stresscalc(c1_n_list(i),c2_sol(m+1+i),j_list(m+1+i), L_list(m+1+i), L_divs, alpha, beta, start_n, N, T)];
    end
    
    t_EM_eval=toc-t_calc; % EM Stress Evaluation Time

    fprintf("EM stress evaluated at t = %d seconds timestep.\n",T);
    fprintf("Time taken for evaluating EM stress at t = %d seconds timestep:%d seconds \n",T,t_EM_eval);

    %% Voiding Location and Voiding Time Calculation
    fprintf('Calculating voiding locations and voiding timesteps...\n');

    m_void_params_list=[];
    for i =1:m
        m_void_params_list=[m_void_params_list void_calc(c1_m_list(i),c2_sol(i),j_list(i), L_list(i), alpha, beta, start_n, N, sigma_lim,T,t_step)];
    end
    
    central_void_params = void_calc(c1_central,c2_sol(m+1),j_list(m+1), L_list(m+1), alpha, beta, start_n, N, sigma_lim,T,t_step);
    
    n_void_params_list=[];
    for i =1:n
        n_void_params_list = [n_void_params_list void_calc(c1_n_list(i),c2_sol(m+1+i),j_list(m+1+i), L_list(m+1+i), alpha, beta, start_n, N, sigma_lim,T,t_step)];
    end
    
    t_total=toc; % Total Time taken
    t_void_calc=t_total-t_calc-t_EM_eval; % Time taken for Voiding Location & Voiding Time Calculation

    fprintf("Voiding Locations and Voiding Timesteps calculated.\n");
    fprintf("Time taken for calculating voiding location and voiding timesteps: %d seconds \n",t_void_calc);
    fprintf("Total time taken: %d seconds \n",t_total);

    fprintf("Storing data for this iteration...\n");

    data(sim_n+1,:)={double(m),double(n),mat2str(L_list),mat2str(W_list),mat2str(H_list),mat2str(j_list),mat2str([c1_m_list c1_central c1_n_list]),mat2str(double(c2_sol)),inj_loc,I_inj,double(t_calc),mat2str(double(m_seg_stress_list)), mat2str(double(central_seg_stress)), mat2str(double(n_seg_stress_list)),double(t_EM_eval),mat2str(double(m_void_params_list)), mat2str(double(central_void_params)), mat2str(double(n_void_params_list)),double(t_void_calc),double(t_total)};
    
    fprintf("Data stored for this iteration.\n");


% end of sim_N loop
end
clc
fprintf('Simulation data stored. Ready to write...\n');
sz_data=size(data);
fprintf('Data size: %d x %d\n',sz_data(1),sz_data(2));

%% Wrting data to .xlsx file
fprintf('Writing data to .xlsx file...\n')

filetype='.xlsx';
filename=sprintf('ClosedFormSoln_withInjection_m1n_Network_N=%d_Data_T=%dK%s',sim_N,temp,filetype);

%.xlsx format
writecell(data, filename);
disp('Data saved in .xlsx format!!!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting
% Calculate the number of rows in the grid (maximum of m or n)
rows = max(m, n);

% Create figure
figure;


% Define the height of each subplot
subplot_height_m = 0.8 / m; % reduce total height to leave space below
subplot_height_n = 0.8 / n; % reduce total height to leave space below

% Plot subplots in the first column
for i = 1:m
    pos = [0.1, 0.1 + (m-i) * subplot_height_m, 0.25, subplot_height_m - 0.12]; % [left, bottom, width, height]
    subplot('Position', pos);
    hold on
    f=0;
    for t=0:t_step:T
        x=linspace(0,L_list(i),L_divs);
        y=stresscalc(c1_m_list(i),c2_sol(i),j_list(i), L_list(i), L_divs, alpha, beta, start_n, N, t);
        plot(x,y);
        idx_above = find(y >= sigma_lim);
        idx_below = find(y <= -sigma_lim);
        plot(x(idx_above), y(idx_above), 'ro', 'MarkerFaceColor', 'r');
        plot(x(idx_below), y(idx_below), 'ro', 'MarkerFaceColor', 'r');
    end
    hold off
    title(sprintf('Left Segment %d\nLength: %d um, Width: %.4f um, Height: %.4f um\nj = %d A/m^2', i, L_list(i)/1e-6,W_list(i)/1e-6,H_list(i)/1e-6,j_list(i)));
    xlim([0 L_list(i)])
end

% Plot the subplot in the second column (centered vertically)
pos = [0.4, 0.3, 0.25, 0.4]; % [left, bottom, width, height]
subplot('Position', pos);
hold on
f=0;
for t=0:t_step:T
    x=linspace(0,L_list(m+1),L_divs);
    y=stresscalc(c1_central,c2_sol(m+1),j_list(m+1), L_list(m+1), L_divs, alpha, beta, start_n, N, t);
    plot(x,y);
    idx_above = find(y >= sigma_lim);
    idx_below = find(y <= -sigma_lim);
    plot(x(idx_above), y(idx_above), 'ro', 'MarkerFaceColor', 'r');
    plot(x(idx_below), y(idx_below), 'ro', 'MarkerFaceColor', 'r');
end
hold off
title(sprintf('Central Segment\nLength: %d um, Width: %.4f um, Height: %.4f um\nj = %d A/m^2', L_list(m+1)/1e-6,W_list(m+1)/1e-6,H_list(m+1)/1e-6,j_list(m+1)));
xlim([0 L_list(m+1)])

% Plot subplots in the third column
for i = 1:n
    pos = [0.7, 0.1 + (n-i) * subplot_height_n, 0.25, subplot_height_n - 0.12]; % [left, bottom, width, height]
    subplot('Position', pos);
    hold on
    f=0;
    for t=0:t_step:T
        x=linspace(0,L_list(m+1+i),L_divs);
        y=stresscalc(c1_n_list(i),c2_sol(m+1+i),j_list(m+1+i), L_list(m+1+i), L_divs, alpha, beta, start_n, N, t);
        plot(x,y);
        idx_above = find(y >= sigma_lim);
        idx_below = find(y <= -sigma_lim);
        plot(x(idx_above), y(idx_above), 'ro', 'MarkerFaceColor', 'r');
        plot(x(idx_below), y(idx_below), 'ro', 'MarkerFaceColor', 'r');
    end
    hold off
    title(sprintf('Right Segment %d\nLength: %d um, Width: %.4f um, Height: %.4f um\nj = %d A/m^2', i, L_list(m+1+i)/1e-6,W_list(m+1+i)/1e-6,H_list(m+1+i)/1e-6,j_list(m+1+i)));
    xlim([0 L_list(m+1+i)])
end


% end of main
end

%% Functions
%% EM Stress Calculation Function
function result=stresscalc(c1,c2,j, L, L_divs, alpha, beta, start_n, N, t)
    x = linspace(0, L, L_divs);
    result=zeros(1,length(x));
    for i=1:length(x)
        syms n
        f1=(((-1)^n-1)/(n^2))*cos(n*pi*x(i)/L)*exp((-alpha*n^2*pi^2*t)/(L^2));
        f2=(((-1)^n-1)/(n^2))*sin(n*pi*x(i)/L)*exp((-alpha*n^2*pi^2*t)/(L^2));
        result(i)=c1*x(i) + c2 - ((2*c1*L)/pi^2)*symsum(f1, n, start_n, N) + (2/pi)*( c2 + (c1*L)/2 )*symsum(f2, n, start_n, N);
    end
end

%% Void Location & Void Time Calculation Function
function result=void_calc(c1, c2, j, L, alpha, beta, start_n, N, sigma_lim, T, t_step)
    y0 = c2;
    yL = c1 * L + c2;
    if abs(y0) > abs(yL)
        max_point_x = 0;
        max_point_y = y0;
    else
        max_point_x = L;
        max_point_y = yL;
    end

    loc=max_point_x;

    syms void_t
    syms n
    f1=(((-1)^n-1)/(n^2))*cos(n*pi*loc/L)*exp((-alpha*n^2*pi^2*void_t)/(L^2));
    f2=(((-1)^n-1)/(n^2))*sin(n*pi*loc/L)*exp((-alpha*n^2*pi^2*void_t)/(L^2));
    f = @(t) double(subs(c1*loc + c2 - ((2*c1*L)/pi^2)*symsum((((-1)^n-1)/(n^2))*cos(n*pi*loc/L)*exp((-alpha*n^2*pi^2*void_t)/(L^2)), n, start_n, N) + (2/pi)*( c2 + (c1*L)/2 )*symsum((((-1)^n-1)/(n^2))*sin(n*pi*loc/L)*exp((-alpha*n^2*pi^2*void_t)/(L^2)), n, start_n, N),void_t,t));
    
    void_time=T;
    for t=0:t_step/10:T
        if f(t)>=sigma_lim
            void_time=t;
            break;
        end
    end
    result=[loc,void_time];
end


