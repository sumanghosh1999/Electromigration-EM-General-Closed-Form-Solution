%% EM Closed-form Solution Evaluation on Binary-Tree Interconnect Network Geometry
    % Author: Suman Ghosh
    % Affiliation: Georgia Institute of Technology
    % Date: October, 2024
    % Description: Generates and saves EM data for several iterations using the closed-form
    %              solution on binary tree network geomtery. System specific
    %              parameters are adjustable inside the main simulation loop.

function main

close all;
clear all;
clc;

% Using default random state
rng('default')
s=rng;

sim_N=50; % number of iterations
data=cell(sim_N+1,14);
data(1,:)={'n','No.of Segments(2^n - 1)','L','W','H','j','c_i,1','c_i,2','Time taken for Parameter Calculation(sec)','EM-stress at t=T','Time taken for EM-stress evaluation(sec)','Void Location & Void Time','Time taken to calculate voiding(sec)','Total Time(sec)'};

for sim_n=1:sim_N
    clc;
    fprintf ("Iteration number: %d\n",sim_n);
    %% Simulation Parameters
    n = randi([3,8]);                                         % height of tree
    
    % Spatial parameters based on geometry
    L_list = (randi([100,500], 1, 2^n-1))*1e-6;               % lengths of (2^n - 1) segments
    W_list = 2.366e-8 * ones(1,2^n-1);                        % widths of (2^n - 1) segments
    H_list = 6.8e-8 * ones(1,2^n-1);                          % heights of (2^n - 1) segments
    L_divs = 10;
    j_list = (randi([0, 1], 1, 2^n-1) * 2 - 1).*(5e7 + (5e9 - 5e7) * rand(1, 2^n-1));   % current densities of (2^n - 1) segments
    
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
    
    c1_start=beta*j_list(1);                       % left terminal node
    c1_end_list=beta*j_list(end-2^(n-1)+1:end);    % right terminal nodes
    
    syms c1_ [1 2^n-1]
    
    c1_(1)=c1_start;
    c1_(end-2^(n-1)+1:end)=c1_end_list;
    
    equations=[];
    
    % Network-segment mapping
    A = adjacency_matrix_binTree(n);
    seg_map=generate_segment_map(2^n-1);
    
    
    for i=2:2^(n-1)
        segments_list=connected_segments(A,i);
        L=[];
        for j=1:length(segments_list)
            segment=segments_list(j,:);
            L=[L get_segment_number(segment(1),segment(2),seg_map)];
        end
        sum_l=0;
        sum_r=0;
        for k=1:length(L)
            sum_l = sum_l + W_list(L(k))*H_list(L(k))*c1_(L(k));
            sum_r = sum_r + W_list(L(k))*H_list(L(k))*j_list(L(k));
        end
        equations=[equations sum_l==beta*sum_r];
    end
    c1_(2:2^(n-1)-1)=struct2array(vpasolve(equations(1:end-1),c1_(2:2^(n-1)-1)));
    
    c1_=double(c1_);
    
    % c_i,2 calculation
    
    syms c2_ [1 2^n-1]
    syms X
    
    equations=[c2_(1)==X];
    
    for i=2:2^(n-1)
        segments_list=connected_segments(A,i);
        L=[];
        for j=1:length(segments_list)
            segment=segments_list(j,:);
            L=[L get_segment_number(segment(1),segment(2),seg_map)];
        end
        equations=[equations c1_(L(1))*L_list(L(1))+c2_(L(1))==c2_(L(2)) c1_(L(1))*L_list(L(1))+c2_(L(1))==c2_(L(3))];
    end
    
    c2_=struct2array(vpasolve(equations,c2_));
    phi_list=c2_-X;
    
    X_sol=-sum(W_list.*H_list.*(0.5.*c1_.*(L_list.^2) + phi_list.*L_list))/sum(W_list.*H_list.*L_list);
        
    c2_=phi_list+X_sol;
    c2_=double(c2_);
    
    t_calc=toc; % Time taken for parameter calculation
    
    fprintf('EM stress parameters calculated.\n');
    fprintf("Time taken for EM stress parameters calculation: %d seconds\n",t_calc);
    
    %% Stress Evaluation
    fprintf ("Evaluating EM stress at t = %d seconds timestep...\n",T);
    stress_list = [];
    for i = 1:2^n-1
        stress_list = [stress_list stresscalc(c1_(i),c2_(i),j_list(i), L_list(i), L_divs, alpha, beta, start_n, N, T)];
    end
    stress_list;
    
    t_EM_eval=toc-t_calc; % EM Stress Evaluation Time
    
    fprintf("EM stress evaluated at t = %d seconds timestep.\n",T);
    fprintf("Time taken for evaluating EM stress at t = %d seconds timestep:%d seconds \n",T,t_EM_eval);
    
    %% Voiding Location and Voiding Time Calculation
    fprintf('Calculating voiding locations and voiding timesteps...\n');
    
    void_params_list=[];
    for i =1:2^n-1
        void_params_list=[void_params_list void_calc(c1_(i),c2_(i),j_list(i), L_list(i), alpha, beta, start_n, N, sigma_lim,T,t_step)];
    end
    
    void_params_list;
    
    t_total=toc; % Total Time taken
    t_void_calc=t_total-t_calc-t_EM_eval; % Time taken for Voiding Location & Voiding Time Calculation
    
    fprintf("Voiding Locations and Voiding Timesteps calculated.\n");
    fprintf("Time taken for calculating voiding location and voiding timesteps: %d seconds \n",t_void_calc);
    fprintf("Total time taken: %d seconds \n",t_total);

    fprintf("Storing data for this iteration...\n");

    data(sim_n+1,:)={double(n),double(2^n - 1),mat2str(L_list),mat2str(W_list),mat2str(H_list),mat2str(j_list),mat2str([c1_]),mat2str(c2_),double(t_calc),mat2str(double(stress_list)),double(t_EM_eval),mat2str(double(void_params_list)),double(t_void_calc),double(t_total)};
    
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
filename=sprintf('ClosedFormSoln_binTree_Network_N=%d_Data_T=%dK%s',sim_N,temp,filetype);

%.xlsx format
writecell(data, filename);
disp('Data saved in .xlsx format!!!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plotting
%generate_plots(n, c1_, c2_, j_list, L_list, L_divs, alpha, beta, sigma_lim, start_n, N, t_step, T)


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


function A = adjacency_matrix_binTree(n)
 % Calculate the total number of nodes
    total_nodes = 2^n - 1;
    
    % Initialize the adjacency matrix
    A = zeros(total_nodes);
    
    % Iterate through each node and connect it to its children
    for i = 1:total_nodes
        left_child = 2*i;
        right_child = 2*i + 1;
        
        % If left child exists, create the connection
        if left_child <= total_nodes
            A(i, left_child) = 1;
            A(left_child, i) = 1;
        end
        
        % If right child exists, create the connection
        if right_child <= total_nodes
            A(i, right_child) = 1;
            A(right_child, i) = 1;
        end
    end
A1=zeros(size(A)+1);
A1(1,2)=1;A1(2,1)=1;
A1(2:end,2:end)=A;
A=A1;
end

function segments = connected_segments(A, i)
    % A: Adjacency matrix (nxn matrix)
    % i: The node index (integer)
    
    % Find all nodes connected to node i
    connected_nodes = find(A(i, :) == 1);  % Find indices where A(i,j) == 1
    
    % Create segments list as pairs of (i, j)
    segments = [repmat(i, length(connected_nodes), 1), connected_nodes'];
    for i=1:length(segments)
        if segments(i,1)>segments(i,2)
            L=segments(i,1);
            segments(i,1)=segments(i,2);
            segments(i,2)=L;
        end
    end
end

function segment_map = generate_segment_map(total_nodes)
    % total_nodes: the total number of nodes in the binary tree
    % segment_map: a matrix of pairs (parent, child) and their corresponding segment number

    segment_map = [];  % Initialize an empty matrix to store segment mappings
    segment_count = 1; % Start segment numbering from 1

    % Iterating over each node, starting from node 1
    for node = 1:total_nodes
        left_child = 2 * node;
        right_child = 2 * node + 1;

        % Checking if left child exists
        if left_child <= total_nodes
            segment_map = [segment_map; node, left_child, segment_count];  % Add the segment mapping
            segment_count = segment_count + 1;  % Increment segment number
        end

        % Check if right child exists
        if right_child <= total_nodes
            segment_map = [segment_map; node, right_child, segment_count];  % Add the segment mapping
            segment_count = segment_count + 1;  % Increment segment number
        end
    end
    segment_map=[[1 2 1];segment_map+1];
end

function k = get_segment_number(i, j, segment_map)
    % i: parent node
    % j: child node
    % segment_map: matrix containing all (i, j, k) mappings
    % k: the segment number corresponding to (i, j)

    % Finding the row in segment_map where the first two columns are i and j
    index = find(segment_map(:,1) == i & segment_map(:,2) == j);

    if ~isempty(index)
        k = segment_map(index, 3);  % Return the segment number k
    else
        k = -1;  % Return -1 if the segment (i,j) is not found
    end
end


function generate_plots(n, c1_, c2_, j_list, L_list, L_divs, alpha, beta, sigma_lim, start_n, N, t_step, T)
    % n: the number of layers
    f=0;
    figure_counter = 1;  % Counter for figure numbers

    for i = 1:n
        num_plots = 2^(i-1);  % Number of plots in the current layer
        figure(figure_counter);  % Create a new figure window
        figure_counter = figure_counter + 1;

        % Number of rows and columns for subplot layout
        rows = ceil(sqrt(num_plots));
        cols = ceil(num_plots / rows);

        for j = 1:num_plots
            subplot(rows, cols, j);  % Create a subplot
            hold on
            for t=0:t_step:T
                x=linspace(0,L_list(f+1),L_divs);
                y=stresscalc(c1_(i),c2_(f+1),j_list(f+1), L_list(f+1), L_divs, alpha, beta, start_n, N, t);
                plot(x,y);
            end
            yline(sigma_lim, '--r', ['y = ', num2str(sigma_lim)], 'LineWidth', 1.5);  % Add the y=constant line
            title(['Segment ', num2str(f+1)]);
            hold off;
            f=f+1;
        end

        % Adjust layout for better visibility
        sgtitle(['Layer ', num2str(i), ' with ', num2str(num_plots), ' segments']);
    end
end


