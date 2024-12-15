Rhoop = 3; % the radius of the hoop
r0 = 1; % the equilibrial length of the springs
kappa = 1; % the spring constant
Nnodes = 21;
A = zeros(Nnodes,Nnodes); % spring adjacency matrix
% vertical springs
for k = 1 : 3
    A(k,k+4) = 1;
end
for k = 5 : 7  
    A(k,k+5) = 1;
end
for k = 10 : 12  
    A(k,k+5) = 1;
end
for k = 15 : 17  
    A(k,k+4) = 1;
end
% horizontal springs
for k = 4 : 7
    A(k,k+1) = 1;
end
for k = 9 : 12  
    A(k,k+1) = 1;
end
for k = 14 : 17  
    A(k,k+1) = 1;
end
% symmetrize
Asymm = A + A';
%indices of nodes on the hoop
ind_hoop = [0,3,8,13,18,19,20,17,12,7,2,1] + 1;
Nhoop = length(ind_hoop);
% indices of free nodes (not attached to the hoop)
ind_free = [4,5,6,9,10,11,14,15,16] + 1;
Nfree = length(ind_free);
% list of springs
[springs_0,springs_1] = ind2sub([Nnodes,Nnodes],find(A));
springs = [springs_0';springs_1'];

Nsprings = length(springs_0);
% maps indices of nodes to their indices in the gradient vector

%% Initialization

% Initial angles for the nodes are uniformly distributed around the range of 2*pi
% startting from theta0 and going counterclockwise
theta0 = 2*pi/3;
theta = theta0 + linspace(0,2*pi,Nhoop+1)';
theta(end) = [];
% Initial positions
pos = zeros(Nnodes,2);
pos(ind_hoop,1) = Rhoop*cos(theta);
pos(ind_hoop,2) = Rhoop*sin(theta);
pos(ind_free,1) = [-1.,0.,1.,-1.,0.,1.,-1.,0.,1.]';
pos(ind_free,2) = [1.,1.,1.,0.,0.,0.,-1.,-1.,-1.]'; 

% Initiallize the vector of parameters to be optimized
vec = [theta;pos(ind_free,1);pos(ind_free,2)]; % a column vector with 30 components

draw_spring_system(pos,springs,Rhoop,ind_hoop,ind_free);

gradient = @(vec)compute_gradient(vec,Asymm,r0,kappa,Rhoop,ind_hoop,ind_free);
func = @(vec)Energy(vec,springs,r0,kappa,Rhoop,ind_hoop,ind_free);

%% optimization

%% Method 1: Gradient Descent
max_iterations = 1000; % maximum number of iterations
learning_rate = 0.02;  % learning rate
convergence_tolerance = 1e-8; % tolerance

current_vector = vec; 
energy_history = zeros(max_iterations, 1); 
gradient_norm_history = zeros(max_iterations, 1); 

for iteration = 1:max_iterations
    current_energy = func(current_vector); 
    current_gradient = gradient(current_vector); 
    energy_history(iteration) = current_energy; 
    gradient_norm_history(iteration) = norm(current_gradient); 
    
    if gradient_norm_history(iteration) < convergence_tolerance 
        break;
    end
    
    current_vector = current_vector - learning_rate * current_gradient; 
end

[theta_gd, pos_gd] = vec_to_pos(current_vector, Rhoop, ind_hoop, ind_free); 

figure;
subplot(2, 1, 1);
plot(1:iteration, energy_history(1:iteration), 'LineWidth', 2); 
xlabel('Iteration');
ylabel('Energy');
title('Energy vs Iteration (Gradient Descent)');
subplot(2, 1, 2);
semilogy(1:iteration, gradient_norm_history(1:iteration), 'LineWidth', 2);
xlabel('Iteration');
ylabel('Gradient Norm'); 
title('Gradient Norm vs Iteration (Gradient Descent)');

draw_spring_system(pos_gd, springs, Rhoop, ind_hoop, ind_free);
title('Spring System (Gradient Descent)');

fprintf('Gradient Descent Results:\n');
fprintf('Final Energy: %f\n', func(current_vector)); 
fprintf('Final Gradient Norm: %.2e\n', norm(gradient(current_vector)));
fprintf('Positions of Nodes:\n');
disp(pos_gd);

%% Method 2: Conjugate Gradient
max_iterations = 1000; % maximum number of iterations
tolerance = 1e-8;      % tolerance for convergence

current_vector_cg = vec; % initialize the vector for conjugate gradient
energy_history_cg = zeros(max_iterations, 1); % preallocate energy array
gradient_norm_history_cg = zeros(max_iterations, 1); % preallocate gradient norm array

current_gradient_cg = gradient(current_vector_cg);
direction = -current_gradient_cg;

for iteration = 1:max_iterations
    current_energy = func(current_vector_cg); % compute energy
    previous_gradient = current_gradient_cg; % store previous gradient
    
    energy_history_cg(iteration) = current_energy; % store energy
    gradient_norm_history_cg(iteration) = norm(current_gradient_cg); % store gradient norm
    
    if gradient_norm_history_cg(iteration) < tolerance % check for convergence
        break; % exit loop if converged
    end
    
    alpha = 1; % initial step size
    c = 1e-4; % constant for Wolfe condition
    rho = 0.9; % reduction factor for step size
    while func(current_vector_cg + alpha * direction) > current_energy + c * alpha * current_gradient_cg' * direction
        alpha = rho * alpha; % reduce step size
    end
    
    current_vector_cg = current_vector_cg + alpha * direction; % update vector
    current_gradient_cg = gradient(current_vector_cg); % compute new gradient
    beta = (current_gradient_cg' * current_gradient_cg) / (previous_gradient' * previous_gradient); % compute beta
    direction = -current_gradient_cg + beta * direction; % update direction
end

[theta_cg, pos_cg] = vec_to_pos(current_vector_cg, Rhoop, ind_hoop, ind_free); % convert vector to positions

figure; % create a new figure
subplot(2, 1, 1); % create first subplot
plot(1:iteration, energy_history_cg(1:iteration), 'LineWidth', 2); % plot energy
xlabel('Iteration'); % x-axis label
ylabel('Energy'); % y-axis label
title('Energy Reduction Over Iterations (Conjugate Gradient)'); % updated title

subplot(2, 1, 2); % create second subplot
semilogy(1:iteration, gradient_norm_history_cg(1:iteration), 'LineWidth', 2); % plot gradient norm
xlabel('Iteration'); % x-axis label
ylabel('Gradient Norm'); % y-axis label
title('Gradient Norm Convergence Over Iterations (Conjugate Gradient)'); % updated title

draw_spring_system(pos_cg, springs, Rhoop, ind_hoop, ind_free); % draw the spring system
title('Final Spring System Configuration After Conjugate Gradient'); % updated title

fprintf('Results from Conjugate Gradient:\n'); % updated print output
fprintf('Final Energy Value: %f\n', func(current_vector_cg)); % updated print output
fprintf('Final Gradient Norm Value: %.2e\n', norm(gradient(current_vector_cg))); % updated print output
fprintf('Final Positions of Nodes:\n'); % updated print output
disp(pos_cg); % display positions

%%
function draw_spring_system(pos,springs,R,ind_hoop,ind_free)
% draw the hoop 
figure;
hold on;
t = linspace(0,2*pi,200);
plot(R*cos(t),R*sin(t),'linewidth',5,'color','r');
% plot springs
Nsprings = size(springs,2);
for k = 1 : Nsprings
    j0 = springs(1,k);
    j1 = springs(2,k);
    plot([pos(j0,1),pos(j1,1)],[pos(j0,2),pos(j1,2)],'linewidth',3,'color','k');
end
% plot nodes
plot(pos(ind_hoop,1),pos(ind_hoop,2),'.','Markersize',100,'Color',[0.5,0,0]);
plot(pos(ind_free,1),pos(ind_free,2),'.','Markersize',100,'Color','k');
set(gca,'Fontsize',20);
daspect([1,1,1]);
end

%% 
function grad = compute_gradient(vec,Asymm,r0,kappa,R,ind_hoop,ind_free)
    [theta,pos] = vec_to_pos(vec,R,ind_hoop,ind_free);    
    Nhoop = length(ind_hoop);
    g_hoop = zeros(Nhoop,1); % gradient with respect to the angles of the hoop nodes
    Nfree = length(ind_free);
    g_free = zeros(Nfree,2); % gradient with respect to the x- and y-components of the free nodes
    for k = 1 : Nhoop
        ind = find(Asymm(ind_hoop(k),:)); % index of the node adjacent to the kth node on the hoop
        rvec = pos(ind_hoop(k),:) - pos(ind,:); % the vector from that adjacent node to the kth node on the hoop
        rvec_length = norm(rvec); % the length of this vector
        g_hoop(k) = (rvec_length - r0)*R*kappa*(rvec(1)*(-sin(theta(k))) + rvec(2)*cos(theta(k)))/rvec_length;
    end
    for k  = 1 : Nfree
        ind = find(Asymm(ind_free(k),:)); % indices of the nodes adjacent to the kth free node
        Nneib = length(ind);
        for j = 1 : Nneib
            rvec = pos(ind_free(k),:) - pos(ind(j),:); % the vector from the jth adjacent node to the kth free node 
            rvec_length = norm(rvec);  % the length of this vector
            g_free(k,:) = g_free(k,:) + (rvec_length - r0)*R*kappa*rvec/rvec_length;
        end
    end
    % return a single 1D vector
    grad = [g_hoop;g_free(:,1);g_free(:,2)];
end

%%
function E = Energy(vec,springs,r0,kappa,R,ind_hoop,ind_free)
    [~,pos] = vec_to_pos(vec,R,ind_hoop,ind_free);
    Nsprings = size(springs,2);
    E = 0.;
    for k =1 : Nsprings
        j0 = springs(1,k);
        j1 = springs(2,k);
        rvec = pos(j0,:) - pos(j1,:);
        rvec_length = norm(rvec);       
        E = E + kappa*(rvec_length - r0)^2;
    end
    E = E*0.5;
end

%%
function [theta,pos] = vec_to_pos(vec,R,ind_hoop,ind_free)
    Nhoop = length(ind_hoop);
    Nfree = length(ind_free);
    Nnodes = Nhoop + Nfree;
    theta = vec(1:Nhoop);
    pos = zeros(Nnodes,2);
    pos(ind_hoop,1) = R*cos(theta);
    pos(ind_hoop,2) = R*sin(theta);
    % positions of the free nodes
    pos(ind_free,1) = vec(Nhoop+1:Nnodes);
    pos(ind_free,2) = vec(Nnodes+1:end); 
end