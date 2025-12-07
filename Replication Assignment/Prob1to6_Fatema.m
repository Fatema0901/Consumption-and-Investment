%%%%%%%%%%%%%%%%%%%%%%%%%%% Bellman equation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V(k,A,eps) = max { V_R(k,A,eps), V_NR(k,A,eps) }
% where
%   V_R (replacement)    = A*eps*k - F + beta * E V(1, A', eps')
%   V_NR(no replacement) = A*eps*k       + beta * E V((1-delta)k, A', eps')
%
% State variables:
%   k   : capital stock, discrete grid {1, (1-delta), (1-delta)^2, ...}
%   A   : aggregate productivity, 2-state Markov
%   eps : idiosyncratic productivity, iid uniform on [0.4,1.6]
%
% Choice variable:
%   z in {0,1}, where z=1 means replace capital (k'=1),
%   and z=0 means no replacement (k'=(1-delta)k).

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problem 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Variable and Parameter Setup and Preallocation

%%% Model Parameters
beta  = 0.9;
lambda = 0.75;          
F     = 0.2;
delta = 0.1;

%%% Aggregate Productivity Grid and Transition
A_vals = [1.25, 0.75];  % [A_H, A_L]
A_n    = length(A_vals);

% Markov transition for A: P(A' = Aj | A = Aj) = 0.9
P_A = [0.9 0.1;
       0.1 0.9];

%%% Idiosyncratic Productivity Grid (eps)
eps_n   = 20;                  % number of gridpoints for eps
eps_min = 0.4;
eps_max = 1.6;
eps_grid = linspace(eps_min, eps_max, eps_n);
pi_eps  = (1/eps_n) * ones(1, eps_n);   % iid uniform distribution

%%% Capital Grid k (vintages of the machine)
% New machine quality is normalized to 1.
% Each period of no replacement multiplies quality by (1-delta).
k_n   = 40;          % number of capital vintages
k_max = 1;           % new machine = 1

k_vec = k_max * (1 - delta).^(0:k_n-1)';   % [1, 0.9, 0.81, 0.729, ...]

%%% Revenue Function A * eps * k as 3D array: k x A x eps
revenue = NaN(k_n, A_n, eps_n);
for ik = 1:k_n
    for iA = 1:A_n
        for ie = 1:eps_n
            revenue(ik,iA,ie) = A_vals(iA) * eps_grid(ie) * k_vec(ik);
        end
    end
end

%% VFI Preallocations and Tolerances

mu      = 1e-6;                      % tolerance for convergence (For 1c)
maxits  = 1e5;                       % maximum number of iterations
count   = 0;
dif     = 1;

% (1b) Initial guess for value function V(k,A,eps)
V  = zeros(k_n, A_n, eps_n);

% (1d) Pre-allocate V_R, V_NR, and policy z(k,A,eps)
V_R  = NaN(k_n, A_n, eps_n);
V_NR = NaN(k_n, A_n, eps_n);
z    = zeros(k_n, A_n, eps_n);       % z=1 if replacement chosen, 0 otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problem 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main VFI Loop (Problen 2a, 2b, 2c)

tic
while dif > mu && count < maxits
    
    % Holding the previous iteration
    V_old = V;
    
    % Loop over current capital and aggregate productivity
    for ik = 1:k_n
        % indices for next-period k under replacement and no replacement
        k_rep_index  = 1;                       % k' = 1
        k_norep_index = min(ik+1, k_n);         % k' = (1-delta)k; truncated at top
        
        for iA = 1:A_n
            
            % ---------- Expected continuation values ----------
            % E V(k',A',eps') for replacement and no replacement.
            % Here eps is iid and independent of A.
            EV_R  = 0;
            EV_NR = 0;
            for iA_next = 1:A_n
                for ie_next = 1:eps_n
                    prob = P_A(iA,iA_next) * pi_eps(ie_next);
                    EV_R  = EV_R  + prob * V_old(k_rep_index , iA_next, ie_next);
                    EV_NR = EV_NR + prob * V_old(k_norep_index, iA_next, ie_next);
                end
            end
            
            % ---------- Current-period values for each eps ----------
            for ie = 1:eps_n
                % (2a) compute V_R^{i+1} and V_NR^{i+1}
                V_R(ik,iA,ie)  = revenue(ik,iA,ie) - F + beta * EV_R;
                V_NR(ik,iA,ie) = revenue(ik,iA,ie)      + beta * EV_NR;
                
                % (2b) update V^{i+1} and policy z^{i+1}
                if V_R(ik,iA,ie) >= V_NR(ik,iA,ie)
                    V(ik,iA,ie) = V_R(ik,iA,ie);
                    z(ik,iA,ie) = 1;
                else
                    V(ik,iA,ie) = V_NR(ik,iA,ie);
                    z(ik,iA,ie) = 0;
                end
            end
            
        end
    end
    
    % (2c) Checking convergence: max max max |V_i - V_{i+1}| < mu
    dif   = max(max(max(abs(V - V_old))));
    count = count + 1;
    
end
toc
fprintf('VFI converged in %d iterations, max diff = %.2e\n', count, dif);


%%%%%%%%%%%%%%%%% Problem 3: Cut-off level of capital %%%%%%%%%%%%%%%%%%%%
% To examine the policy function z(k,A,eps), I am looking for a cut-off level of capital k_bar such that for all k <= k_bar, the firm always replaces capital (z=1 for every A and eps).

% For each k index, checking if z(k,:,:)=1 for all A and eps
always_replace = squeeze( all( all( z==1 , 2) , 3) );  % k_n x 1 logical

% Largest k index where capital is ALWAYS replaced
cut_index = find(always_replace, 1, 'last');

if isempty(cut_index)
    warning('No cut-off found where replacement is always chosen at all states.');
else
    fprintf('Cut-off index for capital replacement: k_index = %d, k = %.4f\n', ...
            cut_index, k_vec(cut_index));
end

% Reducing the k dimension to 1 + cut_index
k_n = 150;                      % more dense grid
k_vec = (1-delta).^(0:k_n-1)';  % geometric decay from 1 downward

% Rebuilding the revenue array on the reduced k-grid
revenue = NaN(k_n, A_n, eps_n);
for ik = 1:k_n
    for iA = 1:A_n
        for ie = 1:eps_n
            revenue(ik,iA,ie) = A_vals(iA) * eps_grid(ie) * k_vec(ik);
        end
    end
end

%% --------- Repeating tasks 1 and 2 on reduced k-grid ------------------%%

% Reset convergence objects
mu      = 1e-6;
maxits  = 1e5;
count   = 0;
dif     = 1;

% New initial guess and preallocation on reduced grid
V   = zeros(k_n, A_n, eps_n);          % first guess for V(k,A,eps)
V_R = NaN(k_n, A_n, eps_n);
V_NR= NaN(k_n, A_n, eps_n);
z   = zeros(k_n, A_n, eps_n);          % policy: 1 = replace, 0 = no replace

tic
while dif > mu && count < maxits
    
    V_old = V;
    
    for ik = 1:k_n
        % next-period k indices under replacement and no replacement
        k_rep_index   = 1;                        % k' = 1
        k_norep_index = min(ik+1, k_n);           % k' = (1-delta)k, truncated
        
        for iA = 1:A_n
            
            % Expected continuation values under both actions
            EV_R  = 0;
            EV_NR = 0;
            for iA_next = 1:A_n
                for ie_next = 1:eps_n
                    prob = P_A(iA,iA_next) * pi_eps(ie_next);
                    EV_R  = EV_R  + prob * V_old(k_rep_index , iA_next, ie_next);
                    EV_NR = EV_NR + prob * V_old(k_norep_index, iA_next, ie_next);
                end
            end
            
            % Current-period values and policy decision for each eps
            for ie = 1:eps_n
                V_R(ik,iA,ie)  = revenue(ik,iA,ie) - F + beta * EV_R;
                V_NR(ik,iA,ie) = revenue(ik,iA,ie)      + beta * EV_NR;
                
                if V_R(ik,iA,ie) >= V_NR(ik,iA,ie)
                    V(ik,iA,ie) = V_R(ik,iA,ie);
                    z(ik,iA,ie) = 1;
                else
                    V(ik,iA,ie) = V_NR(ik,iA,ie);
                    z(ik,iA,ie) = 0;
                end
            end
            
        end
    end
    
    dif   = max(max(max(abs(V - V_old))));
    count = count + 1;
    
end
toc
fprintf('VFI on reduced k-grid converged in %d iterations, max diff = %.2e\n', count, dif);

%%%%%%%%%%%% Problem 4:Plotting Policy Function using Spy Plot %%%%%%%%%%%%
figure(3); clf
spy(z(:,:,round(eps_n/2)));          % middle productivity draw ε
xlabel('A index'); ylabel('k index');
title('Spy Plot of Replacement Policy (1 = Replace, 0 = Keep)'); 


%%%%%%%%%%%%%%%%%%%%%% Problem 5: Hazard Function %%%%%%%%%%%%%%%%%%%%%%%%%
% H(age, A) = Prob( replace | time since last replacement = age, A )
% Previously: k_vec(i) = (1-delta)^(i-1), so:
%   age = 0  -> k = 1          -> index ik = 1
%   age = 1  -> k = (1-delta)  -> index ik = 2
%   ...
% Following the paper Cooper, Haltiwanger & Power (1999) and plot ages 1,...,8 (one period or more since last replacement)
max_age = min(8, k_n-1);         % can't exceed available k-grid
age_grid = 1:max_age;

% Preallocating hazards: column vectors of length max_age
H_high = zeros(max_age,1);       % A = high state
H_low  = zeros(max_age,1);       % A = low state

% Identifying which index in A_vals is high/low
% (I built A_vals = [1.25, 0.75]; so index 1 = high, 2 = low)
iA_high = 1;
iA_low  = 2;

for age = 1:max_age
    ik = age + 1;   % capital index corresponding to this age
    
    % Taking mean over epsilon dimension of the 0/1 replacement decision
    H_high(age) = mean( squeeze( z(ik, iA_high, :) ) );
    H_low(age)  = mean( squeeze( z(ik, iA_low , :) ) );
end

%%%%%%%%%%%%% Plot: Hazard Function of Capital Replacement %%%%%%%%%%%%%%%%
figure(5); clf
plot(age_grid, H_high, '-','LineWidth',1.5); hold on
plot(age_grid, H_low , '--','LineWidth',1.5);
hold off

xlabel('Time Since Last Replacement');
ylabel('Probability of Replacement');
title('Theoretical Hazard for Machine Replacement');
legend('High State','Low State','Location','SouthEast');
grid on
axis([1 max_age 0 1]);   % y-axis between 0 and 1 as in the paper 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problem 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I am plotting about 40 periods of:
%   - output of the firm
%   - capital stock of the firm

sims = 200;                  % total simulation length (I will plot first 40)
Tplot = 40;

% ---- Simulate aggregate state A_t (2-state Markov) ----
% A index: 1 = high, 2 = low
mc_A   = dtmc(P_A);              % Markov chain object (like your prof's code)
A_ind  = simulate(mc_A, sims);   % sims x 1 vector of {1,2}

% ---- Simulate idiosyncratic shock eps_t (i.i.d. uniform over grid) ---- %
eps_ind = randi(eps_n, sims, 1);  % each is an index 1,...,eps_n

% ---- Containers for endogenous variables ----
k_ind   = zeros(sims,1);     % index of k on grid
k_path  = zeros(sims,1);     % actual capital level
y_path  = zeros(sims,1);     % output A*eps*k
z_path  = zeros(sims,1);     % replacement decision

% Initial conditions: start with new machine and some A,eps draw
k_ind(1)  = 1;               % k = 1 → brand new
k_path(1) = k_vec(k_ind(1));

% ---- Simulation loop ---- %
for t = 1:(sims-1)
    
    ik = k_ind(t);
    iA = A_ind(t);
    ie = eps_ind(t);
    
    % Current shocks
    A_t   = A_vals(iA);
    eps_t = eps_grid(ie);
    
    % Policy: replace (1) or not (0)
    z_t = z(ik, iA, ie);
    z_path(t) = z_t;
    
    % Output in period t (same revenue function used in VFI)
    y_path(t) = A_t * eps_t * k_vec(ik);
    
    % Law of motion for capital index:
    % if replace → k' = 1 (index 1)
    % if no replace → k' = (1-delta)*k, which corresponds to next index
    if z_t == 1
        k_ind(t+1) = 1;                     % replacement: reset to new
    else
        k_ind(t+1) = min(ik + 1, k_n);      % no replacement: age machine
    end
    
    % Record next-period capital level
    k_path(t+1) = k_vec(k_ind(t+1));
end

% Last period output 
ik_last = k_ind(sims);
A_last  = A_vals(A_ind(sims));
eps_last= eps_grid(eps_ind(sims));
y_path(sims) = A_last * eps_last * k_vec(ik_last);

%% ------ Plots: Output and Capital for One Firm (first 40 periods) -----%%
figure(6); clf

subplot(2,1,1)
plot(1:Tplot, y_path(1:Tplot), 'LineWidth',1.5);
xlabel('Time');
ylabel('Output');
title('Output Path for One Simulated Firm');
grid on

subplot(2,1,2)
plot(1:Tplot, k_path(1:Tplot), 'LineWidth',1.5);
xlabel('Time');
ylabel('Capital Stock');
title('Capital Path for One Simulated Firm');
grid on
