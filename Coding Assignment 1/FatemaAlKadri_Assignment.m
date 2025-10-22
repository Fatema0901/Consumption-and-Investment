clear; clc; close all;
rng(2);

%% -----------------------------
%  1) PARAMETERS (given)
% ------------------------------
beta  = 0.96;      % discount factor
gamma = 1.3;       % CRRA risk aversion
r     = 0.04;      % interest rate
rho   = 0.9;       % persistence of income y
sigma = 0.04;      % standard deviation of income shocks (eps ~ N(0,sigma^2))

%% -----------------------------------------------------------
%  2) STATE SPACES: assets 'a' (choose) and income 'y' (Tauchen)
% ------------------------------------------------------------
Na   = 75;         % number of asset grid points 
amin = 0.0;        % borrowing constraint (cannot go below 0)
amax = 50.0;       % upper bound
agrid = linspace(amin, amax, Na)';   % assets as a COLUMN vector for easy indexing

Ny = 7;                        % income grid size 
mu = 0; m = 3;                 % Tauchen inputs: mean 0, span +/- m*std of stationary y
% Tauchen returns: ygrid (levels of AR(1) state) and P (Ny x Ny transition probs)
[ygrid, P] = Tauchen(Ny, mu, rho, sigma, m);

% Here exp(y) is actual income level in the budget constraint.
ey = exp(ygrid);               % pointwise income levels for each y-state

%% ----------------------------------------------------
%  3) UTILITY helper as an inline function handle
% ----------------------------------------------------
% CRRA utility with gamma != 1. Return -Inf if c<=0 to enforce c>0 constraint.
u = @(c) ((c).^(1-gamma)) ./ (1-gamma);
u_neg_inf = -1.0e10;   % large negative fallback in case I avoid -Inf numerics

%% ----------------------------------------------------
%  4) INITIALIZE VALUE FUNCTION AND POLICY STORAGE
% ----------------------------------------------------
V      = zeros(Na, Ny);      % initial guess: zeros everywhere 
V_old  = V;                   % previous iterate
polAix = ones(Na, Ny);       % indices of optimal a' choices (store argmax)
polA   = agrid(polAix);      % levels of a' (will overwrite after convergence)
polC   = zeros(Na, Ny);      % consumption policy c(a,y)

tol      = 1e-9;             % tolerance on norm difference
max_iter = 5000;             % safety cap (per hint: prevent infinite loop)
iter     = 0;
diffV    = inf;

%% ----------------------------------------------------
%  5) VALUE FUNCTION ITERATION (Bellman)
% ----------------------------------------------------
% Bellman: V(a,y) = max_{a' in grid} { u(c) + beta * E_y'[ V(a',y') | y ] }
% where c = a + exp(y) - a'/(1+r), and c>0.

% Precompute the expectation term structure:
% For each y (row j), and for each a' index k, we need E[V(a'_k, y')|y_j]
% That is sum_{y'} P(j,y') * V_old(k,y'). I am computing inside the loop to stay transparent.

while diffV > tol && iter < max_iter
    V_old = V;                 % keep a copy to measure convergence
    iter  = iter + 1;

    % Loop over current states (a_i, y_j)
    for j = 1:Ny              % income state index
        for i = 1:Na          % asset state index

            % Candidate values for all possible a' choices from the grid:
            % Compute consumption for each candidate a'_k
            c_all = agrid(i) + ey(j) - agrid./(1+r);   % vectorized over k = 1..Na

            % Enforce c>0. For invalid c, assign very low utility.
            util = u(c_all);
            util(c_all <= 0) = u_neg_inf;  % using big negative so max ignores them

            % Expected continuation value for each candidate a'_k:
            % For fixed current y=j, the expectation over y' is row j of P.
            EV = (V_old * P(j, :)')';  % WRONG SHAPE if used directly; keep it explicit below
            % To stay crystal clear, compute EV for each k one-by-one:
            EV_k = zeros(Na,1);
            for k = 1:Na
                % dot of row j in P with row k of V_old
                EV_k(k) = P(j, :) * V_old(k, :)';
            end

            % Right-hand side of Bellman for each a'_k
            RHS = util + beta * EV_k;

            % Take the best a' (argmax over k)
            [V(i,j), polAix(i,j)] = max(RHS);
        end
    end

    % Convergence check: norm of the change in V (vector 2-norm; hint allows 'norm')
    diffV = norm(V - V_old);
    if mod(iter,50)==0
        % Light progress ping so I can tell in class that it’s iterating:
        disp(['Iter ', num2str(iter), '  ||V - V_old|| = ', num2str(diffV)]);
    end
end

if iter == max_iter
    disp('WARNING: Reached max_iter before hitting tolerance.');
else
    disp(['Converged in ', num2str(iter), ' iterations.  Diff = ', num2str(diffV)]);
end

% Recover policy levels for a' and c using the argmax indices
polA = agrid(polAix);                          % a'(a,y)
for j = 1:Ny
    for i = 1:Na
        polC(i,j) = agrid(i) + ey(j) - polA(i,j)/(1+r);  % c = a + exp(y) - a'/(1+r)
    end
end

%% ----------------------------------------------------
%  6) PLOTS (Value function and Policies)
% ----------------------------------------------------
% (a) Plot the converged value function V(a,y) for all y.
figure;
hold on;
for j = 1:Ny
    plot(agrid, V(:,j), 'LineWidth', 1.2);  % explain: 'plot' draws lines
end
xlabel('Assets a'); ylabel('Value V(a,y)');
title('Converged Value Function for all y-states');
legend(cellstr("y="+string(round(ygrid,3))),'Location','southeast');
grid on; hold off;

% (b) Policy for a'(a,y)
figure;
hold on;
for j = 1:Ny
    plot(agrid, polA(:,j), 'LineWidth', 1.2);
end
xlabel('Assets a'); ylabel('Policy a''(a,y)');
title('Policy Function for Savings a'' across y-states');
legend(cellstr("y="+string(round(ygrid,3))),'Location','southeast');
grid on; hold off;

% (c) Policy for c(a,y)
figure;
hold on;
for j = 1:Ny
    plot(agrid, polC(:,j), 'LineWidth', 1.2);
end
xlabel('Assets a'); ylabel('Consumption c(a,y)');
title('Policy Function for Consumption across y-states');
legend(cellstr("y="+string(round(ygrid,3))),'Location','southeast');
grid on; hold off;

%% ----------------------------------------------------
%  7) SIMULATION (1000 periods; discard first 500)
% ----------------------------------------------------
T = 1000;
a_path  = zeros(T,1);      % simulated assets (on-grid, via policy)
y_path  = zeros(T,1);      % simulated continuous y via AR(1)
iy_path = zeros(T,1);      % index of nearest ygrid each period
c_path  = zeros(T,1);

% Initial states (start from mid asset grid and mean income)
iy_path(1) = round((Ny+1)/2);
y_path(1)  = ygrid(iy_path(1));
ia         = round((Na+1)/2);      % start in the middle of the asset grid
a_path(1)  = agrid(ia);

% Generate shocks eps_t ~ N(0, sigma^2) using randn (standard normal).
% (Not in the hint list; explanation: randn draws i.i.d. N(0,1) numbers.)
eps = sigma * randn(T,1);

for t = 1:T-1
    % Income AR(1): y_{t+1} = rho*y_t + eps_{t+1}
    y_path(t+1) = rho * y_path(t) + eps(t+1);

    % Map continuous y_t to the nearest index on ygrid (nearest-neighbor projection)
    [~, iy] = min(abs(ygrid - y_path(t)));
    iy_path(t) = iy;

    % Given current (a_index = ia, y_index = iy), choose next a' via policy index
    ia_next = polAix(ia, iy);
    a_next  = agrid(ia_next);

    % Consumption from budget: c_t = a_t + exp(y_t) - a_{t+1}/(1+r)
    c_t = agrid(ia) + exp(y_path(t)) - a_next/(1+r);
    % Enforce positivity in simulation (should already be ensured by policy):
    if c_t <= 0
        c_t = NaN; % flag if ever violated (should not)
    end

    % Record and advance:
    c_path(t) = c_t;
    a_path(t+1) = a_next;
    ia = ia_next;
end
% Fill last-period indices cleanly
[~, iy] = min(abs(ygrid - y_path(T))); iy_path(T) = iy;
c_path(T) = agrid(ia) + exp(y_path(T)) - a_path(T)/(1+r);

% Discard burn-in
burn = 500;
y_sim  = y_path(burn+1:end);
a_sim  = a_path(burn+1:end);
c_sim  = c_path(burn+1:end);

% Tile-style time-series view using tiledlayout (explained since not in hint list):
% 'tiledlayout' organizes multiple subplots neatly; 'nexttile' selects each panel.
figure;
tiledlayout(3,1);
nexttile; plot(exp(y_sim)); ylabel('exp(y)'); title('Simulated Income (levels)');
nexttile; plot(a_sim);      ylabel('a');      title('Simulated Assets');
nexttile; plot(c_sim);      ylabel('c');      title('Simulated Consumption');
xlabel('Time (post-burn-in)');

% Standard deviation of simulated consumption:
std_c = std(c_sim, 'omitnan');
disp(['Std dev of simulated consumption (post-burn-in): ', num2str(std_c)]);

%% ----------------------------------------------------
%  8) NOTES FOR PART (d) 
% ----------------------------------------------------
% (a) If borrowing constraint were zero (i.e., allow a<0), consumption smoothing improves,
%     so simulated c typically becomes *less* volatile => std(c) tends to fall.
% (b) If risk aversion gamma doubled, the agent dislikes variation in c more strongly,
%     saving policy smooths c more aggressively => std(c) tends to fall.
% (c) If natural interest rate doubles, saving becomes more attractive; policy may
%     tilt toward higher a', which can dampen c fluctuations for precautionary reasons,
%     but transition dynamics can raise c sensitivity. Qualitatively: std(c) ambiguous,
%     often *slightly lower* in stationary distribution with higher buffer-stock saving.
% (d) If income volatility doubles (sigma up), precautionary saving rises, but shocks
%     hitting the budget are larger; absent perfect insurance, std(c) tends to *rise*.
