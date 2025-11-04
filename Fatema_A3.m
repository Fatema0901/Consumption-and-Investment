clear; clc;
%% ------------------- PARAMETERS ----------------------------------------
beta   = 0.96;           % Discount factor
gamma  = 3.5;            % Risk aversion
r      = 0.04;           % Interest rate
rho    = 0.92;           % Persistence of income
sigma  = 0.05;           % Std. dev. of income shocks
pu     = 0.05;           % Prob. of unemployment from any employment state
pe     = 0.75;           % Prob. of re-employment to lowest income state

%% ------------------- INCOME GRID ---------------------------------------
Y_n    = 5;              % Number of employed income states
sd_y   = 2.5 * sigma;    % Number of std deviations around mean

% Tauchen discretization (or Rouwenhorst if you prefer)
[Y, P] = Tauchen(Y_n, 0, rho, sigma, sd_y);

% Add unemployment state: y = -Inf (exp(y)=0)
Y      = [-Inf; Y(:)];
Y_n    = length(Y);

% Augment transition matrix
P_u = zeros(Y_n);        % Initilizing initialize an empty (K+1)×(K+1) transition matrix full of zeros
P_u(1,1)   = 1 - pe;     % Prob of staying unemployed in the next period (1-p_e) 
P_u(1,2)   = pe;         % Prob of re-employment
for i = 2:Y_n            % i runs over from 2 to Y_n = All employed income grid points
    P_u(i,1)     = pu;   % Prob of becoming unemployed
    P_u(i,2:end) = (1 - pu) * P(i-1,:); 
end
P = P_u; clear P_u; 

% Income in levels
Y_exp = exp(Y); 
Y_exp(1) = 0; % unemployed income = 0

%% ------------------- ASSET GRID ----------------------------------------
a_n   = 500;                   
a_max = 4 * max(Y_exp);
a_min = 0;                     % no borrowing (tight limit)
A     = linspace(a_min, a_max, a_n)';

%% ------------------- MAIN EGM ITERATION --------------------------------
tol     = 1e-6;         % tolerance value
maxits  = 1e4;          % maximum iteration number 
count   = 0;            % iteration counter
diff    = 1;            % initial gap measure 

% Initial guess: consume all cash-on-hand (on the w=A grid)
C       = repmat(A,1,Y_n);                 % c(w,y)
eps_c   = 1e-8;                            % to keep c>0

% Helpers
inv_u_p = @(x) x.^(-1/gamma);           % helps Euler equation to solve for consumption

tic
while diff > tol && count < maxits      % Keeps iterating untill reaches maximum iteration

    % --------- Expectation of marginal utility (fully vectorized) ---------
    % Next-period cash-on-hand for every (a', y'):
    %   w' = a' + exp(y')
    w_next = A + Y_exp';                   % next-period cash-on-hand value

    % c'(w', y') for all columns y' at once (loop only on y', small)
    C_next = zeros(a_n, Y_n);              % Stores next-period consumption for every combination of future assets and income.
    for yp = 1:Y_n
        % Interpolate along assets dimension; vectorized over all a'
        C_next(:,yp) = interp1(A, C(:,yp), w_next(:,yp), 'pchip', 'extrap');
    end

    % Expected marginal utility across y' for each current y
    MU_next = max(C_next, eps_c).^(-gamma);      % [a_n x Y_n]
    Emu     = MU_next * P.';                     % [a_n x Y_n]   Emu(:,y)=Σ_{y'} MU_next(:,y')*P(y,y')'

    % ---------- Euler inversion and endogenous grid ----------
    C_euler = inv_u_p( beta*(1+r) * Emu );       % [a_n x Y_n]
    % w = c + a'/(1+r)  for each column y (broadcasting a' across y)
    W_endo  = C_euler + A/(1+r);                 % Endogenous grid of cash-on-hand value

    % ---------- Interpolate back to exogenous w-grid (here using A) -------
    C_new = C;                                   % New C consumptoin matrix like old dimension
    for y = 1:Y_n                                % 
        w_vec = real(W_endo(:,y));
        c_vec = real(C_euler(:,y));
        [w_sorted, idx] = sort(w_vec);
        c_sorted = c_vec(idx);
        % Borrowing-limit branch: for w below smallest endogenous point
        % c = w - a_min/(1+r). Here a_min= A(1).
        c_bc = @(w) max(w - A(1)/(1+r), eps_c);

        % Interpolate onto the target grid A (serving as w-grid)
        c_interp = interp1(w_sorted, c_sorted, A, 'pchip', 'extrap');

        % Guard the far left with borrowing-limit rule
        left_mask = A < w_sorted(1);
        if any(left_mask)
            c_interp(left_mask) = c_bc(A(left_mask));
        end

        % Enforce positivity
        C_new(:,y) = max(c_interp, eps_c); 
    end

    % ---------- Convergence and update ----------
    diff  = max(max(abs(C_new - C))); % measures how much the consumption function has changed between the current iteration and the previous one
    C     = C_new;                    % the updated consumption matrix after the latest EGM loop
    count = count + 1;                % Increment iteration counter
end
toc
fprintf('EGM converged in %d iterations, max diff = %.2e\n', count, diff);

%% ----------------------- RESULTS (part a) ------------------------------
% Consumption policy on the cash-on-hand grid (here we used A as w-grid)
figure(1); clf
plot(A, C, 'LineWidth', 1.5);
xlabel('Cash-on-hand  w'); ylabel('Consumption  c(w,y)');
title('Consumption Policy Function by Income State');
lgd = arrayfun(@(k) sprintf('y_{%d}',k-1), 1:Y_n, 'UniformOutput', false);
lgd{1} = 'Unemployed';
legend(lgd, 'Location','SouthEast'); grid on 

%% ----------------------- RESULTS (part b) ------------------------------
% Plotting c(w,y) for each income state on the cash-on-hand grid w(a,y)=a+exp(y)
% Note: For each y, the x-axis grid is w_y = A + Y_exp(y)

figure(2); clf; hold on
for y = 1:Y_n
    w_y = A + Y_exp(y);               % cash-on-hand grid for state y
    % ensures increasing x for plotting 
    [w_sorted, idx] = sort(w_y);
    c_sorted = C(idx, y);
    plot(w_sorted, c_sorted, 'LineWidth', 1.5);
end
% 45-degree reference: c = w (only for visual comparison on a common w range)
wmin = min(A + Y_exp(1));                   % usually unemployment state
wmax = max(A + Y_exp(end));
wref = linspace(wmin, wmax, 200)';
plot(wref, wref, ':', 'LineWidth', 1.0);    % reference only

xlabel('Cash-on-hand  w(a,y) = a + e^{y}')
ylabel('Consumption  c(w,y)')
title('Consumption Policy Function by Income State on the w(a,y) Grid')
lgd = arrayfun(@(k) sprintf('y_{%d}',k-1), 1:Y_n, 'UniformOutput', false);
lgd{1} = 'Unemployed';
legend([lgd, {'c=w (reference)'}], 'Location','SouthEast');
grid on; box on
hold off
%% ----------------------- SIMULATION (part c) ----------------------------
T = 5000;                         % number of periods to simulate
a_sim = zeros(T,1);               % assets
y_sim = zeros(T,1);               % income state index
c_sim = zeros(T,1);               % consumption
rng(123);                         % for reproducibility

% start from middle income state and zero assets
y_sim(1) = ceil(Y_n/2);
a_sim(1) = 0;

for t = 1:T-1
    % current state
    y_idx = y_sim(t);

    % compute current cash-on-hand
    w_t = a_sim(t) + Y_exp(y_idx);

    % interpolate consumption for given (w_t, y_t)
    c_t = interp1(A + Y_exp(y_idx), C(:,y_idx), w_t, 'pchip', 'extrap');
    c_sim(t) = max(c_t, eps_c);

    % savings choice from budget constraint
    a_next = (w_t - c_sim(t)) * (1+r);
    a_sim(t+1) = min(max(a_next, A(1)), A(end));  % stay within grid bounds

    % income transition
    P_cum = cumsum(P(y_idx,:));
    y_sim(t+1) = find(rand <= P_cum, 1);
end
c_sim(T) = interp1(A + Y_exp(y_sim(T)), C(:,y_sim(T)), ...
                   a_sim(T) + Y_exp(y_sim(T)), 'pchip', 'extrap');

% Construct income in levels
y_level = Y_exp(y_sim);

%% ----------------------- CORRELOGRAM -----------------------------------
maxlag = 4;
[ccf,lags] = xcorr(c_sim - mean(c_sim), y_level - mean(y_level), maxlag, 'coeff');

figure(3); clf
stem(lags, ccf, 'filled');
xlabel('Lags'); ylabel('Correlation');
title('Correlogram: Consumption vs Income (up to 4 lags)');
grid on; box on;

%% ----------------------- PLOTTING CORRELOGRAM (part d) ----------------------------
figure(4); clf
bar(lags, ccf, 'FaceColor',[0.3 0.5 0.8]); hold on
plot([0 0], ylim, 'k--'); % vertical line at zero lag
xlabel('Lag');
ylabel('Correlation Coefficient');
title('Correlogram: Consumption vs Income (up to 4 lags)');
grid on; box on;
