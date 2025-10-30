clear; clc; close all;

%% -----------------------------
% PARAMETERS 
% ------------------------------
beta  = 0.96; % discount factor
r     = 0.04;    % interest rate
rho   = 0.9;     % wage persistence
sigma = 0.15;    % std dev of wage shock
wbar  = 2.5;     % long-run mean wage
phi   = 2;       % Frisch elasticity
n_target = 40/168;  % target steady-state hours (≈ 0.238)

%% ----------------------------------------------------------
% 1) Caliberating "Omega" so that steady-state n = 40/168 at wbar
%    From intratemporal FOC: n*(w) = (w / Omega)^phi  (clipped to [0,1])
%    => Omega = wbar / (n_target)^(1/phi)
% -----------------------------------------------------------
Omega = wbar / (n_target)^(1/phi);

%% ----------------------------------------------------------
% 2) Discretizing STATE SPACES
%    (a) Wage grid via Tauchen
%    (b) Asset grid with natural borrowing limit
% -----------------------------------------------------------

% (a) Wage grid: Ny points around wbar, span m * stationary std
Ny = 7;                                         % number of wage states
m  = 3;                                         % span +/- 3 stationary std devs
[wgrid, P] = Tauchen(Ny, wbar, rho, sigma, m);  % wgrid: Nyx1 column vector, P: NyxNy matrix

% (b) Natural borrowing limit x:
%     Bank assumes worst-case income every period: lowest wage and max hours.
%     With no consumption, present value of income floor = w_min * 1 / r.
w_min = min(wgrid);     % Minimum wage from the wage grid
n_max = 1;              % Maximum hours worked
%     Natural debt limit x = - (w_min * n_max) / r, with n_max = 1.
x     = - w_min / r;     % a' > x  (putting a grid starting from x)

% Asset grid (choose upper bound large enough to avoid top binding)
Na    = 10;              % number of asset points
a_min = x;
a_max = 50;
agrid = linspace(a_min, a_max, Na)';   % column vector of current assets a

%% ----------------------------------------------------------
% 3) LABOR SUPPLY: closed-form from FOC, n*(w) = (w/Omega)^phi, clipped
%    Clip to [0,1] because hours are a fraction of available time.
% -----------------------------------------------------------
n_star = (wgrid / Omega) .^ phi;        % Nyx1 vector
n_star = min(max(n_star, 0), 1);        % enforce 0 <= n <= 1

% Precomputing the disutility term for each wage state:
% D(w) = Omega * n^(1+1/phi) / (1+1/phi), used inside the log.
D = Omega * (n_star .^ (1 + 1/phi)) / (1 + 1/phi);   % Nyx1

%% ----------------------------------------------------------
% 4) UTILITY HELPER (CRRA in effective consumption with log form)
%    u(c,n) = log( c - D(w) ), where D(w) computed above per w-state.
% Composing the "inside of log" manually to check positivity and avoid NaNs
% Because Infeasible c -> big negative utility.
% -----------------------------------------------------------
u_neg = -1.0e10;  % numerical -Inf

%% ----------------------------------------------------------
% 5) POLICY-FUNCTION ITERATION (PFI) SCAFFOLD
%    Maintaining:
%      V      : Na x Ny value function
%      polAix : Na x Ny indices for optimal a' (on the a-grid)
%      polA   : Na x Ny levels for a' (agrid(polAix))
%      polC   : Na x Ny optimal consumption
%    PFI loop alternates:
%       (i) Policy Improvement: (for each w, choosing a' maximizing U + EV)
%      (ii) Howard Policy Evaluation: (holding a'(a,w) fixed, updating V)
% -----------------------------------------------------------
V      = zeros(Na, Ny);        % initial guess for V
V_old  = V;                    % storage for convergence check
polAix = repmat((1:Na)', 1, Ny); % replicate matrix (1:Na)' 1 time vertically & Ny time horizontally
% Here I am assuming households choose their next period'd asset equal to current period'd asset
polA   = agrid(polAix);        % optimal a' levels
polC   = zeros(Na, Ny);        % optimal consumtion level (will fill after convergence)

iter     = 0;                  % outer iteration counter
diffV    = inf;                % initialize difference to enter loop
tol      = 1e-9;               % convergence tolerance on ||V - V_old||
max_iter = 1000;               % safety cap on outer PFI iterations
H_eval   = 20;                 % Howard evaluation steps per policy update

%% ----------------------------------------------------------
% 6) VALUE FUNCTIN ITERATION
% ----------------------------------------------------------
% Precomputing utility U(a,a',w): log( (1+r)(a + w*n(w)) - a' - D(w) )
U = -inf(Na,Na,Ny);                                % This creates 3D array for 
    % 1. number of current asset grid points 
    % 2. number of possible future asset choices
    % 3. number of wage states
for iw = 1:Ny                                      % Starts wage loop over each discrete wage state w_i 
    res = (1+r)*(agrid + wgrid(iw)*n_star(iw));    % Total Resource before choosing a'
    cons = res - agrid';                           % Consumption = Total Resource - Next Period Asset(a')
    good = cons > D(iw);                           % Only Taking the Positive Consumptions
    tmp  = -inf(Na,Na);                            % Creating a temporary matrix for asset choices 
    tmp(good) = log(cons(good) - D(iw));           % u(c,n) with positivity check
    U(:,:,iw) = tmp;                               % Stores the result for i-th wage state & adds to the Utility array
end


for iter = 1:max_iter          % Start iterating from 1 to maximum (1000)
    EV = beta*(V*P');          % Expected Discounted Continuation Value

    % Policy Improvement: for each w, choosing a' maximizing U + EV
    V_new = zeros(Na,Ny);      % Placeholder for updated value function
    for iw = 1:Ny              % Loop over each wage state w_i
        [V_new(:,iw), polAix(:,iw)] = max( U(:,:,iw) + EV(:,iw)', [], 2 );
    end

    % Howard Policy Evaluation: hold a'(a,w) fixed, update V
    for h = 1:H_eval           % Repeat policy evaluation H_eval = 20 times before updating the policy again.
        for iw = 1:Ny
            ap  = polAix(:,iw);                           % chosen a' indices
            c   = (1+r)*(agrid + wgrid(iw)*n_star(iw)) - agrid(ap); % implied c(a,w)
            V_new(:,iw) = log(c - D(iw)) + beta*( V_new(ap,:)*P(iw,:)' );
        end
    end

    % Convergence
    if max(abs(V_new(:) - V(:))) < tol % finds the largest absolute change across all grid points,
                                       %if value function stops changing by more than tolerance level, exit loop
        
        V = V_new; break;              % if the condition is true, the value function has converged 
    end
    V = V_new;                         % if the condition is not true goes to new next iteration (iter + 1)
end
%% ----------------------------------------------------------
% Final answer for (a)
% ----------------------------------------------------------
polA = agrid(polAix);                                      % (1) a'(a,w)
polC = (1+r)*(agrid + wgrid'.*n_star') - polA;             % (2) c(a,w)
V    = V;                                                  % (3) V(a,w) (already converged)

%% ----------------------------------------------------------
% Answer for (b): Plot V(a,w) for all wage states
% -----------------------------------------------------------
figure('Name','Value Function by Wage State'); 
tiledlayout('flow');
for iw = 1:Ny
    nexttile;
    plot(agrid, V(:,iw), 'LineWidth', 1.25);
    xlabel('a');
    ylabel(sprintf('V(a,w_%d)', iw));
    title(sprintf('w state %d', iw));
    grid on;
end 

%% ----------------------------------------------------------
% Simulations (Answer for (c))
% -----------------------------------------------------------
sims = 1000;
y_sim = simulate(dtmc(P), sims);     % Markov chain of wage indices
a_index = round(Na/2) * ones(sims+1,1);  % starts from mid asset index

% Preallocate
c_sim = zeros(sims,1);
a_sim = zeros(sims+1,1);
n_sim = zeros(sims,1);

for t = 1:sims
    iw = y_sim(t);                          % current wage state index
    n_sim(t) = n_star(iw);                  % labor at that wage
    c_sim(t) = (1+r)*(agrid(a_index(t)) + wgrid(iw)*n_sim(t)) ...
                - agrid( polAix(a_index(t), iw) );   % consumption
    a_index(t+1) = polAix(a_index(t), iw);           % next asset index
    a_sim(t+1) = agrid(a_index(t+1));                % actual asset level
end

% ---- Plots ----
figure;
subplot(3,1,1)
plot(sims/2+1:sims, wgrid(y_sim(sims/2+1:sims)))
xlabel('Time'); ylabel('Wage'); title('Simulated Wage')

subplot(3,1,2)
plot(sims/2+1:sims, c_sim(sims/2+1:sims))
xlabel('Time'); ylabel('Consumption'); title('Simulated Consumption')

subplot(3,1,3)
plot(sims/2+1:sims, a_sim(sims/2+1:sims))
xlabel('Time'); ylabel('Assets'); title('Simulated Assets')


%% ----------------------------------------------------------
% (d) Standard deviation of simulated hours and qualitative discussion
% -----------------------------------------------------------

% Dropping first half of simulation as burn-in
n_trim = n_sim(sims/2+1:sims);

% Computing standard deviation of hours worked
std_n = std(n_trim);
fprintf('Standard deviation of simulated n (post burn-in): %.6f\n', std_n);

disp('Comparative statics on std(n):')
disp(' (1) Relax borrowing constraint -> std(n) ↓')
% Agents can borrow when income is low, so consumption smoothing becomes
% easier. Std(n) decreases. 
disp(' (2) Double risk aversion -> std(n) ↑')
% Households dislike consumption risk more strongly. Std(n) increases.
disp(' (3) Double Frisch elasticity -> std(n) ↑')
% Labor supply becomes more responsive to changes in the after-tax wage or
% consumption–leisure tradeoff. Hours fluctuate more in response to wage shocks. 
% Std(n) increases.
disp(' (4) Double wage volatility -> std(n) ↑')
% Wage shocks directly introduce more variation in income. Both consumption and hours become more volatile.
% Std(n) increases.
