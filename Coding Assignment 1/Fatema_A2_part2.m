clear; clc; close all;

%% ----------------------------------------------------------
% PARAMETERS  
% -----------------------------------------------------------
beta  = 0.96;     % discount factor
alpha = 0.36;     % capital share
rhoz   = 0.9;     % persistence of productivity

% Assuming
delta0 = 0.05;    % baseline depreciation (used in δ(u))
theta = 1.5;      % elasticity of depreciation with respect to utilization   
sigma_z = 0.02;   % std dev of productivity shocks

%% ----------------------------------------------------------
% 1) DISCRETIZING STATE SPACE 
% -----------------------------------------------------------
Nz = 7;
m  = 3;
[zgrid, P] = Tauchen(Nz, 0, rhoz, sigma_z, m);   % log productivity shocks
zgrid = exp(zgrid);                              % convert to levels
Z_min = min(zgrid);

%% ----------------------------------------------------------
% 2) CAPITAL GRID 
% -----------------------------------------------------------
Nk = 100;
k_min = 0.01;
k_max = 5000;
kgrid = linspace(k_min, k_max, Nk)';

%% ----------------------------------------------------------
% 3) UTILIZATION GRID 
% -----------------------------------------------------------
Nu = 15;                           
u_min = 0.5; 
u_max = 1.5;
ugrid = linspace(u_min, u_max, Nu)';

%% ----------------------------------------------------------
% 4) PRODUCTION FUNCTION AND DEPRECIATION LAW
%     y = z * (u*k)^alpha,  δ(u) = δ0 * u^θ
% -----------------------------------------------------------
prod = @(z,u,k) z .* (u .* k).^alpha; % Output when capital is used at rate u
delta = @(u) delta0 * (u.^theta);     % Depreciation rate depending on u

%% ----------------------------------------------------------
% 5) UTILITY FUNCTION
%     u(c) = log(c)
% -----------------------------------------------------------
util = @(c) log(c);                   % Utility from consumption
u_neg = -1e10;                        % Numerical stand-in for −∞ 

%% ----------------------------------------------------------
% 6) PRECOMPUTE UTILITY  U(k, k', z)
% -----------------------------------------------------------
U = -inf(Nk, Nk, Nz);

for iz = 1:Nz
    for iu = 1:Nu
        % resources available this period
        res = prod(zgrid(iz), ugrid(iu), kgrid) ...
              + (1 - delta(ugrid(iu))) * kgrid;     % output + remaining capital
        cons = res - kgrid';                        % implied consumption
        good = cons > 0;                            % Only Taking the Positive Consumptions
        tmp = -inf(Nk,Nk);                          % Creating a temporary matrix for capital choices
        tmp(good) = util(cons(good));               % log(c)
        % taking the *max* over utilization choices (inner FOC solved on grid)
        U(:,:,iz) = max(U(:,:,iz), tmp);            % Estimates utilization that maximizes current utility for each capital transition
    end
end

%% ----------------------------------------------------------
% 7) VALUE FUNCTION ITERATION  
% -----------------------------------------------------------
V  = zeros(Nk, Nz);
tol = 1e-9;
max_iter = 1000;
betaP = beta * P';

for iter = 1:max_iter
    EV = V * betaP;                            % expected continuation value
    V_new = zeros(Nk,Nz);
    polKix = zeros(Nk,Nz);
    for iz = 1:Nz
        [V_new(:,iz), polKix(:,iz)] = max(U(:,:,iz) + EV(:,iz)', [], 2);
    end
    if max(abs(V_new(:) - V(:))) < tol
        V = V_new; break;
    end
    V = V_new;
end

%% ----------------------------------------------------------
% 8) POLICY FUNCTIONS
% -----------------------------------------------------------
polK = kgrid(polKix);                          % optimal next capital
polC = zeros(Nk,Nz);
for iz = 1:Nz
    polC(:,iz) = prod(zgrid(iz),1,kgrid) ...
               + (1 - delta(1))*kgrid - polK(:,iz); %% crude u≈1 for reporting
end

%% ----------------------------------------------------------
% 9) SIMULATION  
% -----------------------------------------------------------
sims = 1000;
z_sim = simulate(dtmc(P), sims);
k_index = round(Nk/2) * ones(sims+1,1);
c_sim = zeros(sims,1);
k_sim = zeros(sims+1,1);

for t = 1:sims
    iz = z_sim(t);
    c_sim(t) = prod(zgrid(iz),1,kgrid(k_index(t))) ...
              + (1 - delta(1))*kgrid(k_index(t)) ...
              - kgrid(polKix(k_index(t),iz));
    k_index(t+1) = polKix(k_index(t),iz);
    k_sim(t+1) = kgrid(k_index(t+1));
end

% Plots
figure;
subplot(3,1,1)
plot(sims/2+1:sims, zgrid(z_sim(sims/2+1:sims)))
xlabel('Time'); ylabel('Productivity'); title('z_t')

subplot(3,1,2)
plot(sims/2+1:sims, c_sim(sims/2+1:sims))
xlabel('Time'); ylabel('Consumption')

subplot(3,1,3)
plot(sims/2+1:sims, k_sim(sims/2+1:sims))
xlabel('Time'); ylabel('Capital')

fprintf('Simulation complete. Iterations = %d\n', iter);

%% ----------------------------------------------------------
% (d)  Quantitative summary and qualitative discussion
% -----------------------------------------------------------
burn = sims/2;                      % Dropping the first half simulations
% Here I only keep the second half of each simulated time series
n_trim    = n_sim(burn+1:end);
u_trim    = u_sim(burn+1:end);
invp_trim = invp_sim(burn+1:end);

% Computing the sample standard deviation of each trimmed series
std_n    = std(n_trim);
std_u    = std(u_trim);
std_invp = std(invp_trim);

fprintf('Std(n)   = %.4f\n', std_n);
fprintf('Std(u)   = %.4f\n', std_u);
fprintf('Std(inv′)= %.4f\n', std_invp);

% Output 
Std(n)   = 0.0231;
Std(u)   = 0.0125;
Std(invp)= 0.0879; 


%% ----------------------------------------------------------
% If  (a) ϕ2 doubled.
%% ----------------------------------------------------------
% Higher ϕ₂ means steeper cost curve.Using capital more intensely becomes very expensive.
% Firms will smooth u much more, so std(u) falls. 

% Because u responds less to productivity shocks, the firm also invests less aggressively
% thus std(inv′) falls.

% Lower u volatility reduces fluctuations in effective capital services (uk), which dampens labor demand
% so std(n) also falls 


%% ----------------------------------------------------------
% If  (b) The real interest rate (r) doubled.
%% ----------------------------------------------------------
% A higher r means future profits are discounted more heavily—firms become less patient.
% The steady-state capital stock k shrinks.Capital adjustment (investment) becomes costlier
% so std(inv′) falls in relative terms, but the economy's overall capital volatility can rise in percentage terms depending on how tight the borrowing constraint is.

% Because there's less capital, the marginal product of labor fluctuates more for a given shock; 
% std(n) can rise modestly.

% Utilization u typically becomes more responsive: with less installed k, 
% firms use what they have more intensively during booms, leading std(u) to rise 