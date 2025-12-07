%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problem 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% Parameters (baseline-style)
beta   = 0.9;
delta  = 0.10;
F      = 0.4;     % fixed replacement cost (slightly higher than given, I did this for getting a smoother behavior)
A      = 1.0;     % normalizing aggregate productivity (no aggregate shocks)
eps_n  = 11;
eps_min = 0.4;
eps_max = 1.6;
eps_grid = linspace(eps_min, eps_max, eps_n);
pi_eps   = (1/eps_n) * ones(1,eps_n);           % iid uniform

%% Capital grid (vintages)
k_n   = 40;
k_max = 1;
k_vec = k_max * (1-delta).^(0:k_n-1)';  

%% Value function iteration (no aggregate shocks, states = (k, eps))
V      = zeros(k_n, eps_n);
V_new  = V;
z      = zeros(k_n, eps_n);  % 1 = replace, 0 = keep

tol    = 1e-6;
maxits = 5000;
diff   = 1;
count  = 0;

% Precompute expectation over eps'
EV_rep = zeros(1,eps_n);     
while diff > tol && count < maxits

    % given current V, computing expected continuation values
    EV_next = V;  % continuation value depends on (k',eps) directly


    % value of starting with new capital next period
    EV_rep  = V(1,:) * pi_eps'; 
    for ik = 1:k_n
        for ie = 1:eps_n

            k  = k_vec(ik);
            ep = eps_grid(ie);
            y  = A * ep * k;

            % ---- Option 1: Keep machine ----
            % next capital = (1-delta)*k, map to nearest grid point
            k_next = (1-delta) * k;
            [~, ik_next] = min(abs(k_vec - k_next));
            V_keep = y + beta * EV_next(ik_next, ie);

            % ---- Option 2: Replace machine ----
            k_rep = k_max;
            profit_rep = A * ep * k_rep - F;
            V_rep = profit_rep + beta * EV_rep;

            % choose best
            if V_rep >= V_keep
                V_new(ik,ie) = V_rep;
                z(ik,ie)     = 1;
            else
                V_new(ik,ie) = V_keep;
                z(ik,ie)     = 0;
            end

        end
    end

    diff   = max(max(abs(V_new - V)));
    V      = V_new;
    count  = count + 1;
end

fprintf('VFI converged in %d iterations, diff = %.2e\n',count,diff);

%% Check policy: Does replacement probability rise with age?
% (not strictly necessary, but good sanity check)
hazard_age = zeros(k_n,1);
for ik = 1:k_n
    % fraction of eps for which replace=1
    hazard_age(ik) = mean(z(ik,:));
end
disp('First 10 age-hazards:');
disp(hazard_age(1:10)');

%% ========= Simulation for Convergence Without Aggregate Shocks =========%%

T   = 50;       % number of periods
Nf  = 2000;     % number of firms

% Random initial machine age for cross-section of firms
k_index = randi([1 k_n], Nf, 1);    

inv_rate = zeros(T,1);

% draw eps shocks (iid uniform)
eps_sim = randsample(1:eps_n, T*Nf, true, pi_eps);
eps_sim = reshape(eps_sim, T, Nf);

for t = 1:T
    invest_count = 0;

    for i = 1:Nf
        
        % policy decision
        action = z(k_index(i), eps_sim(t,i));

        % if replace = 1
        if action == 1
            invest_count = invest_count + 1;
            k_index(i) = 1;      % reset capital
        else
            % keep machine, go to next lower capital index
            k_index(i) = min(k_index(i)+1, k_n);
        end
    end

    inv_rate(t) = invest_count / Nf;
end

%% ======================  Plot Result =======================
figure;
plot(inv_rate,'LineWidth',2)
ylim([0 0.8])
title('Convergence Without Aggregate Shocks (Problem 7)')
xlabel('Period')
ylabel('Investment Rate')
grid on
