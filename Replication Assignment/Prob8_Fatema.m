%%%%%%%%%%%%%%%%%%%% Problem 8 – With Aggregate Shocks %%%%%%%%%%%%%%%%%%%%

clear; clc;

%% 1. Parameters (calibrated for clear fluctuation like Cooper, Haltiwanger & Power (1999) Figure 4)

beta  = 0.92;
delta = 0.10;         
F     = 0.30;          % Higher fixed cost → boom/bust replacement cycles

% Aggregate productivities
A_states = [0.80, 1.20];     % strong recession vs strong boom
A_n = length(A_states);

% Transition matrix 
P_A = [0.90 0.10;
       0.10 0.90];

% Idiosyncratic shocks
eps_n  = 11;
eps_grid = linspace(0.4,1.6,eps_n);
pi_eps = (1/eps_n) * ones(1,eps_n);   % i.i.d uniform

% Capital grid
k_n = 40;
k_vec = (1-delta).^(0:k_n-1)';

%% 2. Value Function Iteration with Aggregate Shocks

V = zeros(k_n,eps_n,A_n);
V_new = V;
z = zeros(k_n,eps_n,A_n);   % 1=keep, 2=replace

tol=1e-6; diff=1; it=0;

while diff>tol && it<3000
    for a=1:A_n
        for ik=1:k_n
            for ie=1:eps_n

                k = k_vec(ik);
                ep = eps_grid(ie);
                A=A_states(a);

                y = A*ep*k;

                % KEEP
                k_next = (1-delta)*k;
                [~,ik2]=min(abs(k_vec-k_next));

                EV_keep = beta*( P_A(a,1)*sum(pi_eps.*V(ik2,:,1)) + ...
                                 P_A(a,2)*sum(pi_eps.*V(ik2,:,2)) );

                % REPLACE
                profit_rep = A*ep*1 - F;
                EV_rep = beta*( P_A(a,1)*sum(pi_eps.*V(1,:,1)) + ...
                                P_A(a,2)*sum(pi_eps.*V(1,:,2)) );

                [V_new(ik,ie,a),choice] = max([EV_keep+y , EV_rep]);
                z(ik,ie,a)=choice;
            end
        end
    end

    diff=max(abs(V_new(:)-V(:)));
    V=V_new;
    it=it+1;
end

disp("VFI Converged.")
disp("Replacement rates in policy function:")
disp("LOW A  = " + mean(z(:,:,1)==2,'all'))
disp("HIGH A = " + mean(z(:,:,2)==2,'all'))

%% 3. Simulation

T = 150;
Nf = 5000;

% Simulate aggregate shocks
A_index = zeros(T,1);
A_index(1)=1;
for t=2:T
    A_index(t) = find(rand < cumsum(P_A(A_index(t-1),:)),1);
end

% Idiosyncratic shocks
eps_sim = randi(eps_n,T,Nf);

% Firms
k_index = ones(Nf,1);
inv_rate = zeros(T,1);

for t=1:T
    invest_count=0;
    Ast=A_index(t);

    for i=1:Nf
        ie=eps_sim(t,i);

        if z(k_index(i),ie,Ast)==2   % Replace decision
            invest_count=invest_count+1;
            k_index(i)=1;
        else
            k_next=(1-delta)*k_vec(k_index(i));
            [~,kk]=min(abs(k_vec-k_next));
            k_index(i)=kk;
        end
    end
    inv_rate(t)=invest_count/Nf;
end

%% 4. Plot — comparable to Cooper, Haltiwanger & Power (1999) Figure 4

figure;
yyaxis left
plot(inv_rate,'LineWidth',1.4)
ylim([0.15 0.45])
ylabel('Investment Rate')

yyaxis right
stairs(A_index,'--','LineWidth',1.2)
ylim([0 2])
yticks([1 2])
yticklabels({'LOW','HIGH'})
ylabel('Aggregate State')

title('Aggregate Fluctuations (Replicating Fig. 4)')
xlabel('Time')
legend('Investment Rate','Aggregate State')
grid on
