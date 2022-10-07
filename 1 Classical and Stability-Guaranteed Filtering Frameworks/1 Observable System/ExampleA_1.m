%%% Example A - 1
%
%   (c) Yirui Cong, created: 01-Oct-2021, last modified: --

close all;
clear;
clc;

rng(1); % For reproducible results


%%  Parameters
kmax = 1000;
T = 1;

A = [1 T; 0 1];
B = [T^2/2; T];
C = [1 0];

n = size(A, 1);
p = size(B, 2);
m = size(C, 1);

%%-  Observability Decomposition
[A_o, B_o, C_o, obsv_flag, P, A_obar, A_21, B_obar, n_o] = obsv_dec(A, B, C);
% obsv_flag: 0 - not detectable, 1 - detectable but not observable, 2 - observable

mu_o = observ_index(A_o, C_o);

n_obar = n - n_o;

r_B_o = rank(B_o);

if obsv_flag == 0
    disp('The pair (A, C) is not detectable!');
elseif obsv_flag == 1
    disp('The pair (A, C) is not observable but detectable.');
    %- Test
    test_error_max = max(max(abs(P \ [A_o zeros(n_o, n_obar); A_21 A_obar] * P - A)));
else
    disp('The pair (A, C) is observable.');
    A_o = A;
    B_o = B;
    C_o = C;
    P = eye(n);
end

%   True initial range
x0RangePara = 1;
x0CenterPara = 2 * ones(n, 1);

%   Initial condition A
x0RangeParaA = 2;
x0CenterParaA = 2 * ones(n, 1);

%   Initial condition B
x0RangeParaB = 4;
x0CenterParaB = zeros(n, 1);

%   Initial condition C (for OIT-inspired filter in Algorithm 2)
x0RangeParaC = 1;
x0CenterParaC = zeros(n, 1);

wkRangePara = 1;
wkCenterPara = zeros(p, 1);
vkRangePara = 1;
vkCenterPara = zeros(m, 1);

delta = 2; % The observation-horizon length is delta + 1.


%%  Simulation
%   True ranges
G_k_prior_total = cell(kmax+1, 1);
theta_k_prior_total = cell(kmax+1, 1);
G_k_posterior_total = cell(kmax+1, 1);
theta_k_posterior_total = cell(kmax+1, 1);

%   Set-membership filter A
G_k_prior_total_A = cell(kmax+1, 1);
theta_k_prior_total_A = cell(kmax+1, 1);
G_k_posterior_total_A = cell(kmax+1, 1);
theta_k_posterior_total_A = cell(kmax+1, 1);

%   Set-membership filter B
G_k_prior_total_B = cell(kmax+1, 1);
theta_k_prior_total_B = cell(kmax+1, 1);
G_k_posterior_total_B = cell(kmax+1, 1);
theta_k_posterior_total_B = cell(kmax+1, 1);

%   Set-membership filter C
G_k_prior_total_C = cell(kmax+1, 1);
theta_k_prior_total_C = cell(kmax+1, 1);
G_k_posterior_total_C = cell(kmax+1, 1);
theta_k_posterior_total_C = cell(kmax+1, 1);

G_k_OIT_total = cell(kmax+1, 1);
theta_k_OIT_total = cell(kmax+1, 1);

kIndexC = 1; % A compensator for the index 0 in matlab: for example, y[0] in matlab is y(1) = y(0 + k_indexC).
kSequence = 0: kmax;

xSequence = zeros(n, kmax+kIndexC);
ySequence = zeros(m, kmax+kIndexC);
wSequence = wkRangePara * 2 * (rand(p, kmax+kIndexC) - 0.5) + wkCenterPara * ones(1, kmax+kIndexC);
vSequence = vkRangePara * 2 * (rand(m, kmax+kIndexC) - 0.5) + vkCenterPara * ones(1, kmax+kIndexC);

G_w = kron(eye(p), [1; -1]);
theta_w = wkRangePara * ones(2*p, 1) + kron(wkCenterPara, [1; -1]);
G_v = kron(eye(m), [1; -1]);
theta_v = vkRangePara * ones(2*m, 1) + kron(vkCenterPara, [1; -1]);

counter = zeros(2, 1);

for k = kSequence
    k
    %%  Realizations of States and Measurements
    if k == 0
        xSequence(:, kIndexC) = x0RangePara * 2 * (rand(n, 1) - 0.5) + x0CenterPara;
    else
        xSequence(:, k+kIndexC) = A * xSequence(:, k-1+kIndexC) + B * wSequence(:, k-1+kIndexC);
    end
    
    ySequence(:, k+kIndexC) = C * xSequence(:, k+kIndexC) + vSequence(:, k+kIndexC);
    
    %%  Initialization & Prediction - Prior
    if k == 0
        %%  Initialization
        %   True prior range
        G_k_prior = kron(eye(n), [1; -1]);
        theta_k_prior = x0RangePara * ones(2*n, 1) + kron(x0CenterPara, [1; -1]);
        
        %   Prior range A
        G_k_prior_A = kron(eye(n), [1; -1]);
        theta_k_prior_A = x0RangeParaA * ones(2*n, 1) + kron(x0CenterParaA, [1; -1]);
        
        %   Prior range B
        G_k_prior_B = kron(eye(n), [1; -1]);
        theta_k_prior_B = x0RangeParaB * ones(2*n, 1) + kron(x0CenterParaB, [1; -1]);
        
        %   Prior range C
        G_k_prior_C = kron(eye(n), [1; -1]);
        theta_k_prior_C = x0RangeParaC * ones(2*n, 1) + kron(x0CenterParaC, [1; -1]);
        
        G_k_prior_total_C{k+kIndexC} = G_k_prior_C;
        theta_k_prior_total_C{k+kIndexC} = theta_k_prior_C;
    else
        %%  Prediction
        %   True prior range
        [G_k_prior, theta_k_prior] = OLSMF_prediction(A, B, G_k_posterior, theta_k_posterior, G_w, theta_w);
        
        %   Prior range A
        [G_k_prior_A, theta_k_prior_A] = OLSMF_prediction(A, B, G_k_posterior_A, theta_k_posterior_A, G_w, theta_w);
        
        %   Prior range B
        [G_k_prior_B, theta_k_prior_B] = OLSMF_prediction(A, B, G_k_posterior_B, theta_k_posterior_B, G_w, theta_w);
        
        %   Prior range C
        if k < mu_o - 1
            [G_k_prior_C, theta_k_prior_C] = OLSMF_prediction(A, B, G_k_posterior_C, theta_k_posterior_C, G_w, theta_w);
        end
    end
    
    %   True prior range
    G_k_prior_total{k+kIndexC} = G_k_prior;
    theta_k_prior_total{k+kIndexC} = theta_k_prior;
    
    %   Prior range A
    G_k_prior_total_A{k+kIndexC} = G_k_prior_A;
    theta_k_prior_total_A{k+kIndexC} = theta_k_prior_A;
    
    %   Prior range B
    G_k_prior_total_B{k+kIndexC} = G_k_prior_B;
    theta_k_prior_total_B{k+kIndexC} = theta_k_prior_B;
    
    %   Prior range C
    if k < mu_o - 1
        G_k_prior_total_C{k+kIndexC} = G_k_prior_C;
        theta_k_prior_total_C{k+kIndexC} = theta_k_prior_C;
    end
    
    %%  Update - Posterior
    %   True posterior range
    [G_k_posterior, theta_k_posterior] = OLSMF_update(C, ySequence(:, k+kIndexC), G_k_prior, theta_k_prior, G_v, theta_v);
    G_k_posterior_total{k+kIndexC} = G_k_posterior;
    theta_k_posterior_total{k+kIndexC} = theta_k_posterior;
    
    %   Posterior range A
    [G_k_posterior_A, theta_k_posterior_A] = OLSMF_update(C, ySequence(:, k+kIndexC), G_k_prior_A, theta_k_prior_A, G_v, theta_v);
    G_k_posterior_total_A{k+kIndexC} = G_k_posterior_A;
    theta_k_posterior_total_A{k+kIndexC} = theta_k_posterior_A;
    
    %   Posterior range B
    [G_k_posterior_B, theta_k_posterior_B] = OLSMF_update(C, ySequence(:, k+kIndexC), G_k_prior_B, theta_k_prior_B, G_v, theta_v);
    G_k_posterior_total_B{k+kIndexC} = G_k_posterior_B;
    theta_k_posterior_total_B{k+kIndexC} = theta_k_posterior_B;
    
    %   Posterior range C
    if k < mu_o - 1
        [G_k_posterior_C, theta_k_posterior_C] = OLSMF_update(C, ySequence(:, k+kIndexC), G_k_prior_C, theta_k_prior_C, G_v, theta_v);
        
        [x, fval, exitflag, output] = linprog(zeros(1, n), G_k_posterior_C, theta_k_posterior_C);
        alpha = x0RangeParaC; % Initial alpha, which can be tuned
        while exitflag == -2 % The posterior range is empty
            alpha = alpha * 3;
            
            G_k_prior_temp = kron(eye(n), [1; -1]);
            theta_k_prior_temp = alpha * ones(2*n, 1);
            [G_k_posterior_temp, theta_k_posterior_temp] = OLSMF_update(C, ySequence(:, 0+kIndexC), G_k_prior_temp, theta_k_prior_temp, G_v, theta_v);
            
            for i = 1: k
                [G_k_prior_temp, theta_k_prior_temp] = OLSMF_prediction(A, B, G_k_posterior_temp, theta_k_posterior_temp, G_w, theta_w);
                [G_k_posterior_temp, theta_k_posterior_temp] = OLSMF_update(C, ySequence(:, i+kIndexC), G_k_prior_temp, theta_k_prior_temp, G_v, theta_v);
            end
            
            [x, fval, exitflag, output] = linprog(zeros(1, n), G_k_posterior_temp, theta_k_posterior_temp);
            G_k_posterior_C = G_k_posterior_temp;
            theta_k_posterior_C = theta_k_posterior_temp;
        end
        
        G_k_posterior_total_C{k+kIndexC} = G_k_posterior_C;
        theta_k_posterior_total_C{k+kIndexC} = theta_k_posterior_C;
    elseif k == mu_o - 1
        G_k_prior_temp = kron(eye(n), [1; -1]);
        theta_k_prior_temp = 1000 * ones(2*n, 1); % Sufficiently large range for this example
        
        [G_k_posterior_temp, theta_k_posterior_temp] = OLSMF_update(C, ySequence(:, 0+kIndexC), G_k_prior_temp, theta_k_prior_temp, G_v, theta_v);
            
        for i = 1: k
            [G_k_prior_temp, theta_k_prior_temp] = OLSMF_prediction(A, B, G_k_posterior_temp, theta_k_posterior_temp, G_w, theta_w);
            [G_k_posterior_temp, theta_k_posterior_temp] = OLSMF_update(C, ySequence(:, i+kIndexC), G_k_prior_temp, theta_k_prior_temp, G_v, theta_v);
        end
        
        G_k_posterior_total_C{k+kIndexC} = G_k_posterior_temp;
        theta_k_posterior_total_C{k+kIndexC} = theta_k_posterior_temp;
    else
        [G_k_prior_C, theta_k_prior_C] = OLSMF_prediction(A, B, G_k_posterior_total_C{k-1+kIndexC}, theta_k_posterior_total_C{k-1+kIndexC}, G_w, theta_w);
        [G_k_posterior_C, theta_k_posterior_C] = OLSMF_update(C, ySequence(:, k+kIndexC), G_k_prior_C, theta_k_prior_C, G_v, theta_v);
        G_k_posterior_total_C{k+kIndexC} = G_k_posterior_C;
        theta_k_posterior_total_C{k+kIndexC} = theta_k_posterior_C;
    end
    
    
    %%  Reduced Observation-Information Tower
    if k >= delta && obsv_flag == 2
        [G_k_OIT_total{k+kIndexC}, theta_k_OIT_total{k+kIndexC}] = OIT(A, B, C, ySequence, k, delta, G_w, theta_w, G_v, theta_v, mu_o);
    end
end

%%  Calculations and Figures
%   True ranges
x_prior_volume = zeros(kmax+kIndexC, 1);
x_posterior_volume = zeros(kmax+kIndexC, 1);
x_posterior_diameter = zeros(kmax+kIndexC, 1);
x_OIT_volume = inf(kmax+kIndexC, 1);
x_OIT_diameter = inf(kmax+kIndexC, 1);

%   Set-membership filter A
x_prior_volume_A = zeros(kmax+kIndexC, 1);
x_posterior_volume_A = zeros(kmax+kIndexC, 1);
x_posterior_diameter_A = zeros(kmax+kIndexC, 1);
x_posterior_estimation_gap_A = zeros(kmax+kIndexC, 1);

%   Set-membership filter B
x_prior_volume_B = zeros(kmax+kIndexC, 1);
x_posterior_volume_B = zeros(kmax+kIndexC, 1);
x_posterior_diameter_B = zeros(kmax+kIndexC, 1);
x_posterior_estimation_gap_B = zeros(kmax+kIndexC, 1);

%   Set-membership filter C
x_prior_volume_C = zeros(kmax+kIndexC, 1);
x_posterior_volume_C = zeros(kmax+kIndexC, 1);
x_posterior_diameter_C = zeros(kmax+kIndexC, 1);
x_posterior_estimation_gap_C = zeros(kmax+kIndexC, 1);


%   Upper bound of diameters (k \geq delta_optimal)
diameter_upper_bound = inf;
d_w = 2 * wkRangePara * sqrt(p);
d_v = 2 * vkRangePara * sqrt(m);

delta_optimal = 0;

for delta_temp = mu_o - 1: 100
    sum_temp_j = 0;
    for j = 0: delta_temp
        sum_temp_l = 0;
        for l = 1: delta_temp - j
            sum_temp_l = sum_temp_l + norm(C / (A^l) * B) * d_w;
        end
        sum_temp_j = sum_temp_j + (d_v + sum_temp_l)^2;
    end
    
    O_delta = C;
    for j = 1: delta_temp
        O_delta = [C / (A^j); O_delta];
    end
    
    diameter_upper_bound_temp = sqrt(sum_temp_j) / min(svd(O_delta));
%     diameter_upper_bound_temp = norm(pinv(O_delta)) * sqrt(sum_temp_j);
    
    if diameter_upper_bound > diameter_upper_bound_temp
        diameter_upper_bound = diameter_upper_bound_temp;
        delta_optimal = delta_temp;
    end
    
%     diameter_upper_bound = min(diameter_upper_bound, diameter_upper_bound_temp);
end

for k = kSequence
    k
    
    %   True ranges
    [vertex_k_prior, nr_prior] = con2vert(G_k_prior_total{k+kIndexC}, theta_k_prior_total{k+kIndexC});
    [CH_k_prior, x_prior_volume(k+kIndexC)] = convhull(vertex_k_prior);
    
    [vertex_k_posterior, nr_posterior] = con2vert(G_k_posterior_total{k+kIndexC}, theta_k_posterior_total{k+kIndexC});
    [CH_k_posterior, x_posterior_volume(k+kIndexC)] = convhull(vertex_k_posterior);
    
    x_posterior_diameter(k+kIndexC) = diameter_conv(vertex_k_posterior);
    
    %   Set-membership filter A
    [vertex_k_prior_A, nr_prior_A] = con2vert(G_k_prior_total_A{k+kIndexC}, theta_k_prior_total_A{k+kIndexC});
    [CH_k_prior_A, x_prior_volume_A(k+kIndexC)] = convhull(vertex_k_prior_A);
    
    [vertex_k_posterior_A, nr_posterior_A] = con2vert(G_k_posterior_total_A{k+kIndexC}, theta_k_posterior_total_A{k+kIndexC});
    [CH_k_posterior_A, x_posterior_volume_A(k+kIndexC)] = convhull(vertex_k_posterior_A);
    
    x_posterior_diameter_A(k+kIndexC) = diameter_conv(vertex_k_posterior_A);
    x_posterior_estimation_gap_A(k+kIndexC) = estimation_gap(vertex_k_posterior_A, vertex_k_posterior);
    
    %   Set-membership filter B
    [vertex_k_prior_B, nr_prior_B] = con2vert(G_k_prior_total_B{k+kIndexC}, theta_k_prior_total_B{k+kIndexC});
    [CH_k_prior_B, x_prior_volume_B(k+kIndexC)] = convhull(vertex_k_prior_B);
    
    [vertex_k_posterior_B, nr_posterior_B] = con2vert(G_k_posterior_total_B{k+kIndexC}, theta_k_posterior_total_B{k+kIndexC});
    [CH_k_posterior_B, x_posterior_volume_B(k+kIndexC)] = convhull(vertex_k_posterior_B);
    
    x_posterior_diameter_B(k+kIndexC) = diameter_conv(vertex_k_posterior_B);
    x_posterior_estimation_gap_B(k+kIndexC) = estimation_gap(vertex_k_posterior_B, vertex_k_posterior);
    
    %   Set-membership filter C
    if k == 0
        [vertex_k_prior_C, nr_prior_C] = con2vert(G_k_prior_total_C{k+kIndexC}, theta_k_prior_total_C{k+kIndexC});
        [CH_k_prior_C, x_prior_volume_C(k+kIndexC)] = convhull(vertex_k_prior_C);
    end
    
    [vertex_k_posterior_C, nr_posterior_C] = con2vert(G_k_posterior_total_C{k+kIndexC}, theta_k_posterior_total_C{k+kIndexC});
    [CH_k_posterior_C, x_posterior_volume_C(k+kIndexC)] = convhull(vertex_k_posterior_C);
    
    x_posterior_diameter_C(k+kIndexC) = diameter_conv(vertex_k_posterior_C);
    x_posterior_estimation_gap_C(k+kIndexC) = estimation_gap(vertex_k_posterior_C, vertex_k_posterior);
    
    if k >= delta && obsv_flag == 2
        [vertex_k_OIT, nr_OIT] = con2vert(G_k_OIT_total{k+kIndexC}, theta_k_OIT_total{k+kIndexC});
        [CH_k_OIT, x_OIT_volume(k+kIndexC)] = convhull(vertex_k_OIT);
        x_OIT_diameter(k+kIndexC) = diameter_conv(vertex_k_OIT);
    end
    
%     if k >= kmax - 5
% %     if k == 500
%         figure,
%         if obsv_flag == 2
%             plot(vertex_k_prior(CH_k_prior, 1), vertex_k_prior(CH_k_prior, 2), vertex_k_posterior(CH_k_posterior, 1), vertex_k_posterior(CH_k_posterior, 2),...
%                vertex_k_OIT(CH_k_OIT, 1), vertex_k_OIT(CH_k_OIT, 2))
%             hold on;
%             plot(vertex_k_prior(:, 1), vertex_k_prior(:, 2), 's', vertex_k_posterior(:, 1), vertex_k_posterior(:, 2), 'o', vertex_k_OIT(:, 1), vertex_k_OIT(:, 2), '*')
%             legend('Prior Range', 'Posterior Range', 'OIT')
%             title(strcat('k = ', num2str(k)))
%         else
%             plot(vertex_k_prior(CH_k_prior, 1), vertex_k_prior(CH_k_prior, 2), vertex_k_posterior(CH_k_posterior, 1), vertex_k_posterior(CH_k_posterior, 2))
%             hold on;
%             plot(vertex_k_prior(:, 1), vertex_k_prior(:, 2), 's', vertex_k_posterior(:, 1), vertex_k_posterior(:, 2), 'o')
%             legend('Prior Range', 'Posterior Range')
%             title(strcat('k = ', num2str(k)))
%         end
%         grid on;
%     end
    
    if k == 0
        figure,
        plot(vertex_k_prior(CH_k_prior, 1), vertex_k_prior(CH_k_prior, 2), vertex_k_prior_A(CH_k_prior_A, 1), vertex_k_prior_A(CH_k_prior_A, 2),...
            vertex_k_prior_B(CH_k_prior_B, 1), vertex_k_prior_B(CH_k_prior_B, 2), vertex_k_prior_C(CH_k_prior_C, 1), vertex_k_prior_C(CH_k_prior_C, 2))
        legend('True prior range', 'Alice', 'Bob', 'Carol')
        grid on;
    end
    
    if k <= 10
        figure,
%         plot(vertex_k_posterior(CH_k_posterior, 1), vertex_k_posterior(CH_k_posterior, 2), vertex_k_posterior_A(CH_k_posterior_A, 1), vertex_k_posterior_A(CH_k_posterior_A, 2),...
%             vertex_k_posterior_B(CH_k_posterior_B, 1), vertex_k_posterior_B(CH_k_posterior_B, 2), vertex_k_posterior_C(CH_k_posterior_C, 1), vertex_k_posterior_C(CH_k_posterior_C, 2),...
%             vertex_k_OIT(CH_k_OIT, 1), vertex_k_OIT(CH_k_OIT, 2))
%         legend('True posterior range', 'Alice', 'Bob', 'Carol', 'OIT')
        plot(vertex_k_posterior(CH_k_posterior, 1), vertex_k_posterior(CH_k_posterior, 2), vertex_k_posterior_A(CH_k_posterior_A, 1), vertex_k_posterior_A(CH_k_posterior_A, 2),...
            vertex_k_posterior_B(CH_k_posterior_B, 1), vertex_k_posterior_B(CH_k_posterior_B, 2), vertex_k_posterior_C(CH_k_posterior_C, 1), vertex_k_posterior_C(CH_k_posterior_C, 2))
        legend('True posterior range', 'Alice', 'Bob', 'Carol')
        grid on;
    end
end

figure,
% plot(kSequence, x_prior_volume, kSequence, x_posterior_volume)
% plot(kSequence, x_prior_volume, kSequence, x_posterior_volume, kSequence, x_OIT_volume)
if obsv_flag == 2
    plot(kSequence, x_posterior_volume, kSequence, x_OIT_volume)
    legend('Posterior range', 'OIT')
else
    plot(kSequence, x_posterior_volume)
    legend('Posterior range')
end
title('Volume')
grid on;

figure,
if obsv_flag == 2
    plot(kSequence, x_posterior_diameter, kSequence, x_OIT_diameter)
    legend('Posterior range', 'OIT')
else
    plot(kSequence, x_posterior_diameter)
    legend('Posterior range')
end
title('Diamter')
grid on;

figure,
if obsv_flag == 2
    plot(kSequence, x_posterior_volume, kSequence, x_posterior_volume_A, kSequence, x_posterior_volume_B, kSequence, x_posterior_volume_C, kSequence, x_OIT_volume)
    legend('True posterior range', 'Alice', 'Bob', 'Carol', 'OIT')
else
    plot(kSequence, x_posterior_volume, kSequence, x_posterior_volume_A, kSequence, x_posterior_volume_B, kSequence, x_posterior_volume_C)
    legend('True posterior range', 'Alice', 'Bob', 'Carol')
end
title('Volume')
grid on;

figure,
if obsv_flag == 2
    plot(kSequence, x_posterior_diameter, kSequence, x_posterior_diameter_A, kSequence, x_posterior_diameter_B, kSequence, x_posterior_diameter_C, kSequence, x_OIT_diameter, [min(kSequence) max(kSequence)], [diameter_upper_bound, diameter_upper_bound])
    legend('True posterior range', 'Alice', 'Bob', 'Carol', 'OIT')
else
    plot(kSequence, x_posterior_diameter, kSequence, x_posterior_diameter_A, kSequence, x_posterior_diameter_B, kSequence, x_posterior_diameter_C)
    legend('True posterior range', 'Alice', 'Bob', 'Carol')
end
title('Diamter')
grid on;

diameter_upper_bound

% diameter_upper_bound_new = bound_diameter_estimate(A, B, C, 2, d_w, d_v)

figure,
plot(kSequence, x_posterior_estimation_gap_A, kSequence, x_posterior_estimation_gap_B, kSequence, x_posterior_estimation_gap_C)
legend('Alice', 'Bob', 'Carol')
title('Estimation Gap')
grid on;