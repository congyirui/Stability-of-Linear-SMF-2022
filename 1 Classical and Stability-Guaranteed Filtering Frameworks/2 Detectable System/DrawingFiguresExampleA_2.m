%%% Drawing Figures

close all;
clear;
clc;

load('experiment20211002')


%%  Ranges
kSequence_selected = [0 1 2 3 6 10];

for k = kSequence_selected
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
    
    if k == 0
        figure,
        plot(vertex_k_prior(CH_k_prior, 1), vertex_k_prior(CH_k_prior, 2), '-s', vertex_k_prior_A(CH_k_prior_A, 1), vertex_k_prior_A(CH_k_prior_A, 2), '-o',...
            vertex_k_prior_B(CH_k_prior_B, 1), vertex_k_prior_B(CH_k_prior_B, 2), '-*', vertex_k_prior_C(CH_k_prior_C, 1), vertex_k_prior_C(CH_k_prior_C, 2), '-^', 'LineWidth', 1.2, 'MarkerSize', 8)
        legend('True prior', 'Alice', 'Bob', 'Carol')
        grid on;
        set(gca,'FontSize',12);
        xlabel('x^{(1)}-axis')
        ylabel('x^{(2)}-axis')
    end
    
    if k > 0
        if k >= delta
            figure,
            plot(vertex_k_posterior(CH_k_posterior, 1), vertex_k_posterior(CH_k_posterior, 2), '-s', vertex_k_posterior_A(CH_k_posterior_A, 1), vertex_k_posterior_A(CH_k_posterior_A, 2), '-o',...
                vertex_k_posterior_B(CH_k_posterior_B, 1), vertex_k_posterior_B(CH_k_posterior_B, 2), '-*', vertex_k_posterior_C(CH_k_posterior_C, 1), vertex_k_posterior_C(CH_k_posterior_C, 2), '-^',...
                'LineWidth', 1.2, 'MarkerSize', 8)
            legend('True conditional range', 'Alice', 'Bob', 'Carol')
            grid on;
            set(gca,'FontSize',12);
            xlabel('x^{(1)}-axis')
            ylabel('x^{(2)}-axis')
        else
            figure,
            plot(vertex_k_posterior(CH_k_posterior, 1), vertex_k_posterior(CH_k_posterior, 2), '-s', vertex_k_posterior_A(CH_k_posterior_A, 1), vertex_k_posterior_A(CH_k_posterior_A, 2), '-o',...
                vertex_k_posterior_B(CH_k_posterior_B, 1), vertex_k_posterior_B(CH_k_posterior_B, 2), '-*', vertex_k_posterior_C(CH_k_posterior_C, 1), vertex_k_posterior_C(CH_k_posterior_C, 2), '-^',...
                'LineWidth', 1.2, 'MarkerSize', 8)
            legend('True conditional range', 'Alice', 'Bob', 'Carol')
            grid on;
            set(gca,'FontSize',12);
            xlabel('x^{(1)}-axis')
            ylabel('x^{(2)}-axis')
        end
    end
end


%%  Maximum Diameters
MD = max(x_posterior_diameter(delta_optimal+kIndexC:end))
MD_A = max(x_posterior_diameter_A(delta_optimal+kIndexC:end))
MD_B = max(x_posterior_diameter_B(delta_optimal+kIndexC:end))
MD_C = max(x_posterior_diameter_C(delta_optimal+kIndexC:end))
diameter_upper_bound

k_draw_max = 50+kIndexC;

figure,
plot(kSequence(1:k_draw_max), x_posterior_diameter(1:k_draw_max), '-s', kSequence(1:k_draw_max), x_posterior_diameter_A(1:k_draw_max), '-o', kSequence(1:k_draw_max), x_posterior_diameter_B(1:k_draw_max), '-*', kSequence(1:k_draw_max), x_posterior_diameter_C(1:k_draw_max), '-^', [min(kSequence(1:k_draw_max)) max(kSequence(1:k_draw_max))], [diameter_upper_bound, diameter_upper_bound], '-', 'LineWidth', 1.2, 'MarkerSize', 6)
% hold on;
% plot([1 1], [0 9], '--')
legend('True conditional range', 'Alice', 'Bob', 'Carol', 'Upper Bound')
grid on;
set(gca,'FontSize',10);
ylim([0 9])
xlabel('Time Step')
ylabel('Diameter')


%%  Estimation Gap
figure,
semilogy(kSequence(1:k_draw_max), x_posterior_estimation_gap_A(1:k_draw_max), '-o', 'Color', [0.85,0.33,0.10], 'LineWidth', 1.2, 'MarkerSize', 6)
hold on;
semilogy(kSequence(1:k_draw_max), x_posterior_estimation_gap_B(1:k_draw_max), '-*', 'Color', [0.93,0.69,0.13], 'LineWidth', 1.2, 'MarkerSize', 6)
semilogy(kSequence(1:k_draw_max), x_posterior_estimation_gap_C(1:k_draw_max), '-^', 'Color', [0.49,0.18,0.56], 'LineWidth', 1.2, 'MarkerSize', 6)
legend('Alice', 'Bob', 'Carol')
title('Estimation Gap')
grid on;
set(gca,'FontSize',10);
xlabel('Time Step')
ylabel('Estimation Gap')


%%  Figure 3
fid = figure;
set(fid,'position',[277.8,161,868,564.8000000000001]);

subplot(2, 2, [1 3])
% plot(kSequence(1:k_draw_max), x_posterior_diameter(1:k_draw_max), '-s', kSequence(1:k_draw_max), x_posterior_diameter_A(1:k_draw_max), '-o', kSequence(1:k_draw_max), x_posterior_diameter_B(1:k_draw_max), '-*', kSequence(1:k_draw_max), x_posterior_diameter_C(1:k_draw_max), '-^', [min(kSequence(1:k_draw_max)) max(kSequence(1:k_draw_max))], [diameter_upper_bound, diameter_upper_bound], '-', 'LineWidth', 1.2, 'MarkerSize', 6)
plot(kSequence(1:k_draw_max), x_posterior_diameter(1:k_draw_max), '-s', kSequence(1:k_draw_max), x_posterior_diameter_A(1:k_draw_max), '-o', kSequence(1:k_draw_max), x_posterior_diameter_B(1:k_draw_max), '-*', kSequence(1:k_draw_max), x_posterior_diameter_C(1:k_draw_max), '-^', 'LineWidth', 1.2, 'MarkerSize', 6)
% hold on;
% plot([1 1], [0 9], '--')
% legend('Actual posterior', 'Alice', 'Bob', 'Carol', 'Upper Bound')
legend('True', 'Alice', 'Bob', 'Carol')
grid on;
set(gca,'FontSize',16);
xlim([0 50])
ylim([0 9])
xlabel('Time Step')
ylabel('Diameter')

subplot(2, 2, [2 4])
semilogy(kSequence(1:k_draw_max), x_posterior_estimation_gap_A(1:k_draw_max), '-o', 'Color', [0.85,0.33,0.10], 'LineWidth', 1.2, 'MarkerSize', 6)
hold on;
semilogy(kSequence(1:k_draw_max), x_posterior_estimation_gap_B(1:k_draw_max), '-*', 'Color', [0.93,0.69,0.13], 'LineWidth', 1.2, 'MarkerSize', 6)
semilogy(kSequence(1:k_draw_max), x_posterior_estimation_gap_C(1:k_draw_max), '-^', 'Color', [0.49,0.18,0.56], 'LineWidth', 1.2, 'MarkerSize', 6)
legend('Alice', 'Bob', 'Carol')
grid on;
set(gca,'FontSize',16);
xlim([0 50])
xlabel('Time Step')
ylabel('Estimation Gap')