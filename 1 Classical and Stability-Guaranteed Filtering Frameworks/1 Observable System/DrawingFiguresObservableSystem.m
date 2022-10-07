%%% Drawing Figures

close all;
clear;
clc;

load('experiment20200211')


%%  Figure 3
k = 500;

[vertex_k_prior, nr_prior] = con2vert(G_k_prior_total{k+kIndexC}, theta_k_prior_total{k+kIndexC});
[CH_k_prior, x_prior_volume(k+kIndexC)] = convhull(vertex_k_prior);

[vertex_k_posterior, nr_posterior] = con2vert(G_k_posterior_total{k+kIndexC}, theta_k_posterior_total{k+kIndexC});
[CH_k_posterior, x_posterior_volume(k+kIndexC)] = convhull(vertex_k_posterior);

[vertex_k_ROIT, nr_ROIT] = con2vert(G_k_ROIT_total{k+kIndexC}, theta_k_ROIT_total{k+kIndexC});
[CH_k_ROIT, x_ROIT_volume(k+kIndexC)] = convhull(vertex_k_ROIT);
x_ROIT_diameter(k+kIndexC) = diameter_conv(vertex_k_ROIT);

figure,
plot(vertex_k_prior(CH_k_prior, 1), vertex_k_prior(CH_k_prior, 2), '-s',...
    vertex_k_posterior(CH_k_posterior, 1), vertex_k_posterior(CH_k_posterior, 2), '-o',...
    vertex_k_ROIT(CH_k_ROIT, 1), vertex_k_ROIT(CH_k_ROIT, 2), '-*', 'LineWidth', 1)
legend('Prior', 'Posterior', 'ROIT')
% title(strcat('k = ', num2str(k)))
grid on;


%%  Table 1
% kSequence_selected = [0 100 200 500 1000];
% x_posterior_diameter_selected = x_posterior_diameter(kSequence_selected + kIndexC)
% x_ROIT_volume_selected = x_ROIT_volume(kSequence_selected + kIndexC)
% x_posterior_volume_selected = x_posterior_volume(kSequence_selected + kIndexC)
% x_ROIT_diameter_selected = x_ROIT_diameter(kSequence_selected + kIndexC)

x_posterior_diameter_averaged = [sum(x_posterior_diameter(1: 251)) / 251; sum(x_posterior_diameter(252: 501)) / 250; sum(x_posterior_diameter(502: 751)) / 250; ...
    sum(x_posterior_diameter(752: 1001)) / 250]
x_ROIT_diameter_averaged = [sum(x_ROIT_diameter(1: 251)) / 251; sum(x_ROIT_diameter(252: 501)) / 250; sum(x_ROIT_diameter(502: 751)) / 250; ...
    sum(x_ROIT_diameter(752: 1001)) / 250]
x_posterior_volume_averaged = [sum(x_posterior_volume(1: 251)) / 251; sum(x_posterior_volume(252: 501)) / 250; sum(x_posterior_volume(502: 751)) / 250; ...
    sum(x_posterior_volume(752: 1001)) / 250]
x_ROIT_volume_averaged = [sum(x_ROIT_volume(1: 251)) / 251; sum(x_ROIT_volume(252: 501)) / 250; sum(x_ROIT_volume(502: 751)) / 250; ...
    sum(x_ROIT_volume(752: 1001)) / 250]

%%  Figure 4
kSequence_selected = [0 3 6];

for k = kSequence_selected
    k
    
    %   True initial prior range
    [vertex_k_prior, nr_prior] = con2vert(G_k_prior_total{k+kIndexC}, theta_k_prior_total{k+kIndexC});
    [CH_k_prior, x_prior_volume(k+kIndexC)] = convhull(vertex_k_prior);
    
    [vertex_k_posterior, nr_posterior] = con2vert(G_k_posterior_total{k+kIndexC}, theta_k_posterior_total{k+kIndexC});
    [CH_k_posterior, x_posterior_volume(k+kIndexC)] = convhull(vertex_k_posterior);
    
    %   Initial prior range A
    [vertex_k_prior_A, nr_prior_A] = con2vert(G_k_prior_total_A{k+kIndexC}, theta_k_prior_total_A{k+kIndexC});
    [CH_k_prior_A, x_prior_volume_A(k+kIndexC)] = convhull(vertex_k_prior_A);
    
    [vertex_k_posterior_A, nr_posterior_A] = con2vert(G_k_posterior_total_A{k+kIndexC}, theta_k_posterior_total_A{k+kIndexC});
    [CH_k_posterior_A, x_posterior_volume_A(k+kIndexC)] = convhull(vertex_k_posterior_A);
    
    %   Initial prior range B
    [vertex_k_prior_B, nr_prior_B] = con2vert(G_k_prior_total_B{k+kIndexC}, theta_k_prior_total_B{k+kIndexC});
    [CH_k_prior_B, x_prior_volume_B(k+kIndexC)] = convhull(vertex_k_prior_B);
    
    [vertex_k_posterior_B, nr_posterior_B] = con2vert(G_k_posterior_total_B{k+kIndexC}, theta_k_posterior_total_B{k+kIndexC});
    [CH_k_posterior_B, x_posterior_volume_B(k+kIndexC)] = convhull(vertex_k_posterior_B);
    
    %   Initial prior range C
    [vertex_k_prior_C, nr_prior_C] = con2vert(G_k_prior_total_C{k+kIndexC}, theta_k_prior_total_C{k+kIndexC});
    [CH_k_prior_C, x_prior_volume_C(k+kIndexC)] = convhull(vertex_k_prior_C);
    
    [vertex_k_posterior_C, nr_posterior_C] = con2vert(G_k_posterior_total_C{k+kIndexC}, theta_k_posterior_total_C{k+kIndexC});
    [CH_k_posterior_C, x_posterior_volume_C(k+kIndexC)] = convhull(vertex_k_posterior_C);
    
    if k == 0
        figure,
        plot(vertex_k_prior(CH_k_prior, 1), vertex_k_prior(CH_k_prior, 2), '-s', vertex_k_prior_A(CH_k_prior_A, 1), vertex_k_prior_A(CH_k_prior_A, 2), '-o',...
            vertex_k_prior_B(CH_k_prior_B, 1), vertex_k_prior_B(CH_k_prior_B, 2), '-*', vertex_k_prior_C(CH_k_prior_C, 1), vertex_k_prior_C(CH_k_prior_C, 2), '-^', 'LineWidth', 1.2, 'MarkerSize', 8)
        legend('True prior', 'Alice', 'Bob', 'Carol')
        grid on;
        set(gca,'FontSize',12);
    end
    
    if k > 0
        [vertex_k_ROIT, nr_ROIT] = con2vert(G_k_ROIT_total{k+kIndexC}, theta_k_ROIT_total{k+kIndexC});
        [CH_k_ROIT, x_ROIT_volume(k+kIndexC)] = convhull(vertex_k_ROIT);
        figure,
        plot(vertex_k_posterior(CH_k_posterior, 1), vertex_k_posterior(CH_k_posterior, 2), '-s', vertex_k_posterior_A(CH_k_posterior_A, 1), vertex_k_posterior_A(CH_k_posterior_A, 2), '-o',...
            vertex_k_posterior_B(CH_k_posterior_B, 1), vertex_k_posterior_B(CH_k_posterior_B, 2), '-*', vertex_k_posterior_C(CH_k_posterior_C, 1), vertex_k_posterior_C(CH_k_posterior_C, 2), '-^',...
            vertex_k_ROIT(CH_k_ROIT, 1), vertex_k_ROIT(CH_k_ROIT, 2), '--', 'LineWidth', 1.2, 'MarkerSize', 8)
        legend('True posterior', 'Alice', 'Bob', 'Carol', 'ROIT')
        grid on;
        set(gca,'FontSize',12);
    end
end