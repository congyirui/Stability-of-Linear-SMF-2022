%%% Drawing Figures

close all;
clear;
clc;

load('experiment_os_20211024')


%%  Figure 5
k_sequence_drawing = k_sequence;

diameters_average_drawing = diameters_average(k_sequence_drawing + kIndexC);
diameters_max_drawing = diameters_max(k_sequence_drawing + kIndexC);
diameters_min_drawing = diameters_min(k_sequence_drawing + kIndexC);

bounds_estimation_gap_average_drawing = bounds_estimation_gap_average(k_sequence_drawing + kIndexC);
bounds_estimation_gap_max_drawing = bounds_estimation_gap_max(k_sequence_drawing + kIndexC);
bounds_estimation_gap_min_drawing = bounds_estimation_gap_min(k_sequence_drawing + kIndexC);


figure,
subplot(2, 1, 1)
shadedata=zeros(2,2*length(k_sequence_drawing));
shadedata(1,1:length(k_sequence_drawing))=k_sequence_drawing;
shadedata(1,length(k_sequence_drawing)+1:2*length(k_sequence_drawing))=k_sequence_drawing(end:-1:1);
shadecolor=[100,149,237]/255;
shadedata(2,1:length(k_sequence_drawing))=diameters_max_drawing;
reversedata=diameters_min_drawing;
shadedata(2,length(k_sequence_drawing)+1:2*length(k_sequence_drawing))=reversedata(end:-1:1);
h=fill(shadedata(1,:)',shadedata(2,:)',shadecolor);
set(h,'LineStyle','none')
hold on;
plot(k_sequence_drawing, diameters_average_drawing, 'LineWidth', 1.2)
set(gca,'FontSize',12);
legend('Diameter range', 'Averaged diameter')
xlabel('Time Step')
ylabel('Diameter')
grid on;

subplot(2, 1, 2)
shadedata=zeros(2,2*length(k_sequence_drawing));
shadedata(1,1:length(k_sequence_drawing))=k_sequence_drawing;
shadedata(1,length(k_sequence_drawing)+1:2*length(k_sequence_drawing))=k_sequence_drawing(end:-1:1);
shadecolor=[100,149,237]/255;
shadedata(2,1:length(k_sequence_drawing))=bounds_estimation_gap_max_drawing;
reversedata=bounds_estimation_gap_min_drawing;
shadedata(2,length(k_sequence_drawing)+1:2*length(k_sequence_drawing))=reversedata(end:-1:1);
h=fill(shadedata(1,:)',shadedata(2,:)',shadecolor);
set(h,'LineStyle','none')
hold on;
plot(k_sequence_drawing, bounds_estimation_gap_average_drawing, 'LineWidth', 1.2)
set(gca,'FontSize',12);
legend('Bound range', 'Averaged bound')
xlabel('Time Step')
ylabel('Estimation Gap Bound')
grid on;