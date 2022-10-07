%%% Drawing Figures

close all;
clear;
clc;

load('experiment20211024')


%%  Averaged Running Times
averaged_running_time_description
averaged_running_time_feasible_point


%%  Figure 5
k_sequence_drawing = k_sequence;

error_2_norm_average_drawing = error_2_norm_average(k_sequence_drawing + kIndexC);
error_2_norm_max_drawing = error_2_norm_max(k_sequence_drawing + kIndexC);
error_2_norm_min_drawing = error_2_norm_min(k_sequence_drawing + kIndexC);

error_inf_norm_average_drawing = error_inf_norm_average(k_sequence_drawing + kIndexC);
error_inf_norm_max_drawing = error_inf_norm_max(k_sequence_drawing + kIndexC);
error_inf_norm_min_drawing = error_inf_norm_min(k_sequence_drawing + kIndexC);


figure,
shadedata=zeros(2,2*length(k_sequence_drawing));
shadedata(1,1:length(k_sequence_drawing))=k_sequence_drawing;
shadedata(1,length(k_sequence_drawing)+1:2*length(k_sequence_drawing))=k_sequence_drawing(end:-1:1);
shadecolor=[100,149,237]/255;
shadedata(2,1:length(k_sequence_drawing))=error_2_norm_max_drawing;
reversedata=error_2_norm_min_drawing;
shadedata(2,length(k_sequence_drawing)+1:2*length(k_sequence_drawing))=reversedata(end:-1:1);
h=fill(shadedata(1,:)',shadedata(2,:)',shadecolor);
set(h,'LineStyle','none')
hold on;
plot(k_sequence_drawing, error_2_norm_average_drawing, 'LineWidth', 1.2)
set(gca,'FontSize',12);
% legend('Distance range', 'Averaged distance')
xlabel('Time Step')
ylabel('Distance')
grid on;

k_sequence_drawing_highlight = k_sequence_drawing(50: end);

error_2_norm_average_drawing_highlight = error_2_norm_average_drawing(k_sequence_drawing_highlight);
error_2_norm_max_drawing_highlight = error_2_norm_max_drawing(k_sequence_drawing_highlight);
error_2_norm_min_drawing_highlight = error_2_norm_min_drawing(k_sequence_drawing_highlight);

axes('position', [0.374999999999998,0.384285714285717,0.499999999999999,0.500000000000003])
box on % put box around new pair of axes
shadedata=zeros(2,2*length(k_sequence_drawing_highlight));
shadedata(1,1:length(k_sequence_drawing_highlight))=k_sequence_drawing_highlight;
shadedata(1,length(k_sequence_drawing_highlight)+1:2*length(k_sequence_drawing_highlight))=k_sequence_drawing_highlight(end:-1:1);
shadecolor=[100,149,237]/255;
shadedata(2,1:length(k_sequence_drawing_highlight))=error_2_norm_max_drawing_highlight;
reversedata=error_2_norm_min_drawing_highlight;
shadedata(2,length(k_sequence_drawing_highlight)+1:2*length(k_sequence_drawing_highlight))=reversedata(end:-1:1);
h=fill(shadedata(1,:)',shadedata(2,:)',shadecolor);
set(h,'LineStyle','none')
hold on;
plot(k_sequence_drawing_highlight, error_2_norm_average_drawing_highlight, 'LineWidth', 1.2)
set(gca,'FontSize',12);
legend('Distance range', 'Averaged distance')
% xlabel('Time Step')
% ylabel('Distance (2-norm)')
axis tight
grid on;


figure,
shadedata=zeros(2,2*length(k_sequence_drawing));
shadedata(1,1:length(k_sequence_drawing))=k_sequence_drawing;
shadedata(1,length(k_sequence_drawing)+1:2*length(k_sequence_drawing))=k_sequence_drawing(end:-1:1);
shadecolor=[100,149,237]/255;
shadedata(2,1:length(k_sequence_drawing))=error_inf_norm_max_drawing;
reversedata=error_inf_norm_min_drawing;
shadedata(2,length(k_sequence_drawing)+1:2*length(k_sequence_drawing))=reversedata(end:-1:1);
h=fill(shadedata(1,:)',shadedata(2,:)',shadecolor);
set(h,'LineStyle','none')
hold on;
plot(k_sequence_drawing, error_inf_norm_average_drawing, 'LineWidth', 1.2)
set(gca,'FontSize',12);
% legend('Distance range', 'Averaged distance')
xlabel('Time Step')
ylabel('Distance (\infty-norm)')
grid on;

k_sequence_drawing_highlight = k_sequence_drawing(50: end);

error_inf_norm_average_drawing_highlight = error_inf_norm_average_drawing(k_sequence_drawing_highlight);
error_inf_norm_max_drawing_highlight = error_inf_norm_max_drawing(k_sequence_drawing_highlight);
error_inf_norm_min_drawing_highlight = error_inf_norm_min_drawing(k_sequence_drawing_highlight);

axes('position', [0.374999999999998,0.384285714285717,0.499999999999999,0.500000000000003])
box on % put box around new pair of axes
shadedata=zeros(2,2*length(k_sequence_drawing_highlight));
shadedata(1,1:length(k_sequence_drawing_highlight))=k_sequence_drawing_highlight;
shadedata(1,length(k_sequence_drawing_highlight)+1:2*length(k_sequence_drawing_highlight))=k_sequence_drawing_highlight(end:-1:1);
shadecolor=[100,149,237]/255;
shadedata(2,1:length(k_sequence_drawing_highlight))=error_inf_norm_max_drawing_highlight;
reversedata=error_inf_norm_min_drawing_highlight;
shadedata(2,length(k_sequence_drawing_highlight)+1:2*length(k_sequence_drawing_highlight))=reversedata(end:-1:1);
h=fill(shadedata(1,:)',shadedata(2,:)',shadecolor);
set(h,'LineStyle','none')
hold on;
plot(k_sequence_drawing_highlight, error_inf_norm_average_drawing_highlight, 'LineWidth', 1.2)
set(gca,'FontSize',12);
legend('Distance range', 'Averaged distance')
% xlabel('Time Step')
% ylabel('Distance (2-norm)')
axis tight
grid on;