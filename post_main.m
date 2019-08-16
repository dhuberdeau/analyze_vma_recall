% Must run R-scripts first.

load('sigmoid_fit_E1.mat')

%% plot dir. err. aligned to minPT, and catch trial err., and pr(corr), and variability
figure;
n_subplots = 4;

subplot(n_subplots,1,1); hold on;
plot([0 0], [-200 200], 'k-')
plot(view_time_all_0, de_all_0, '.', 'Color', RED_COLOR) %type 0
plot(view_time_all_1, de_all_1, '.', 'Color', GREEN_COLOR) %type 1
plot(view_time_all_2, de_all_2, '.', 'Color', BLUE_COLOR) %type 2
plot([-.5 .5], [30 30], '-', 'Color', [.5 .5 .5]);
plot([-.5 .5], -[30 30], '-', 'Color', [.5 .5 .5]);
axis([-0.5 0.5 -200 200]);
set(gca,'fontsize',18)
ylabel('Directional error', 'Fontsize', 12);

subplot(n_subplots,1,2); hold on;
plot([0 0], [-200 200], 'k-')
plot(view_time_all_3, de_all_3, 'o', 'Color', GREEN_COLOR) %type 3 (catch trials, direct cue)
plot(view_time_all_4, de_all_4, 'o', 'Color', BLUE_COLOR) %type 4 (catch trials, symbol cue) 
plot([-.5 .5], [30 30], '-', 'Color', [.5 .5 .5]);
plot([-.5 .5], -[30 30], '-', 'Color', [.5 .5 .5]);
axis([-0.5 0.5 -200 200]);
set(gca,'fontsize',18)
ylabel('Directional error','Fontsize', 12);

subplot(n_subplots,1,3); hold on;
plot([0 0], [0 1], 'k-')
plot(x_ind, p_all_0, '.-', 'Color', RED_COLOR);
plot(x_ind, p_all_1, '.-', 'Color', GREEN_COLOR);
plot(x_ind, p_all_2, '.-', 'Color', BLUE_COLOR);
plot(x_ind, p_all_3, 'o-', 'Color', GREEN_COLOR);
plot(x_ind, p_all_4, 'o-', 'Color', BLUE_COLOR);
axis([-0.5 0.5 -0.05 1.05]);
set(gca,'fontsize',18)
plot([-.5 .5], [.25 .25], '-', 'LineWidth', 2, 'Color', [.5 .5 .5]);
plot([-.5 .5], [30 30], '--', 'LineWidth', 2, 'Color', [.5 .5 .5]);
plot([-.5 .5], -[30 30], '--', 'LineWidth', 2, 'Color', [.5 .5 .5]);
ylabel('Probability correct', 'Fontsize', 12)

% subplot(4,1,4); hold on;
% plot([0 0], [0 12], 'k-')
% plot(x_ind, var_all_0, 'r.-');
% plot(x_ind, var_all_1, 'g.-');
% plot(x_ind, var_all_2, 'b.-');

if n_subplots > 3
    % Plot the sigmoid fits if desired:
    subplot(n_subplots,1,4); hold on
    plot(mat_out.t, mat_out.y_0, 'Color', RED_COLOR, 'LineWidth', 2)
    plot(mat_out.t, mat_out.y_1, 'Color', GREEN_COLOR, 'LineWidth', 2)
    plot(mat_out.t, mat_out.y_2, 'Color', BLUE_COLOR, 'LineWidth', 2)
    plot(mat_out.t, mat_out.y_3, 'k--', 'LineWidth', 2)
    axis([-0.5 0.5 -0.05 1.05]);
    set(gca,'fontsize',18)
    ylabel('Probability correct', 'Fontsize', 12);
    desired_pos = [552 538 500 600];
else
    desired_pos = [552 538 500 458];
end
xlabel('Preparation time (sec)');

ff = gcf;
set(ff, 'Position', desired_pos);
set(ff, 'PaperOrientation', 'landscape')
saveas(ff, 'Aggregated_catchtrials.pdf');



%% Run Rscript for sigmoid analysis:
experiment_names = {'3T', '4T', '6T', '3Tb'};
quarter_index = 0:4;
plot_index = [1:5; 6:10; 11:15; 16:20]';
f_sig1 = figure;
for i_grp = 1:4
    for i_qrt = 1:5
        load(['sigmoid_fit_E2_', experiment_names{i_grp}, '_', num2str(quarter_index(i_qrt)), '.mat']);

        subplot(4,5,plot_index(i_qrt, i_grp)); hold on;

        % important note: 
        % y_1 -> type 1
        % y_2 -> type 2
        % y_3 -> catch trials (type 3 or 4)
        if i_qrt == 1
            plot(sig_out.t, sig_out.y_0, 'LineWidth', 2, 'Color', RED_COLOR);
        else
            plot(sig_out.t, sig_out.y_1, 'LineWidth', 2, 'Color', GREEN_COLOR);
            plot(sig_out.t, sig_out.y_2, 'LineWidth', 2, 'Color', BLUE_COLOR);
            plot(sig_out.t, sig_out.y_3, 'k', 'LineWidth', 2);
        end
    end
end
set(f_sig1, 'Position', [72 549 818 406]);
set(f_sig1, 'PaperOrientation', 'landscape')
saveas(f_sig1, 'Groups_sigmoid_fit.pdf');