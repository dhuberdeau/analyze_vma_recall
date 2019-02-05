% David Huberdeau 01.23.2019
%
% Fit a q-learning working memory reinforcement learning model (qwmrl) from
% data from participants.
%
%

% must have run analyze_statLearn_v2.m

N_valid = 285;
% Experiment 1:
actions = mov_class_all{1}(1:N_valid, :);
s_input = target_all{1}(1:N_valid, :);

beta_0 = [0.6; 20; 1; 45; 0.25; .95];

p_hat = nan(6, size(actions,2));
for i_sub = 1:size(actions, 2)
    p_hat(:,i_sub) = mle(s_input(:, i_sub), 'logpdf', @q_learn_rlwm_pdf, 'start', beta_0);
end

%%
% Experiment 2
actions = mov_class_all{2}(1:N_valid, :);
s_input = target_all{2}(1:N_valid, :);

beta_0 = [0.6; 20; 1; 45; 0.25; .95];

p_hat = nan(6, size(actions,2));
for i_sub = 1:size(actions, 2)
    p_hat(:,i_sub) = mle(s_input(:, i_sub), 'logpdf', @q_learn_rlwm_pdf, 'start', beta_0);
end