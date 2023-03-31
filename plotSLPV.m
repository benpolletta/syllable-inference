function [sentence, vowel_lh, likelihood, posterior] = plotSLPV(method, vowel_window)

if nargin < 1, method = ''; end
if nargin < 2, vowel_window = []; end
if isempty(vowel_window), vowel_window = 10000; end

%% Getting sentence.

sentence = getSentence;

figure()

subplot(3, 1, 1)

plotSentence(gca, sentence)

vowel_lh = calc_vowel_likelihood(sentence, vowel_window);
scaled_vl = (vowel_lh.vowel_likelihood - min(vowel_lh.vowel_likelihood))*(length(vowel_lh.feature_names)/range(vowel_lh.vowel_likelihood));
plot(gca, vowel_lh.time, scaled_vl, 'w', 'LineWidth', 0.5)

D2_vowel_lh = diff(diff(vowel_lh.vowel_likelihood));
scaled_dvl = (D2_vowel_lh - min(D2_vowel_lh))*(length(vowel_lh.feature_names)/range(D2_vowel_lh));
plot(gca, vowel_lh.indicator_time, scaled_dvl, 'w--', 'LineWidth', 0.5)

vowel_ind = vowel_lh.v_indicator;
scaled_vi = (vowel_ind - min(vowel_ind))*(length(vowel_lh.feature_names)/range(vowel_ind));
plot(gca, vowel_lh.indicator_time, scaled_vi, 'y', 'LineWidth', 0.5)

window_length = 500; % 2501;

[for_plot(1:2)] = deal(sentence);

%% Calculating likelihood.
likelihood = calc_phone_likelihood_new(sentence, window_length);

for_plot(1).time = likelihood.time;
for_plot(1).input_vec = likelihood.phone_likelihood;
for_plot(1).feature_names = likelihood.phones;

subplot(3, 1, 2)

plotSentence(gca, for_plot(1))

hold on

scaled_tl = (likelihood.trans_likelihood - min(likelihood.trans_likelihood))*(length(likelihood.phones)/range(likelihood.trans_likelihood));
plot(gca, likelihood.time, scaled_tl, 'w', 'LineWidth', 0.5)

scaled_t2 = (likelihood.delta_likelihood - min(likelihood.delta_likelihood))*(length(likelihood.phones)/range(likelihood.delta_likelihood));
plot(gca, likelihood.time, scaled_t2, 'y', 'LineWidth', 0.5)

scaled_ie = (likelihood.inv_entropy - min(likelihood.inv_entropy))*(length(likelihood.phones)/range(likelihood.inv_entropy));
plot(gca, likelihood.time, scaled_ie, 'w--', 'LineWidth', 0.5)

%% Calculating posterior.
posterior = calc_phone_posterior_new(likelihood, method);

for_plot(2).time = likelihood.time;
for_plot(2).input_vec = nanunitsum(posterior.phone_posterior);
for_plot(2).feature_names = likelihood.phones;

subplot(3, 1, 3)

plotSentence(gca, for_plot(2))

hold on

scaled_t3 = (posterior.trans_posterior - min(posterior.trans_posterior))*(length(likelihood.phones)/range(posterior.trans_posterior));
plot(gca, likelihood.time, scaled_t3, 'w', 'LineWidth', 0.5)

scaled_h = (posterior.hazard_est - min(posterior.hazard_est))*(length(likelihood.phones)/range(posterior.hazard_est));
plot(gca, likelihood.time, scaled_h, 'y', 'LineWidth', 0.5)

scaled_hv = (posterior.hazard_var - min(posterior.hazard_var))*(length(likelihood.phones)/range(posterior.hazard_var));
plot(gca, likelihood.time, scaled_hv, 'y--', 'LineWidth', 0.5)

scaled_ie = (posterior.inv_entropy - min(posterior.inv_entropy))*(length(likelihood.phones)/range(posterior.inv_entropy));
plot(gca, likelihood.time, scaled_ie, 'w--', 'LineWidth', 0.5)

tree = split(sentence.filename, '/');

saveas(gcf, sprintf('plotSLPV_%s%s.fig', tree{end}, method))

end


function plotSentence(plot_axis, sentence)

fields = fieldnames(sentence);

for f = 1:length(fields)
    
    eval(sprintf('%s = sentence.%s;', fields{f}, fields{f}))
    
end

feature_dim = length(feature_names);

if size(input_vec, 2) ~= feature_dim
    input_vec = input_vec';
end

imagesc(plot_axis, time, 1:feature_dim, input_vec')%*diag(1./max(feature_mat)))')
axis xy

hold on

plot(repmat(phone_transition_times(:), 1, 2)', repmat([0; size(input_vec, 1)], 1, length(phone_transition_times(:))), 'w', 'LineWidth', .5)
    
xticks(phone_transition_times(1, :) + diff(phone_transition_times)/2)
xticklabels(phone_sequence)
xtickangle(45)

xlim([min(min(phone_transition_times)), max(max(phone_transition_times))])

step = floor(feature_dim/10);
yticks(1:step:feature_dim)
yticklabels(feature_names(1:step:feature_dim))

title(join(sentence.word_sequence, ' '))

colorbar

end