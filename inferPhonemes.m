function inferPhonemes

sentence = getSentence;

methods = {'', 'meanpdf', 'prodpdf'};
num_methods = length(methods);

figure()

subplot(num_methods + 1, 1, 1)

plotSentence(gca, sentence)

window_length = 2501;

for m = 1:num_methods
    
    likelihood = calc_phone_likelihood(sentence, window_length, methods{m});
    
    subplot(num_methods + 1, 1, m + 1)
    
    imagesc(likelihood.time, likelihood.phones, likelihood.likelihood)
    
end

end


function plotSentence(plot_axis, sentence)

fields = fieldnames(sentence);

for f = 1:length(fields)
    
    eval(sprintf('%s = sentence.%s;', fields{f}, fields{f}))
    
end

feature_dim = size(input_vec, 2);

imagesc(plot_axis, time, 1:feature_dim, input_vec')%*diag(1./max(feature_mat)))')

hold on

plot(plot_axis, repmat(phone_transition_times, 2, 1), repmat([0; size(input_vec, 1)], 1, length(phone_transition_times)), 'w', 'LineWidth', .5)

xticks(phone_transition_times(1:(end - 1))+diff(phone_transition_times)/2)
xticklabels(phone_sequence)
xtickangle(45)

xlim([min(phone_transition_times), max(phone_transition_times)])

yticks(1:feature_dim)
yticklabels(feature_names)

axis xy

colorbar

end