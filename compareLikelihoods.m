function likelihoods = compareLikelihoods

sentence = getSentence;

functions = {@calc_phone_likelihood_new}; % , @calc_phone_transition_likelihood};
methods = {''}; %,''}; %, 'meanpdf', 'prodpdf'};
num_methods = length(methods);

figure()

subplot(num_methods + 1, 1, 1)

plotSentence(gca, sentence)

window_length = 500; % 2501;

[likelihoods(1:num_methods)] = deal(sentence);

for m = 1:num_methods
    
    likelihood = feval(functions{m}, sentence, window_length, methods{m});
    
    likelihoods(m).time = likelihood.time;
    likelihoods(m).input_vec = likelihood.likelihood;
    if isfield(likelihood, 'transition_list')
        likelihoods(m).feature_names = likelihood.transition_list(:);
    else
        likelihoods(m).feature_names = likelihood.phones;
    end
    
    subplot(num_methods + 1, 1, m + 1)
    
    plotSentence(gca, likelihoods(m))

    if isfield(likelihood, 'trans_likelihood')

        hold on

        scaled_tl = (likelihood.trans_likelihood - min(likelihood.trans_likelihood))*(length(likelihood.phones)/range(likelihood.trans_likelihood));

        plot(gca, likelihood.time, scaled_tl, 'w', 'LineWidth', 0.5)

    end

    if isfield(likelihood, 'delta_likelihood')
        
        scaled_t2 = (likelihood.delta_likelihood - min(likelihood.delta_likelihood))*(length(likelihood.phones)/range(likelihood.delta_likelihood));

        plot(gca, likelihood.time, scaled_t2, 'y', 'LineWidth', 0.5)

    end

    if isfield(likelihood, 'inv_entropy')
        
        scaled_ie = (likelihood.inv_entropy - min(likelihood.inv_entropy))*(length(likelihood.phones)/range(likelihood.inv_entropy));

        plot(gca, likelihood.time, scaled_ie, 'w--', 'LineWidth', 0.5)

    end
    
%     likelihood = calc_phone_transition_likelihood(sentence, window_length, methods{m});
%     
%     likelihoods(2*m).time = likelihood.time;
%     likelihoods(2*m).input_vec = likelihood.likelihood;
%     likelihoods(2*m).feature_names = likelihood.phones;
%     
%     subplot(2*num_methods + 1, 1, 2*m + 1)
%     
%     plotSentence(gca, likelihoods(2*m))
    
end

tree = split(sentence.filename, '/');

saveas(gcf, sprintf('compareLikelihoods_%s', tree{end}))

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

plot(plot_axis, repmat(phone_transition_times, 2, 1), repmat([0; size(input_vec, 1)], 1, length(phone_transition_times)), 'w', 'LineWidth', .5)

xticks(phone_transition_times(1:(end - 1))+diff(phone_transition_times)/2)
xticklabels(phone_sequence)
xtickangle(45)

xlim([min(phone_transition_times), max(phone_transition_times)])

step = floor(feature_dim/10);
yticks(1:step:feature_dim)
yticklabels(feature_names(1:step:feature_dim))

title(join(sentence.word_sequence, ' '))

colorbar

end