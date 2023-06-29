function [t_count, prob] = wordTransitions

SI = (1:6300)/6300;

time = 0:.1:10000;

onset_time = 1000;

name = 'wordTransitions';
if length(SI) ~= 6300
    name = sprintf('%s_%dsentences', name, length(SI));
end

wsp_map = load('word2sylb2phone_bysentence.mat');
wsp_map = wsp_map.results;

word_data = load('wordCanonData.mat');
word_list = word_data.id;

num_words = length(word_list);

[w1_count, w2_count] = deal(zeros(num_words, 1));
t_count = zeros(num_words);

% phone2feature_data = load('phones2features.mat');
% phone_list = phone2feature_data.phone_list;
% phones2features = phone2feature_data.feature_mat;
   
timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

file_list_id = fopen([timit_dir, 'DOC/allphonelist_filenames.txt'], 'r');
file_list = textscan(file_list_id, '%s');
fclose(file_list_id);

file_list = file_list{1};

for s = 1:length(SI)

    words = wsp_map(s).canonical_word_cell;
    
    %% Computing transition counts.
    
    for w = 1:(length(words) - 1)
        
        this_index = find(strcmpi(word_list, words{w}));
        next_index = find(strcmpi(word_list, words{w + 1}));
        w1_count(this_index) = w1_count(this_index) + 1;
        w2_count(next_index) = w2_count(next_index) + 1;
        t_count(next_index, this_index) = t_count(next_index, this_index) + 1;
        
    end

    % results(syl) = struct('transitions', {transitions}); %, 'transition_features', transition_features);
    
end

% t_count = sparse(t_count);
num_transitions=sum(sum(t_count));
fprintf('Number of transitions: %d.\n', num_transitions)

prob = t_count/num_transitions;
w1_prob = w1_count/sum(w1_count);
w2_prob = w2_count/sum(w2_count);
t_prob = diag(1./w1_prob)*prob;

save([name, '.mat'], 'word_list', 't_count', 'prob', 'w1_prob', 'w2_prob', 't_prob')


plotTransitions([name, '_Counts'], 'Word Transition Counts', word_list, t_count)

plotTransitions([name, '_Prob'], 'Word Transition Probabilities (col.)', word_list, nanunitsum(t_count), 'exp')

end

function plotTransitions(name, this_title, units, transitions, scaling)
if nargin < 5, scaling = ''; end

% if sort_option
%     [prob, sort_order] = sort(prob, 'descend');
%     ids = ids(sort_order);
% 
%     these_x_ticks = int32(logspace(0, log10(length(ids)), 20));
% else
%     x_tick_step = round(length(ids)/20);
%     these_x_ticks = 1:x_tick_step:length(ids);
% end
% x_tick_labels = ids(these_x_ticks);

%% Plotting histograms.

figure

x = 1:length(units);

[row, col, val] = find(transitions);

if strcmp(scaling, 'sqrt')
    sizes = sqrt(val); colors = val;
elseif strcmp(scaling, 'log')
    sizes = log(val); colors = val;
elseif strcmp(scaling, 'exp')
    sizes = 20.^val; colors = val;
else
    sizes = val; colors = val;
end
scatter(row, col, sizes, colors, 'filled', 'LineWidth', .25)

%colormap('hot')

tick_step = round(length(units)/10);
these_ticks = 1:tick_step:length(units);

set(gca, 'YTick', these_ticks, 'YTickLabel', units(these_ticks))

set(gca, 'XTick', these_ticks, 'XTickLabel', units(these_ticks))

xtickangle(-45)

nochange_colorbar(gca)

axis tight

title(this_title)

saveas(gcf, [name, '.fig'])

save_as_pdf(gcf, name)

end