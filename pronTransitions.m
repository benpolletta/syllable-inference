function [pron_t_count, t_prob] = pronTransitions

SI = (1:6300)/6300;

time = 0:.1:10000;

onset_time = 1000;

name = 'pronTransitions';
if length(SI) ~= 6300
    name = sprintf('%s_%dsentences', name, length(SI));
end

wsp_map = load('word2sylb2phone_bysentence.mat');
wsp_map = wsp_map.results;

all_pron_list = cat(2, wsp_map.word_cell)';
all_word_list = cat(1, wsp_map.canonical_word_cell);
all_pair_list = cellfun(@(x, y) strjoin({x,y}, '='), all_word_list, all_pron_list, 'unif', 0);

list = unique(all_pron_list);
list(cellfun(@isempty, list)) = [];
list{end + 1} = '#';
pair_list = unique(all_pair_list);
pair_list{end + 1} = '#=#';

num_prons = length(list);
num_pairs = length(pair_list);

[p1_count, p2_count] = deal(zeros(num_prons, 1));
[wp1_count, wp2_count] = deal(zeros(num_pairs, 1));
pron_t_count = zeros(num_prons);
pair_t_count = zeros(num_pairs);

% phone2feature_data = load('phones2features.mat');
% phone_list = phone2feature_data.phone_list;
% phones2features = phone2feature_data.feature_mat;
   
timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

file_list_id = fopen([timit_dir, 'DOC/allphonelist_filenames.txt'], 'r');
file_list = textscan(file_list_id, '%s');
fclose(file_list_id);

file_list = file_list{1};

for s = 1:length(SI)

    prons = wsp_map(s).word_cell';
    words = wsp_map(s).canonical_word_cell;
    pairs = cellfun(@(x, y) strjoin({x,y}, '='), words, prons, 'unif', 0);
    pairs(2:(end + 1)) = pairs;
    pairs([1 end + 1]) = {'#=#', '#=#'};
    prons = cellfun(@(x) extractAfter(x, '='), pairs, 'UniformOutput', false);
    
    %% Computing transition counts.
    
    for p = 1:(length(pairs) - 1)
        
        this_index = find(strcmpi(list, prons{p}));
        next_index = find(strcmpi(list, prons{p + 1}));
        p1_count(this_index) = p1_count(this_index) + 1;
        p2_count(next_index) = p2_count(next_index) + 1;
        pron_t_count(next_index, this_index) = pron_t_count(next_index, this_index) + 1;
        
        this_index = find(strcmpi(pair_list, pairs{p}));
        next_index = find(strcmpi(pair_list, pairs{p + 1}));
        wp1_count(this_index) = wp1_count(this_index) + 1;
        wp2_count(next_index) = wp2_count(next_index) + 1;
        pair_t_count(next_index, this_index) = pair_t_count(next_index, this_index) + 1;
        
    end

    % results(syl) = struct('transitions', {transitions}); %, 'transition_features', transition_features);
    
end

% t_count = sparse(t_count);
num_transitions=sum(sum(pair_t_count));
fprintf('Number of transitions: %d.\n', num_transitions)

prob = pron_t_count/num_transitions;
first_prob = p1_count/sum(p1_count);
second_prob = p2_count/sum(p2_count);
t_prob = diag(1./first_prob)*prob;

pair_prob = pair_t_count/num_transitions;
wp1_prob = wp1_count/sum(wp1_count);
wp2_prob = wp2_count/sum(wp2_count);
pair_t_prob = diag(1./wp1_prob)*pair_prob;

save([name, '.mat'], 'list', 'pron_t_count', 'prob', 'first_prob', 'second_prob', 't_prob', 'pair_list', 'pair_t_count', 'pair_prob', 'wp1_prob', 'wp2_prob', 'pair_t_prob')

plotTransitions([name, '_pronCounts'], 'Pronunciation Transition Counts', list, pron_t_count)

plotTransitions([name, '_pronProb'], 'Pronunciation Transition Probabilities (col.)', list, t_prob, 'exp')

plotTransitions([name, '_pairCounts'], '(Word, Pron.) Pair Transition Counts', pair_list, pair_t_count)

plotTransitions([name, '_pairProb'], '(Word, Pron.) Pair Transition Probabilities (col.)', pair_list, pair_t_prob, 'exp')

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