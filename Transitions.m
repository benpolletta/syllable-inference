function [t_count, prob] = Transitions(unit, tsylb_option)
if nargin < 1, unit = ''; end
if isempty(unit), unit = 'sylb'; end
if nargin < 2, tsylb_option = []; end
if isempty(tsylb_option), tsylb_option = 0; end

SI = (1:6300)/6300;

time = 0:.1:10000;

onset_time = 1000;

name = [unit, 'Transitions'];
if length(SI) ~= 6300
    name = sprintf('%s_%dsentences', name, length(SI));
end

wsp_map = load('word2sylb2phone_bysentence.mat');
wsp_map = wsp_map.results;

data = load([unit, 'Data.mat']);
unit_list = data.id;

num_units = length(unit_list);

t_count = zeros(num_units);

% phone2feature_data = load('phones2features.mat');
% phone_list = phone2feature_data.phone_list;
% phones2features = phone2feature_data.feature_mat;
   
timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

file_list_id = fopen([timit_dir, 'DOC/allphonelist_filenames.txt'], 'r');
file_list = textscan(file_list_id, '%s');
fclose(file_list_id);

file_list = file_list{1};

field = [unit, '_cell'];

for s = 1:length(SI)

    these_units = wsp_map(s).(field);
    
    %% Computing transition counts.
    
    for u = 1:(length(these_units) - 1)
        
        this_index = find(strcmpi(unit_list, these_units{u}));
        next_index = find(strcmpi(unit_list, these_units{u + 1}));
        t_count(next_index, this_index) = t_count(next_index, this_index) + 1;
        
    end

    % results(syl) = struct('transitions', {transitions}); %, 'transition_features', transition_features);
    
end

% t_count = sparse(t_count);
num_transitions=sum(sum(t_count));
fprintf('Number of transitions: %d.\n', num_transitions)

prob = t_count/num_transitions;
prob_cols = nanunitsum(t_count);
prob_rows = nanunitsum(t_count')';

save([name, '.mat'], 'unit_list', 't_count', 'prob', 'prob_rows', 'prob_cols')


plotTransitions([name, '_Counts'], [unit, ' Transition Counts'], unit_list, t_count)

saveas(gcf, [name, '_counts.fig'])

plotTransitions([name, '_Prob'], [unit, ' Transition Probabilities (col.)'], unit_list, prob_cols, 'log')

saveas(gcf, [name, '_prob_cols.fig'])

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

if strcmp(scaling, 'log')
    sizes = log(val); colors = log(val);
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

title(this_title)

saveas(gcf, [name, '.fig'])

save_as_pdf(gcf, name)

end