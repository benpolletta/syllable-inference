function [count, results] = phoneTransitions(SI)

if nargin < 1, SI = []; end
if isempty(SI), SI = (1:6300)/6300; end

time = 0:.1:10000;

onset_time = 1000;

name = 'phoneTransitions';
if length(SI) ~= 6300
    name = sprintf('%s_%dsentences', name, length(SI));
end

[tsylb_phonemes, class_indicator, class_names] = getTIMITphones(1);

num_phones = length(tsylb_phonemes);

count = zeros(num_phones);

% phone2feature_data = load('phones2features.mat');
% phone_list = phone2feature_data.phone_list;
% phones2features = phone2feature_data.feature_mat;
   
timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

file_list_id = fopen([timit_dir, 'DOC/allphonelist_filenames.txt'], 'r');
file_list = textscan(file_list_id, '%s');
fclose(file_list_id);

file_list = file_list{1};

for s = 1:length(SI)
    
    file_index = SI(s);
    
    if isfloat(file_index) && file_index <= 1 && file_index > 0
        file_index = round(file_index*length(file_list));
    end
    
    file_name = file_list{file_index};
    
    %% Retrieving phonemes and their start and end times.
    
    phone_filename = [timit_dir, file_name, '.PHN'];
    fid = fopen(phone_filename, 'r');
    phone_data = textscan(fid, '%s');
    fclose(fid);
    phone_data = phone_data{1};
    phone_data = reshape(phone_data, 3, length(phone_data)/3);
    
    phones = phone_data(3, :);
    
%     phone_features = cell2mat(cellfun(@(x) phones2features(find(strcmpi(x, phone_list), 1), :), phones, 'unif', 0));
    
    transitions = cellfun(@(x, y) sprintf('%s2%s', x, y), phones(1:(end - 1)), phones(2:end), 'unif', 0);
    
%     transition_features = diff(phone_features);
    
    for p = 1:(length(phones) - 1)
        
        this_index = find(strcmp(tsylb_phonemes, phones{p}));
        next_index = find(strcmp(tsylb_phonemes, phones{p + 1}));
        count(next_index, this_index) = count(next_index, this_index) + 1;
        
    end

    results(s) = struct('transitions', {transitions}); %, 'transition_features', transition_features);
    
end

% transition_vec = cat(2, results.transitions);
% 
% [transition_index, transition_id] = findgroups(transition_vec);

num_transitions=sum(sum(count));
fprintf('Number of transitions: %d.\n', num_transitions)

class_count = (class_indicator')*count*class_indicator;
class_prob = class_count/num_transitions;

prob = count/num_transitions;

save([name, '.mat'], 'results', 'tsylb_phonemes', 'count', 'prob', 'class_names', 'class_count', 'class_prob')

plotTransitions([name, '_Counts'], 'Phoneme Transition Counts', tsylb_phonemes, count)

plotTransitions([name, '_Prob'], 'Phoneme Transition Probabilities (col.)', tsylb_phonemes, nanunitsum(count))

plotTransitions([name,'_CountsByClass'], 'Class Transition Counts', class_names, class_count)

plotTransitions([name,'_ProbByClass'], 'Class Transition Probabilities (col.)', class_names, nanunitsum(class_count))

end

function plotTransitions(name, this_title, timit_phonemes, phone_transitions)

%% Plotting histograms.

figure

x = 1:length(timit_phonemes);

imagesc(x, x, phone_transitions)

%colormap('hot')

set(gca, 'YTick', 1:length(timit_phonemes), 'YTickLabel', timit_phonemes)

set(gca, 'XTick', 1:length(timit_phonemes), 'XTickLabel', timit_phonemes)

xtickangle(-45)

nochange_colorbar(gca)

title(this_title)

saveas(gcf, [name, '.fig'])

save_as_pdf(gcf, name)

end