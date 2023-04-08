function results = syllablePhones(SI)

name = sprintf('LengthHistograms_%dsentences', length(SI));

for s = 1:length(SI)
    
    sentence_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';
    
    file_list_id = fopen([sentence_dir, 'wavFileList.txt'], 'r');
    file_list = textscan(file_list_id, '%s');
    fclose(file_list_id);
    
    file_list = file_list{1};
    
    file_index = SI(s);
    
    if isfloat(file_index) && file_index < 1 && file_index > 0
        file_index = round(file_index*length(file_list));
    end
    
    wavfile_name = file_list{file_index};
    file_name = extractBefore(wavfile_name, '.WAV');
    
    %% Retrieving syllable boundary_times.
    
    tsylb_filename = [sentence_dir, file_name, '.TSYLB'];
    fid = fopen(tsylb_filename, 'r');
    tsylb_indices = textscan(fid, '%d');
    fclose(fid);
    tsylb_indices = tsylb_indices{1};
    syllable_lengths = diff(tsylb_indices/16);
    
    %% Retrieving phonemes and their start and end times.
    
    phone_filename = [sentence_dir, file_name, '.tsylbPHN'];
    fid = fopen(phone_filename, 'r');
    phone_data = textscan(fid, '%s');
    fclose(fid);
    phone_data = phone_data{1};
    phone_data = reshape(phone_data, 3, length(phone_data)/3);
    
    phones = phone_data(3, :);
    
    phone_indices = cellfun(@str2num, phone_data(1:2, :));
    mid_phone_indices = mean(phone_indices);
    
    phone_times = (phone_indices/16)';
    phone_lengths = diff(phone_times, [], 2);
    
    syllables = cell(size(syllable_lengths));
    
    for t = 1:length(syllable_lengths)
        
        syl_phone_indicator = mid_phone_indices > tsylb_indices(t) & mid_phone_indices < tsylb_indices(t + 1);
        
        syllables{t} = phones(syl_phone_indicator);
        
    end
    
    %% Saving results.
    
    results(s) = struct('syllable_lengths', syllable_lengths, 'phone_lengths', phone_lengths, 'syllables', {syllables}, 'phones', {phones});
    
end

all_syllables = cat(1, results.syllables);

all_syllables = cellfun(@(x) strjoin(x, '/'), all_syllables, 'unif', 0);

all_syl_lengths = cat(1, results.syllable_lengths);

%% Collecting phone lengths across sentences.

%%% Grouping by phonemes.

[group_index, group_id] = findgroups(phone_vec);

%%% Converting from group_id order to timit_phonemes order.

timit2group = cellfun(@(x) strcmp(x, group_id), timit_phonemes, 'unif', 0);
timit2group_mat = cat(2, timit2group{:});

present_phonemes = timit_phonemes;
present_phonemes(sum(timit2group_mat) == 0) = [];

present2group = cellfun(@(x) strcmp(x, group_id), present_phonemes, 'unif', 0);
present2group_mat = cat(2, present2group{:});

[present_map, ~] = find(present2group_mat');
phone_index = present_map(group_index);

%%% Finding means across phonemes.

phone_spect_mean = splitapply(@mean, phone_spect_vec, phone_index);

% Grouping by phoneme class.

[timit_map, ~] = find(timit2group_mat');
[class_map, ~] = find(class_indicator');

class_index = class_map(timit_map(group_index));

[classLength_counts, classLength_bins] = splitapply(@(x) histcounts(x, no_bins, 'Normalization', 'probability'), phone_length_vec, class_index);

classLength_bin_centers = classLength_bins(:, 1:(end - 1)) + diff(classLength_bins, [], 2)/2;
classLength_bin_centers = [classLength_bin_centers(:, 1) - mean(diff(classLength_bins, [], 2), 2),...
    classLength_bin_centers, classLength_bin_centers(:, end) + mean(diff(classLength_bins, [], 2), 2)];

classLength_counts = [zeros(size(classLength_counts(:, 1))), classLength_counts, zeros(size(classLength_counts(:, end)))];

plotHistograms(['class', name], class_names(1:(end - 1)), classLength_bin_centers(1:(end - 1), :)', classLength_counts(1:(end - 1), :)')

% group_indicator_mat = zeros(length(group_index), length(group_id));
% gim_linear = sub2ind(size(group_indicator_mat), 1:length(group_index), group_index');
% group_indicator_mat(gim_linear) = 1;
% 
% % group_class_perm = timit2group_mat*double(class_indicator);
% 
% class_indicator_mat = group_indicator_mat*class_indicator; %group_class_perm;
% 
% [class_indicator_rows, class_indicator_cols] = find(class_indicator_mat);
% 
% class_index = nan(size(class_indicator_rows));
% class_index(class_indicator_rows) = class_indicator_cols;

save([name, '.mat'], 'results', 'phoneLength_counts', 'phoneLength_bin_centers', 'classLength_counts', 'classLength_bin_centers')

end

function plotHistograms(name, class_names, hist_bin_centers, hist_counts)

%% Plotting histograms.

figure

x = cumsum(ones(length(class_names), size(hist_bin_centers, 1)))';

pcolor(hist_bin_centers, x, hist_counts)

colormap('hot')

set(gca, 'YTick', 1:length(class_names), 'YTickLabel', class_names)

shading interp

% xtickangle(-45)

% set(gca, 'YTick', 1:size(hists, 2), 'YTickLabel', hist_labels)

% ytickangle(-90)

nochange_colorbar(gca)

% rows = size(hists, 2);
% 
% ha = tight_subplot(rows, 1);
% 
% for i = 1:rows
%     
%     axes(ha(i))
%     
%     bar(hists(:, i))
%     
%     box off
%     
%     set(gca, 'XTickLabel', '')
%     
%     if i == rows
%         
%         set(gca, 'XTickLabel', class_names)
%         
%         xtickangle(-45)
%         
%     end
%     
%     ylabel(hist_labels{i}, 'Rotation', 0)
%     
% end

saveas(gcf, [name, '_BPhistograms.fig'])

save_as_pdf(gcf, [name, '_BPhistograms'])

figure

colors = flipud(hsv(length(class_names)));

set(gca, 'NextPlot', 'add', 'ColorOrder', colors)

plot(hist_bin_centers, hist_counts, 'LineWidth', 2)

legend(class_names)

saveas(gcf, [name, '_BPhistograms_line.fig'])

save_as_pdf(gcf, [name, '_BPhistograms_line'])

end