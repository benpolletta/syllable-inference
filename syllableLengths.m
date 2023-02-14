function results = syllableLengths(SI)

name = sprintf('sylLengthHistograms_%dsentences', length(SI));

for s = 1:length(SI)
    
    sentence_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';
    
    file_list_id = fopen([sentence_dir, 'wavFileList.txt'], 'r');
    file_list = textscan(file_list_id, '%s');
    fclose(file_list_id);
    
    file_list = file_list{1};
    
    file_index = SI(s);
    
    if isfloat(file_index) && file_index <= 1 && file_index > 0
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

%% Collecting syllable lengths across sentences.

all_syllables = cat(1, results.syllables);

syl_phoneno_vec = cellfun(@length, all_syllables);

syl_firstphone_vec = cellfun(@(x) x{1}, all_syllables, 'unif', 0);

syl_string_vec = cellfun(@(x) strjoin(x, '/'), all_syllables, 'unif', 0);

syl_duration_vec = cat(1, results.syllable_lengths);

%% Plotting durations by syllable and syllable length (i.e., number of phones).

syl_grouping_vars = {syl_string_vec, syl_phoneno_vec, syl_firstphone_vec};

syl_grouping_labels = {'syllable', 'numPhones', 'firstPhone'};

sylLengthCell = cell(2,length(syl_grouping_vars));

for g = 1:length(syl_grouping_vars)
    
    [syl_index, syl_id] = findgroups(syl_grouping_vars{g});
    
    %%% Getting statistics.
    
    count{g} = splitapply(@(x) length(x), syl_duration_vec, syl_index);
    prob{g} = count/sum(count);
    stats{g} = splitapply(@(x) [mean(x), std(x), quantile(x, [.5, .25, .75])], syl_duration_vec, syl_index);
    cdf{g} = splitapply(@(x){[sort(x), linspace(0, 1, length(x))']}, syl_duration_vec, syl_index);
    
    %%% Finding histograms.
    
    no_bins = ceil(sqrt(length(SI)));
    
    [hists, bins] = splitapply(@(x) histcounts(x, no_bins, 'Normalization', 'probability'), syl_duration_vec, syl_index);
    
    bin_centers = bins(:, 1:(end - 1)) + diff(bins, [], 2)/2;
    bin_centers = [bin_centers(:, 1) - mean(diff(bins, [], 2), 2),...
        bin_centers, bin_centers(:, end) + mean(diff(bins, [], 2), 2)];
    
    hists = [zeros(size(hists(:, 1))), hists, zeros(size(hists(:, end)))];
    
    plotHistograms([syl_grouping_labels{g}, name], syl_id, bin_centers', hists')
    
    sylLengthCell{1, g} = bin_centers;
    sylLengthCell{2, g} = hists;
    
end

save([name, '.mat'], 'results', 'count', 'prob', 'stats', 'cdf', 'sylLengthCell')

end

function plotHistograms(name, class_names, hist_bin_centers, hist_counts)

%% Plotting histograms.

figure

x = cumsum(ones(length(class_names), size(hist_bin_centers, 1)))';

pcolor(hist_bin_centers, x, hist_counts)

colormap('hot')

set(gca, 'YTick', 1:length(class_names), 'YTickLabel', class_names)

shading interp

title('Syllable Duration Distribution')

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

saveas(gcf, [name, '.fig'])

save_as_pdf(gcf, name)

figure

colors = flipud(hsv(length(class_names)));

set(gca, 'NextPlot', 'add', 'ColorOrder', colors)

plot(hist_bin_centers, hist_counts, 'LineWidth', 2)

if isnumeric(class_names)
    class_legend = mat2cell(class_names, ones(size(class_names, 1), 1), ones(size(class_names, 2), 1));
    class_legend = cellfun(@num2str, class_legend, 'unif', 0);
    class_names = class_legend;
end

legend(class_names)

axis tight

title('Syllable Duration Distribution')

saveas(gcf, [name, '_line.fig'])

save_as_pdf(gcf, [name, '_line'])

end