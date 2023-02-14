function results = syllableData(SI)

if nargin < 1, SI = []; end
if isempty(SI), SI = (1:6300)/6300; end

if length(SI) ~= 6300
    name = sprintf('Data_%dsentences', length(SI));
else
    name = 'Data';
end

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
    syllable_lengths = double(diff(tsylb_indices/16));
    norm_syllable_lengths = syllable_lengths/mean(syllable_lengths);
    
    %% Retrieving phonemes and their start and end times.
    
    phone_filename = [sentence_dir, file_name, '.PHN'];
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
        
        this_syl_indicator = mid_phone_indices > tsylb_indices(t) & mid_phone_indices < tsylb_indices(t + 1);
        
        syllables{t} = phones(this_syl_indicator);
        
    end
    
    %% Saving results.
    
    results(s) = struct('syllable_lengths', syllable_lengths, 'norm_syllable_lengths', norm_syllable_lengths,...
        'phone_lengths', phone_lengths, 'syllables', {syllables}, 'phones', {phones});
    
end

%% Collecting syllable lengths across sentences.

all_syllables = cat(1, results.syllables);

phoneno_vec = cellfun(@length, all_syllables);

firstphone_vec = cellfun(@(x) x{1}, all_syllables, 'unif', 0);

string_vec = cellfun(@(x) strjoin(x, '/'), all_syllables, 'unif', 0);

duration_vec = cat(1, results.syllable_lengths);

norm_duration_vec = cat(1, results.norm_syllable_lengths);

%% Plotting durations by syllable and syllable length (i.e., number of phones).

no_bins = ceil(sqrt(length(SI)));

vecs = {duration_vec, norm_duration_vec};

syl_grouping_vars = {string_vec, phoneno_vec, firstphone_vec};

syl_grouping_labels = {'', 'PhoneNum', 'FirstPhone'};

vec_labels = {'syl', 'normSyl'};

for v = 1:length(vecs)
    
    for g = 1:length(syl_grouping_vars)
        
        [index, id] = findgroups(syl_grouping_vars{g});
        
        fname = [vec_labels{v}, syl_grouping_labels{g}, name];
        
        %%% Getting statistics & computing histograms.
        
        [count, prob, stats, cdf, hist, bins, bin_centers] = calcStats(vecs{v}, index, id, no_bins, fname);
        
        figure()
        plot(1:length(id), prob/sum(prob), 'LineWidth', 2, 'Color', 'k')
        set(gca, 'XTick', 1:length(id), 'XTickLabel', id)
        xtickangle(45)
        axis tight
        title('Syllable Distribution')
        ylabel('Probability')
        
        saveas(gcf, [fname, '.fig'])
        
        plotHistograms(fname, id, bin_centers', hist')
        
    end
    
end

end

function [count, prob, stats, cdf, hist, bins, bin_centers] = calcStats(vec, index, id, no_bins, fname)

count = splitapply(@(x) length(x), vec, index);
prob = count/sum(count);
stats = splitapply(@(x) [mean(x), std(x), quantile(x, [.5, .25, .75])], vec, index);
cdf = splitapply(@(x){[sort(x), linspace(0, 1, length(x))']}, vec, index);

[hist, bins] = splitapply(@(x) histcounts(x, no_bins, 'Normalization', 'probability'), vec, index);

bin_centers = bins(:, 1:(end - 1)) + diff(bins, [], 2)/2;
bin_centers = [bin_centers(:, 1) - mean(diff(bins, [], 2), 2),...
    bin_centers, bin_centers(:, end) + mean(diff(bins, [], 2), 2)];

hist = [zeros(size(hist(:, 1))), hist, zeros(size(hist(:, end)))];

save([fname, '.mat'], 'id', 'count', 'prob', 'stats', 'cdf', 'hist', 'bins', 'bin_centers')

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