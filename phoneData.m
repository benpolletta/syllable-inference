function results = phoneData(SI)

if nargin < 1, SI = []; end
if isempty(SI), SI = (1:6300)/6300; end

time = 0:.1:10000;

onset_time = 1000;

name = 'Data';
if length(SI) ~= 6300
    name = sprintf('%s_%dsentences', name, length(SI));
end

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
    
    %% Retrieving phonemes and their start and end times.
    
    phone_filename = [sentence_dir, file_name, '.PHN'];
    fid = fopen(phone_filename, 'r');
    phone_data = textscan(fid, '%s');
    fclose(fid);
    phone_data = phone_data{1};
    phone_data = reshape(phone_data, 3, length(phone_data)/3);
    
    phones = phone_data(3, :);
    
    phone_indices = cellfun(@str2num, phone_data(1:2, :));
    phone_times = ((phone_indices/16 + onset_time)/1000)';
    
    phone_lengths = diff(phone_times, [], 2); 
    
    %% Retrieving syllable boundary_times & normalizing phone lengths.
    
    tsylb_filename = [sentence_dir, file_name, '.TSYLB'];
    fid = fopen(tsylb_filename, 'r');
    tsylb_indices = textscan(fid, '%d');
    fclose(fid);
    tsylb_indices = tsylb_indices{1};
    syllable_times = (tsylb_indices/16 + onset_time)/1000;
    syllable_lengths = diff(syllable_times);
    
    norm_phone_lengths = phone_lengths/mean(syllable_lengths);
    
    %% Saving results.
    
    results(s) = struct('phone_lengths', phone_lengths, 'phones', {phones},...
        'syllable_lengths', syllable_lengths, 'norm_phone_lengths', norm_phone_lengths);
    
end

phone_length_vec = cat(1, results.phone_lengths);

norm_phone_length_vec = cat(1, results.norm_phone_lengths);

phone_vec = cat(2, results.phones)';

%% Getting list of phonemes used in TIMIT.

[timit_phonemes, class_indicator, class_names] = getTIMITphones;

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

%%% Grouping by phoneme class.

[timit_map, ~] = find(timit2group_mat');
[class_map, ~] = find(class_indicator');

class_index = class_map(timit_map(group_index));

no_bins = ceil(sqrt(length(SI)));

vecs = {phone_length_vec, norm_phone_length_vec};
indices = {phone_index, class_index};
ids = {present_phonemes, class_names};
vec_labels = {'phone', 'normPhone'};
index_labels = {'', 'Class'};
no_skipped = [3 1];

for v = 1:length(vecs)
    for i = 1:length(indices)
        
        fname = [vec_labels{v}, index_labels{i}, name];
        [count, prob, stats, cdf, hist, bins, bin_centers] = calcStats(vecs{v}, indices{i}, ids{i}, no_bins, fname);
        
        figure()
        plot(1:length(ids{i}), prob/sum(prob), 'LineWidth', 2, 'Color', 'k')
        set(gca, 'XTick', 1:length(ids{i}), 'XTickLabel', ids{i})
        xtickangle(45)
        axis tight
        title('Phoneme Distribution')
        ylabel('Probability')
        
        saveas(gcf, [fname, '.fig'])
        
        plotHistograms(fname, ids{i}(1:(end - no_skipped(i))), bin_centers(1:(end - no_skipped(i)), :)', hist(1:(end - no_skipped(i)), :)')
        
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

function plotHistograms(fname, id, bin_centers, hist)

%% Plotting hists.

figure

x = cumsum(ones(length(id), size(bin_centers, 1)))';

pcolor(bin_centers, x, hist)

colormap('hot')

set(gca, 'YTick', 1:length(id), 'YTickLabel', id)

shading interp

title('Phoneme Duration Distribution')

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

saveas(gcf, [fname, '.fig'])

save_as_pdf(gcf, fname)

figure

colors = flipud(hsv(length(id)));

set(gca, 'NextPlot', 'add', 'ColorOrder', colors)

plot(bin_centers, hist, 'LineWidth', 2)

legend(id)

axis tight

title('Phoneme Duration Distribution')

saveas(gcf, [fname, '_line.fig'])

save_as_pdf(gcf, [fname, '_line'])

end