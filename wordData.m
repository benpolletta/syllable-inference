function results = wordData(SI)

if nargin < 1, SI = []; end
if isempty(SI), SI = (1:6300)/6300; end

time = 0:.1:10000;

onset_time = 1000;

name = 'Data';
if length(SI) ~= 6300
    name = sprintf('%s_%dsentences', name, length(SI));
end 

timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

timit_filenames_file = [timit_dir, 'DOC/allphonelist_filenames.txt'];
filename_fid = fopen(timit_filenames_file, 'r');
file_list = textscan(filename_fid, '%s', 'Delimiter', '\n');
fclose(filename_fid);
file_list = file_list{1};

%results(length(SI)) = struct();

for s = 1:length(SI)
    
    file_index = SI(s);
    
    if isfloat(file_index) && file_index < 1 && file_index > 0
        file_index = round(file_index*length(file_list));
    end
    
    file_name = file_list{file_index};
    
    %% Retrieving words and their start and end times.
    
    word_filename = [timit_dir, file_name, '.WRD'];
    fid = fopen(word_filename, 'r');
    word_data = textscan(fid, '%s');
    fclose(fid);
    word_data = word_data{1};
    word_data = reshape(word_data, 3, length(word_data)/3);
    
    words = word_data(3, :);
    
    word_indices = cellfun(@str2num, word_data(1:2, :))';
    word_times = (word_indices/16 + onset_time);
    
    word_lengths = diff(word_times, [], 2); 
    
    %% Retrieving phonemes and their start and end times.
    
    phone_filename = [sentence_dir, file_name, '.PHN'];
    fid = fopen(phone_filename, 'r');
    phone_data = textscan(fid, '%s');
    fclose(fid);
    phone_data = phone_data{1};
    phone_data = reshape(phone_data, 3, length(phone_data)/3);
    
    phones = phone_data(3, :);
    
    phone_indices = cellfun(@str2num, phone_data(1:2, :))';
    %mid_phone_indices = mean(phone_indices);
    
    phone_times = (phone_indices/16)';
    phone_lengths = diff(phone_times, [], 2);
    
    %% Retrieving & writing pronunciations.

    pron_filename = [sentence_dir, file_name, '.WRDpron'];
    fid = fopen(pron_filename, 'w');

    pronunciations = cell(size(word_lengths));
    
    for w = 1:length(word_lengths)
        
        this_word_indicator = any(phone_indices > word_indices(w, 1) & phone_indices < word_indices(w, 2), 2);
        
        this_pronunciation = phones(this_word_indicator);

        pronunciations{w} = {this_pronunciation};

        format = join(repmat({'%s\t'}, 1, length(this_pronunciation)), '');
        format = [format{1}, '\n'];

        fprintf(fid, format, this_pronunciation{:})
        
    end

    fclose(fid)
    
    %% Retrieving syllable boundary_times & normalizing word lengths.
    
    sylb_filename = [sentence_dir, file_name, '.SYLB'];
    fid = fopen(sylb_filename, 'r');
    sylb_indices = textscan(fid, '%d %d %s');
    fclose(fid);
    sylb_indices = sylb_indices{1};
    syllable_lengths = diff(sylb_indices/16);
    mid_syl_indices = sylb_indices(1:(end - 1)) + diff(sylb_indices)/2;
    
    norm_word_lengths = word_lengths/mean(syllable_lengths);
    
    syllabifications = cell(size(word_lengths));
    
    for w = 1:length(word_lengths)
        
        this_word_indicator = mid_syl_indices > word_indices(w, 1) & mid_syl_indices < word_indices(w, 2);
        
        this_word_index = find(this_word_indicator);
        
        num_syllables(w) = length(this_word_index);
        
        syllables = cell(num_syllables(w), 1);
        
        for syl = 1:num_syllables(w)
            
            syl_index = this_word_index(syl);
            
            this_syl_indicator = any(phone_indices > sylb_indices(syl_index) & phone_indices < sylb_indices(syl_index + 1), 2);
        
            syllables{syl} = phones(this_syl_indicator);
            
        end
        
        syllabifications{w} = {syllables};
        
    end
    
    %% Saving results.
    
    results(s) = struct('word_lengths', word_lengths, 'words', {words}, 'pronunciations', {pronunciations},...
        'syllable_lengths', syllable_lengths, 'norm_word_lengths', norm_word_lengths, 'syllabifications', {syllabifications});
    
end

word_length_vec = cat(1, results.word_lengths);

norm_word_length_vec = cat(1, results.norm_word_lengths);

word_vec = cat(2, results.words)';

pronunciation_vec_cell = cat(1, results.pronunciations);
pronunciation_vec_string = cellfun(@(x) strjoin(x{1}, '/'), pronunciation_vec_cell, 'unif', 0);


%% Collecting word lengths across sentences.

%%% Grouping by words.

[word_index, word_id] = findgroups(word_vec);

% Finding number of pronunciations.

pronunciation_map_cell = splitapply(@(x) {x}, pronunciation_vec_cell, word_index);
pronunciation_map_string = splitapply(@(x) {x}, pronunciation_vec_string, word_index);
[pronunciation_index, unique_pronunciations] = cellfun(@findgroups, pronunciation_map_string, 'unif', 0);
num_pronunciations = cellfun(@length, unique_pronunciations);

%%% Grouping by phoneme class.

[timit_map, ~] = find(timit2group_mat');
[class_map, ~] = find(class_indicator');

class_index = class_map(timit_map(word_index));

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