function results = wordData(SI, reuse_results)

global onset_time

if nargin < 1, SI = []; end
if isempty(SI), SI = (1:6300)/6300; end
if nargin < 2, reuse_results = []; end
if isempty(reuse_results), reuse_results = true; end

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

wsp_map = load('word2sylb2phone_bysentence.mat');
wsp_map = wsp_map.results;

dict = load('word2sylb2phone_DICT.mat');
dict = dict.results;
dict_words = {dict(:).word};
% dict_words = cellfun(@(x) erase(x, '-'), dict_words, 'unif', 0);

if reuse_results & exist('wordData_results.mat') == 2
    
    word_data = load('wordData_results.mat');

    results = word_data.results;
    % wsp_map = word_data.wsp_map;

else

%     [ambiguous_words, ambiguous_sentences, resolutions] = deal({});

    for s = 1:length(SI)

        file_index = SI(s);

        if isfloat(file_index) && file_index < 1 && file_index > 0
            file_index = round(file_index*length(file_list));
        end

        file_name = file_list{file_index};

        %% Retrieving words, phonemes, and syllables and their start and end times.

        word_filename = [timit_dir, file_name, '.WRD'];
        [words, word_times, word_durations] = getUnits(word_filename);

        phone_filename = [timit_dir, file_name, '.tsylbPHN'];
        [phones, phone_times, phone_durations] = getUnits(phone_filename);

        sylb_filename = [timit_dir, file_name, '.SYLB'];
        [syllables, sylb_times, sylb_durations] = getUnits(sylb_filename);

        %% Normalizing word lengths.

        norm_word_durations = double(word_durations)/mean(double(sylb_durations));

        %% Writing pronunciations.

        pron_filename = [timit_dir, file_name, '.WRDpron'];
        fid = fopen(pron_filename, 'w');

        these_pronunciations = wsp_map(file_index).word_cell;
        canonical_words = wsp_map(file_index).canonical_word_cell;

        for_print = mat2cell(word_times, ones(size(word_times, 1), 1), [1 1]);
        for_print(:, end + 1) = these_pronunciations(:);
        for_print = for_print';

        fprintf(fid, '%d %d %s\n', for_print{:});

        fclose(fid);

        %% Saving results.

        results(s) = struct('word_durations', word_durations, 'words', {words}, 'canonical_words', {canonical_words}, 'pronunciations', {these_pronunciations},...
            'sylb_durations', sylb_durations, 'norm_word_durations', norm_word_durations); % 'syllabifications', {syllabifications});

    end

    save('wordData_results.mat', 'results') % , 'ambiguous_words', 'ambiguous_sentences', 'resolutions') %, 'wsp_map')

end

word_length_vec = double(cat(1, results.word_durations));

norm_word_length_vec = double(cat(1, results.norm_word_durations));

word_sylb_num_vec = cat(2, wsp_map.word_sylb_num)';

word_phone_num_vec = cat(2, wsp_map.word_phone_num)';

word_vec = cat(1, results.words);

cword_vec = cat(1, results.canonical_words);

pron_vec = cat(2, wsp_map.word_cell)';

pair_vec = cellfun(@(x, y) strjoin({x,y}, '='), cword_vec, pron_vec, 'unif', 0);

no_bins = ceil(sqrt(length(SI)));

%% Plotting word duration distribution.

index = ones(size(word_length_vec));
id = {'Word Duration'};
fname = ['wordDuration', name, '_prob'];

[count, prob, stats, cdf, hist, bins, bin_centers, hist_cell, bins_cell, bin_centers_cell] = calcStats(word_length_vec, index, id, no_bins, fname);

plotIndividualizedHistograms(fname, id, bin_centers_cell', hist_cell')

%% Collecting word lengths across sentences.

%%% Calculating grouping variables.

% Grouping by words.

[word_index, word_id] = findgroups(word_vec);
[canonical_word_index, canonical_word_id] = findgroups(cword_vec);
[pron_index, pron_id] = findgroups(pron_vec);
[pair_index, pair_id] = findgroups(pair_vec);

% Grouping by length in syllables.

sylb_num_id = mat2cell(min(word_sylb_num_vec):max(word_sylb_num_vec), 1, ones(1, range(word_sylb_num_vec) + 1));
sylb_num_id = cellfun(@num2str, sylb_num_id, 'unif', 0);

% Grouping by length in phones.

phone_num_id = mat2cell(min(word_phone_num_vec):max(word_phone_num_vec), 1, ones(1, range(word_phone_num_vec) + 1));
phone_num_id = cellfun(@num2str, phone_num_id, 'unif', 0);

%% Finding number of pronunciations per word, & distribution over pronunciations.

%pron_map_cell = splitapply(@(x) {x}, pronunciation_vec_cell, word_index);
pronunciation_map = splitapply(@(x) {x}, pron_vec, canonical_word_index);
[pronunciation_index, unique_pronunciations] = cellfun(@findgroups, pronunciation_map, 'unif', 0);
num_pronunciations = cellfun(@length, unique_pronunciations);

save('wordPronunciations.mat', 'pronunciation_map', 'pronunciation_index', 'unique_pronunciations', 'num_pronunciations')

index = ones(size(num_pronunciations));
id = {'Number of Pronunciations'};
fname = ['numPronunciations', name, '_prob'];

[count, prob, stats, cdf, hist, bins, bin_centers, hist_cell, bins_cell, bin_centers_cell] = calcStats(num_pronunciations, index, id, no_bins, fname);

plotIndividualizedHistograms(fname, id, bin_centers_cell', hist_cell')

%%% Computing_stats.

vecs = {word_length_vec, norm_word_length_vec};
indices = {word_index, canonical_word_index, pron_index, pair_index, word_sylb_num_vec + 1, word_phone_num_vec};
ids = {word_id, canonical_word_id, pron_id, pair_id, sylb_num_id, phone_num_id};
sort_option = {1, 1, 1, 1, 0, 0};
vec_labels = {'word', 'normWord'};
index_labels = {'', 'Canon', 'Pron', 'Pair', 'NumSylbs', 'NumPhones'};

for v = 1:length(vecs)
    for i = 1:length(indices)
        
        fname = [vec_labels{v}, index_labels{i}, name];
        
        [count, prob, stats, cdf, hist, bins, bin_centers, hist_cell, bins_cell, bin_centers_cell] = calcStats(vecs{v}, indices{i}, ids{i}, no_bins, fname);
        
        plotProbability(prob, ids{i}, sort_option{i})
        
        
        saveas(gcf, [fname, '_prob.fig'])

        % save_as_pdf(gcf, [fname, '_prob'])
        
        plotHistograms(fname, ids{i}, bin_centers', hist')

        plotIndividualizedHistograms(fname, ids{i}, bin_centers_cell, hist_cell)
        
    end
end

end

function [units, times, lengths] = getUnits(filename)

global onset_time

fid = fopen(filename, 'r');
data = textscan(fid, '%d%d%s'); % '%s');
fclose(fid);
% data = data{1};
% data = reshape(data, 3, length(data)/3);

units = data{3}; % (3, :);

indices = cell2mat(data(1:2)); % cellfun(@str2num, data(1:2, :))';
times = (indices/16 + onset_time);

lengths = diff(times, [], 2);

end

function [count, prob, stats, cdf, hist, bins, bin_centers, hist_cell, bins_cell, bin_centers_cell] = calcStats(vec, index, id, no_bins, fname)

count = splitapply(@(x) length(x), vec, index);
prob = count/sum(count);
stats = splitapply(@(x) [mean(x), std(x), quantile(x, [.5, .25, .75])], vec, index);
cdf = splitapply(@(x) {[sort(x), linspace(0, 1, length(x))']}, vec, index);

[hist, bins] = splitapply(@(x) histcounts(x, no_bins, 'Normalization', 'probability'), vec, index); % no_bins, 'Normalization', 'probability'), vec, index);

bin_centers = bins(:, 1:(end - 1)) + diff(bins, [], 2)/2;
bin_centers = [bin_centers(:, 1) - mean(diff(bins, [], 2), 2),...
    bin_centers, bin_centers(:, end) + mean(diff(bins, [], 2), 2)];

hist = [zeros(size(hist(:, 1))), hist, zeros(size(hist(:, end)))];

[hist_cell, bins_cell] = arrayfun(@(i) histcounts(vec(index == i), 'Normalization', 'probability', 'BinMethod', 'sqrt'), min(index):max(index), 'UniformOutput', false);

bin_centers_cell = cellfun(@(x) x(:, 1:(end - 1)) + diff(x, [], 2)/2, bins_cell, 'unif', 0);
bin_centers_cell = cellfun(@(x,y) [x(:, 1) - mean(diff(y, [], 2), 2), x, x(:, end) + mean(diff(y, [], 2), 2)], bin_centers_cell, bins_cell, 'unif', 0);

hist_cell = cellfun(@(x) [zeros(size(x(:, 1))), x, zeros(size(x(:, end)))], hist_cell, 'unif', 0);

save([fname, '.mat'], 'id', 'count', 'prob', 'stats', 'cdf', 'hist', 'bins', 'bin_centers', 'hist_cell', 'bins_cell', 'bin_centers_cell')

end

function plotProbability(prob, ids, sort_option)

if sort_option
    [prob, sort_order] = sort(prob, 'descend');
    ids = ids(sort_order);

    these_x_ticks = int32(logspace(0, log10(length(ids)), 20));
else
    x_tick_step = round(length(ids)/20);
    these_x_ticks = 1:x_tick_step:length(ids);
end
x_tick_labels = ids(these_x_ticks);

figure()

if sort_option
    loglog(1:length(ids), prob/sum(prob), 'LineWidth', 2, 'Color', 'k')
else
    plot(1:length(ids), prob/sum(prob), 'LineWidth', 2, 'Color', 'k')
end
set(gca, 'XTick', these_x_ticks, 'XTickLabel', x_tick_labels)
xtickangle(45)
axis tight
title('Word Distribution')
ylabel('Probability')

end

function plotIndividualizedHistograms(fname, id, bin_centers, hist)

figure
hold on

colors = flipud(hsv(length(id)));

set(gca, 'NextPlot', 'add', 'ColorOrder', colors)

for i = 1:length(id)

    plot(bin_centers{i}, hist{i}, 'LineWidth', 2)

end

legend(id)

axis tight

title('Word Duration Distribution')

saveas(gcf, [fname, '_line_ind.fig'])

end

function plotHistograms(fname, id, bin_centers, hist)

%% Plotting hists.

figure

x = cumsum(ones(length(id), size(bin_centers, 1)))';

pcolor(bin_centers, x, hist)

colormap('hot')

set(gca, 'YTick', 1:length(id), 'YTickLabel', id)

shading interp

title('Word Duration Distribution')

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

saveas(gcf, [fname, '_duration.fig'])

% save_as_pdf(gcf, [fname, '_duration'])

if length(id) < 100

    figure

    colors = flipud(hsv(length(id)));

    set(gca, 'NextPlot', 'add', 'ColorOrder', colors)

    plot(bin_centers, hist, 'LineWidth', 2)

    legend(id)

    axis tight

    title('Word Duration Distribution')

    saveas(gcf, [fname, '_line.fig'])

    % save_as_pdf(gcf, [fname, '_line'])

end

end