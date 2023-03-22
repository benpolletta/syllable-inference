function results = sylbData(SI, reuse_results)

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
%dict_words = cellfun(@(x) erase(x, '-'), dict_words, 'unif', 0);
dict_sylbs = {dict(:).sylb_cell};

word_results = load('wordData_results.mat');
word_results = word_results.results;
canonical_words = {word_results(:).canonical_words};

ncsn_fid = fopen([timit_dir, 'DOC/non_canonical_sylb_num.txt'], 'w');

if reuse_results & exist('sylbData_results.mat') == 2

    sylb_data = load('sylbData_results.mat');

    results = sylb_data.results;

else

    for s = 1:length(SI)

        file_index = SI(s);

        if isfloat(file_index) && file_index < 1 && file_index > 0
            file_index = round(file_index*length(file_list));
        end

        file_name = file_list{file_index};

        %% Retrieving words, phonemes, and syllables and their start and end times.

        these_cw = canonical_words{s};

        word_filename = [timit_dir, file_name, '.WRD'];
        [words, word_times, word_durations] = getUnits(word_filename);

        phone_filename = [timit_dir, file_name, '.tsylbPHN'];
        [phones, phone_times, phone_durations] = getUnits(phone_filename);

        sylb_filename = [timit_dir, file_name, '.SYLB'];
        [sylbs, sylb_times, sylb_durations] = getUnits(sylb_filename);

        syllabic_rate = 1000/mean(sylb_durations);

        word_sylbs = wsp_map(s).word_sylb_cell;

        canonicalSylbs = cell(size(word_sylbs));

        for w = 1:length(words)

            %% Finding canonical pronunciations.

            dictIndex = strcmp(these_cw{w}, dict_words);
            dictSylbs = dict_sylbs{dictIndex};

            if length(dictSylbs) == length(word_sylbs{w})

                canonicalSylbs{w} = dictSylbs;

            elseif all(strcmp(dictSylbs, 'y/ih/r'))

                canonicalSylbs{w} = word_sylbs{w};
                
            else

                fprintf(ncsn_fid, '%s %s %s\n', file_name, strjoin(dictSylbs(:), '*'), strjoin(word_sylbs{w}, '*'));

            end

        end

        canonicalSylbs = cat(2, canonicalSylbs{:})';

        %% Saving results.

        results(s) = struct('sylb_durations', sylb_durations, 'syllabic_rate', syllabic_rate,...
            'phone_durations', phone_durations, 'canonicalSylbs', {canonicalSylbs}, 'sylbs', {sylbs}, 'phones', {phones});

    end

    save('sylbData_results.mat', 'results')

end

fclose(ncsn_fid);

%% Collecting syllable lengths across sentences.

all_sylbs = cat(1, results.sylbs);

phoneno_vec = cellfun(@length, all_sylbs);

firstphone_vec = cellfun(@(x) x{1}, all_sylbs, 'unif', 0);

string_vec = cellfun(@(x) strjoin(x, '/'), all_sylbs, 'unif', 0);

duration_vec = cat(1, results.sylb_lengths);

norm_duration_vec = cat(1, results.norm_sylb_lengths);

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

function plotHistograms(name, class_names, hist_bin_centers, hist_counts)

%% Plotting histograms.

figure

x = cumsum(ones(length(class_names), size(hist_bin_centers, 1)))';

pcolor(hist_bin_centers, x, hist_counts)

colormap('hot')

set(gca, 'YTick', 1:length(class_names), 'YTickLabel', class_names)

shading interp

title('Syllable Duration Distribution')

nochange_colorbar(gca)

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