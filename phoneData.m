function results = phoneData(SI, tsylb_flag)

if nargin < 1, SI = []; end
if isempty(SI), SI = (1:6300)/6300; end
if nargin < 2, tsylb_flag = []; end
if isempty(tsylb_flag), tsylb_flag = 0; end

time = 0:.1:10000;

onset_time = 1000;

name = 'Data';
if length(SI) ~= 6300
    name = sprintf('%s_%dsentences', name, length(SI));
end
  
timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

file_list_id = fopen([timit_dir, 'DOC/allphonelist_filenames.txt'], 'r');
file_list = textscan(file_list_id, '%s');
fclose(file_list_id);
    
file_list = file_list{1};

for s = 1:length(SI)
    
    file_index = SI(s);
    
    if isfloat(file_index) && file_index < 1 && file_index > 0
        file_index = round(file_index*length(file_list));
    end
    
    file_name = file_list{file_index};
    % file_name = extractBefore(file_name, '.WAV');
    
    %% Retrieving phonemes and their start and end times.
    
    if tsylb_flag
        phone_filename = [timit_dir, file_name, '.tsylbPHN'];
    else
        phone_filename = [timit_dir, file_name, '.PHN'];
    end
    fid = fopen(phone_filename, 'r');
    phone_data = textscan(fid, '%s');
    fclose(fid);
    phone_data = phone_data{1};
    phone_data = reshape(phone_data, 3, length(phone_data)/3);
    
    phones = phone_data(3, :);
    
    phone_indices = cellfun(@str2num, phone_data(1:2, :));
    phone_times = ((phone_indices/16 + onset_time)/1000)';
    
    phone_durations = diff(phone_times, [], 2); 
    
    %% Retrieving syllable boundary_times & normalizing phone lengths.
    
    sylb_filename = [timit_dir, file_name, '.SYLB'];
    fid = fopen(sylb_filename, 'r');
    tsylb_indices = textscan(fid, '%d %d %s');
    fclose(fid);
    tsylb_indices = tsylb_indices{1};
    syllable_times = (tsylb_indices/16 + onset_time)/1000;
    syllable_durations = diff(syllable_times);
    
    norm_phone_durations = phone_durations/mean(syllable_durations);
    
    %% Saving results.
    
    results(s) = struct('phone_durations', phone_durations, 'phones', {phones},...
        'syllable_durations', syllable_durations, 'norm_phone_durations', norm_phone_durations);
    
end

length_vec = cat(1, results.phone_durations);

norm_length_vec = cat(1, results.norm_phone_durations);

phone_vec = cat(2, results.phones)';

no_bins = ceil(sqrt(length(SI)));

%% Plotting overall phoneme duration distribution.

index = ones(size(length_vec));
id = {'Phoneme Duration'};
fname = ['phoneDuration', name, '_prob'];

[count, prob, stats, cdf, hist, bins, bin_centers, hist_cell, bins_cell, bin_centers_cell] = calcStats(length_vec, index, id, no_bins, fname);


plotIndividualizedHistograms(fname, id, bin_centers_cell', hist_cell')

%% Getting list of phonemes used in TIMIT.

[tsylb_phonemes, class_indicator, class_names] = getPhones(tsylb_flag); % The argument retrieves tsylb-compatible phonemes.

%% Collecting phone lengths across sentences.

%%% Grouping by phonemes.

[group_index, group_id] = findgroups(phone_vec);

%%% Converting from group_id order to timit_phonemes order.

get2group = cellfun(@(x) strcmp(x, group_id), tsylb_phonemes, 'unif', 0);
get2group_mat = cat(2, get2group{:});

present_phonemes = tsylb_phonemes;
if ~isempty(find(sum(get2group_mat') == 2))
    repeat_index = find(strcmp(present_phonemes, group_id(sum(get2group_mat') == 2)), 1, 'last');
else
    repeat_index = [];
end
zero_index = find(sum(get2group_mat) == 0);
present_phonemes([repeat_index, zero_index]) = [];
class_indicator([repeat_index, zero_index], :) = 0;

present2group = cellfun(@(x) strcmp(x, group_id), present_phonemes, 'unif', 0);
present2group_mat = cat(2, present2group{:});

[present_map, ~] = find(present2group_mat');
phone_index = present_map(group_index);

%%% Grouping by phoneme class.

%[tsylb_map, ~] = find(tsylb2group_mat');
[class_map, ~] = find(class_indicator');

class_index = class_map(present_map(group_index));

vecs = {length_vec, norm_length_vec};
indices = {phone_index, class_index};
ids = {present_phonemes, class_names};
if tsylb_flag
    vec_labels = {'tsylbPhone', 'normTsylbPhone'};
else
    vec_labels = {'phone', 'normPhone'};
end
index_labels = {'', 'Class'};
no_skipped = [3 1];
sort_option = [1 1];

for v = 1:length(vecs)
    for i = 1:length(indices)
        
        fname = [vec_labels{v}, index_labels{i}, name];
        [count, prob, stats, cdf, hist, bins, bin_centers, hist_cell, bins_cell, bin_centers_cell] = calcStats(vecs{v}, indices{i}, ids{i}, no_bins, fname);
        
        plotProbability(prob, ids{i}, sort_option(i))
        
        saveas(gcf, [fname, '.fig'])
        
        plotHistograms(fname, ids{i}(1:(end - no_skipped(i))), bin_centers(1:(end - no_skipped(i)), :)', hist(1:(end - no_skipped(i)), :)')
        
    end
end

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
end

if length(ids) > 27
    x_tick_step = round(length(ids)/27);
    these_x_ticks = 3:x_tick_step:length(ids);
    x_tick_labels = ids(these_x_ticks);
else
    x_ticks = 1:length(ids)
    x_tick_labels = ids;
end

figure()

plot(1:length(ids), prob/sum(prob), 'LineWidth', 2, 'Color', 'k')
set(gca, 'XTick', these_x_ticks, 'XTickLabel', x_tick_labels)
xtickangle(45)
axis tight
title('Phoneme Distribution')
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

title('Phoneme Duration Distribution')

saveas(gcf, [fname, '_line.fig'])

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