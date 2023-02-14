function results = phoneSpectrograms(SI)

name = sprintf('Spectrograms_%dsentences', length(SI));

cochlear_freqs = 440 * 2 .^ ((-30:97)/24 - 1);
% sampling_freq = 16*10^3;

% if exist([name, '.mat'], 'file') == 2
%     
%     results = load([name, '.mat']);
%     
%     results = results.results;
%     
%     mean_phone_spect = results.mean_phone_spect;
%     
% else
    
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
        
        %% Retrieving speech spectrogram.
        
        if exist([sentence_dir, file_name, '_wt.mat'], 'file') == 2
            
            sentence_spectrogram = load([sentence_dir, file_name, '_wt.mat']);
            sentence_spectrogram = sentence_spectrogram.ws;
            % freqs = sentence_spectrogram.freqs;
            
        else
            
            this_wavfile = audioread([sentence_dir, wavfile_name]);
            % wavfile_length = 100*length(this_wavfile)/16;
            
            loadload, paras(1) = 1;% paras(4) = 0;
            sentence_spectrogram = wav2aud(this_wavfile, paras);
            save([sentence_dir, file_name, '_aud.mat'], 'sentence_spectrogram', 'this_wavfile', 'cochlear_freqs', 'sampling_freq');
            
            % no_cycles = linspace(7, 49, length(cochlear_freqs));
            % sentence_spectrogram = wavelet_spectrogram(this_wavfile, sampling_freq, cochlear_freqs, no_cycles, 1, [sentence_dir, file_name, '_wt.mat']);
            
        end
        
        %% Retrieving phonemes and their start and end times.
        
        phone_filename = [sentence_dir, file_name, '.tsylbPHN'];
        fid = fopen(phone_filename, 'r');
        phone_data = textscan(fid, '%s');
        fclose(fid);
        phone_data = phone_data{1};
        phone_data = reshape(phone_data, 3, length(phone_data)/3);
        
        phones = phone_data(3, :);
        
        phone_indices = cellfun(@str2num, phone_data(1:2, :));
        phone_spec_indices = max(1, round(phone_indices/(length(this_wavfile)/size(sentence_spectrogram, 1))));
        phone_times = (phone_indices/16)';
        
        phone_lengths = diff(phone_times, [], 2);
        
        %% Computing wavelet transform.
        
        no_phones = length(phones);
        
        phone_spect = cell(no_phones, 1);
        
        for p = 1:length(phones)
            
            phone_spect{p} = sentence_spectrogram(phone_spec_indices(1, p):phone_spec_indices(2, p), :);
            
        end
        
        % plotSpectrograms(cochlear_freqs, phones, phone_spect)
        
        %% Saving results.
        
        results(s) = struct('phone_lengths', phone_lengths, 'phone_spect', {phone_spect}, 'phones', {phones});
        
    end
    
    phone_length_vec = cat(1, results.phone_lengths);
    
    phone_vec = cat(2, results.phones)';
    
    phone_spect_vec = cat(1, results.phone_spect);
    
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
    
    %%% Finding mean and plotting.
    
    indexes = {phone_index, class_index};
    
    index_labels = {present_phonemes, class_names};
    
    methods = {'Fill', 'Stretch'};
    
    statistics = {'Mean', 'CV'};
    
    mean_phone_spect = cell(length(indexes), length(methods), length(statistics));
    
    for i = 1:length(indexes)
        
        for m = 1:length(methods)
            
            for s = 1:length(statistics)
                
                %%% Finding means across phonemes.
                
                mean_phone_spect{i, m, s} = splitapply(@(x) {avgSpectrogram(x, methods{m}, statistics{s})}, phone_spect_vec, indexes{i});
                
            end
            
        end
        
    end
    
    save([name, '.mat'], '-v7.3', 'results', 'mean_phone_spect', 'indexes', 'index_labels')
    
% end

%% Plotting.

index_names = {'phone', 'class'};

methods = {'Fill', 'Stretch'};

statistics = {'Mean', 'CV'};

colorbar_flags = {false, true};

for i = 1:length(indexes)
    
    for m = 1:length(methods)
        
        for s = 1:length(statistics)
            
            plotSpectrograms(cochlear_freqs, index_labels{i}, mean_phone_spect{i, m, s},...
                [index_names{i}, methods{m}, statistics{s}, name], colorbar_flags{s})
            
        end
        
    end
    
end

end

function plotSpectrograms(cochlear_freqs, phones, phone_spect, name, colorbar_flag)

if nargin < 4, name = []; end

no_phones = length(phone_spect);

phone_lengths = cellfun(@(x) size(x, 1), phone_spect);

figure

[rows, cols] = subplot_size(no_phones);

ax = tight_subplot(rows, cols);

for p = 1:no_phones
    
    axes(ax(p))
    
    imagesc(1:phone_lengths(p), 1:length(cochlear_freqs), phone_spect{p}')
    
    set(gca, 'YTick', 1:16:length(cochlear_freqs), 'YTickLabels', round(cochlear_freqs(1:16:end)/1000, 2))
    
    axis xy
    
    xlim([0, min(phone_lengths(p), 500)])
    
    if colorbar_flag
    
        nochange_colorbar(gca)
        
    end
    
    title(phones{p})
    
end

if ~isempty(name)
    
    saveas(gcf, [name, '.fig'])
    
    save_as_pdf(gcf, name)
    
end

end

function spect_stat = avgSpectrogram(phone_spect, method, statistic)

if nargin < 2, method = []; end
if isempty(method), method = 'fill'; end
if nargin < 3, statistic = []; end
if isempty(statistic), statistic = 'mean'; end

switch method
    
    case 'Fill'
        
        phone_length = cellfun(@(x) size(x, 1), phone_spect);
        
        expanded_phone_spect = cell(size(phone_spect));
        
        expanded_phone_spect(:) = {zeros(max(phone_length), size(phone_spect{1}, 2))};
        
        for p = 1:length(expanded_phone_spect)
            
            [r, c] = size(phone_spect{p});
            
            expanded_phone_spect{p}(1:r, 1:c) = phone_spect{p};
            
        end

    case 'Stretch'
        
        spect_xy = cellfun(@(x) {1:size(x, 1), 1:size(x,2)}, phone_spect, 'unif', false);
        
        spect_interp = cellfun(@(x, y) griddedInterpolant(x, y), spect_xy, phone_spect, 'unif', false);
        
        norm_xy = cellfun(@(x) {1:(size(x,1) - 1)/99:size(x, 1), 1:(size(x,2) - 1)/99:size(x,2)}, phone_spect, 'unif', false);
        
        expanded_phone_spect = cellfun(@(x, y) x(y), spect_interp, norm_xy, 'unif', false);
        
end

switch statistic
    
    case 'Mean'
        
        spect_stat = nanmean(cat(3, expanded_phone_spect{:}), 3);
        
    case 'STD'
        
        spect_stat = nanstd(cat(3, expanded_phone_spect{:}), [], 3);
        
    case 'CV'
        
        spect_mean = nanmean(cat(3, expanded_phone_spect{:}), 3);
        
        spect_std = nanstd(cat(3, expanded_phone_spect{:}), [], 3);
        
        spect_stat = spect_std./spect_mean;
        
end

end