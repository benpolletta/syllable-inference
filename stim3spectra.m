function stim3spectra %(band_indices)
% band_indices: vector of integers between 1 and 128, indicating which
% cochlear bands to include in the analysis.

speeds = [17.5 16.0127 14.0343 11.0707 5];
fs = 44.1*10^3;

if exist('stimlist3spec_data.mat') ~=2
    
    time = 0.01:0.01:6000;
    num_channels = 16;
    channels = getCochlearBands(num_channels)';
    onset = 1000;
    norm = 0;
    
    num_bands = 1;
    bands = getCochlearBands(num_bands)';
    % bands = cochlearBands(:,  band_indices); num_bands = length(bands);
    
    nsltools_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/nsltools/'; addpath(genpath(nsltools_dir))
    sentence_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/stimlist3/';
    
    file_list_id = fopen([sentence_dir, 'wavFileList.txt'], 'r');
    file_list = textscan(file_list_id, '%s');
    fclose(file_list_id);
    file_list = file_list{1}; num_files = length(file_list);
    
    no_speeds = num_files/60;
    pow_sum = nan(60, no_speeds, num_bands);
    
    speed_indices = ceil((1:num_files)/60);
    stim_indices = mod(1:num_files, 60);
            
    sentence = getStimulusForSpec(3, num_channels, 1/num_files, .5, bands(1, 1), bands(2, 1), norm);
    size1 = size(sentence.cochlear, 1);
    
    [~, freqs1] = pmtm(sentence.cochlear, [], [], fs); % pmtm(decimate(decimate(sentence.cochlear, 10), 10), [], [], 10^3);
    max_freq = 100;
    freqs = freqs1(freqs1 <= max_freq);
    
    all_pow = nan(num_files, num_channels);
    all_spec = nan(sum(freqs1 <= max_freq), num_files, num_channels);
    
    for f = 1:20:num_files
        
        sentence = getStimulusForSpec(3, num_channels, f/num_files, .5, bands(1, 1), bands(2, 1), norm);
        
        tic
        
        all_pow(f, :) = sqrt(sum(abs(sentence.cochlear).^2));
        
        toc
        
        tic
        
        this_size = size(sentence.cochlear, 1);
        
        if this_size == size1
            
            for c = 1:num_channels
                
                [this_pmtm, this_freqs] = pmtm(abs(sentence.cochlear(:, c)), [], [], fs); % pmtm(decimate(decimate(sentence.cochlear(:, c), 10), 10), [], [], 10^3);
                
                all_spec(:, f, c) = this_pmtm(this_freqs <= max_freq);
                
            end
            
%         else
%             
%             for c = 1:num_channels
%             
%                 sentence.cochlear = decimate(sentence.cochlear(:, c), this_size/size1);
%                 
%                 [this_pmtm, this_freqs] = pmtm(decimate(decimate(sentence.cochlear(:, c), 10), 10), [], [], 10^3);
%                 
%                 all_spec(:, f, c) = this_pmtm(this_freqs <= max_freq);   
%                 
%             end
            
        end
        
        toc
        
    end
    
    all_norm_pow = nanunitsum(all_pow, 2);
    
    all_pow = reshape(all_pow, [60, no_speeds, num_channels]);
    all_norm_pow = reshape(all_norm_pow, [60, no_speeds, num_channels]);
    
    all_spec = reshape(all_spec, [length(freqs), 60, no_speeds, num_channels]);
    
    mean_pow = squeeze(nanmean(all_pow));
    
    std_pow = squeeze(nanstd(all_pow));
    
    mean_norm_pow = squeeze(nanmean(all_norm_pow));
    
    std_norm_pow = squeeze(nanstd(all_norm_pow));
    
    mean_spec = squeeze(nanmean(all_spec, 2));
    
    std_spec = squeeze(nanstd(all_spec, [], 2));
    
    save('stimlist3spec_data.mat', 'all_pow', 'all_norm_pow', 'all_spec',...
        'mean_pow', 'std_pow', 'mean_norm_pow', 'std_norm_pow', 'mean_spec', 'std_spec', 'freqs', 'channels')

else
    
    spec_data = load('stimlist3spec_data.mat');
    
    channels = spec_data.channels;
    freqs = spec_data.freqs;
    
    mean_pow = spec_data.mean_pow;
    std_pow = spec_data.std_pow;
    mean_norm_pow = spec_data.mean_norm_pow;
    std_norm_pow = spec_data.std_norm_pow;
    mean_spec = spec_data.mean_spec;
    std_spec = spec_data.std_spec;
    
end  

figure

size_std_pow = size(std_pow);
bounds = nan(size_std_pow(1), 2, size_std_pow(2));

for i = 1:2
    bounds(:, i, :) = reshape(std_pow, [size_std_pow(1), 1, size_std_pow(2)]);
end

speed_labels = cellfun(@num2str, mat2cell(speeds, 1, ones(1, size(speeds, 2))), 'UniformOutput', 0);

subplot(2, 1, 1)
    
boundedline(mean(channels)', mean_pow', permute(bounds, [3, 2, 1]))
legend(speed_labels)
set(gca, 'XScale', 'log')
axis tight

for i = 1:2
    bounds(:, i, :) = reshape(std_norm_pow, [size_std_pow(1), 1, size_std_pow(2)]);
end

subplot(2, 1, 2)
    
boundedline(mean(channels)', mean_norm_pow', permute(bounds, [3, 2, 1]))
legend(speed_labels)
set(gca, 'XScale', 'log')
axis tight

saveas(gcf, 'stimlist3_mean_pow.fig')

size_std = size(std_spec);
bounds = nan(size_std(1), 2, size_std(2), size_std(3));

for i = 1:2
    bounds(:, i, :, :) = reshape(std_spec, [size_std(1), 1, size_std(2:3)]);
end

figure

channel_labels = cellfun(@num2str, mat2cell(mean(channels), 1, ones(1, size(channels, 2))), 'UniformOutput', 0);

for s = 1:length(speeds)
   
    subplot(length(speeds), 2, s)
    
    imagesc(freqs, channels, squeeze(mean_spec(:, s, :)))
    %boundedline(freqs, squeeze(mean_spec(:, s, :)), squeeze(bounds(:, :, s, :)))
    
    title(sprintf('Speed = %0.2g', speeds(s)))
    
    %legend(channel_labels)
    
end

saveas(gcf, 'stimlist3_mean_spec.fig')

end