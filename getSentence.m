function sentence = getSentence(time, file_index, tsylb_option, no_channels, onset, transition_factor, lowfreq, highfreq, norm)

if nargin == 0, time = (0:.01:6000)/1000; fs = 10^5;
else, fs = length(time)/max(time); end
if nargin < 2, file_index = []; end
if nargin < 3, tsylb_option = []; end
if isempty(tsylb_option), tsylb_option = 1; end
if nargin < 4, no_channels = []; end
if isempty(no_channels), no_channels = 1; end
if nargin < 5, onset = []; end
if isempty(onset), onset = 0; end
if nargin < 6, transition_factor = []; end
if isempty(transition_factor), transition_factor = .1; end
if nargin < 7, lowfreq = []; end
if isempty(lowfreq), lowfreq = 0; end
if nargin < 8, highfreq = []; end
if isempty(highfreq), highfreq = inf; end
if nargin < 9, norm = []; end
if isempty(norm), norm = 0; end

nsltools_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/nsltools/'; addpath(genpath(nsltools_dir))
panphon_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/panphon/panphon/data/'; addpath(genpath(panphon_dir))
sentence_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

file_list_id = fopen([sentence_dir, 'wavFileList.txt'], 'r');
file_list = textscan(file_list_id, '%s');
fclose(file_list_id);

file_list = file_list{1};

if isempty(file_index)
    file_index = round(rand*length(file_list));
elseif isfloat(file_index) && file_index < 1 && file_index > 0
    file_index = round(file_index*length(file_list));
end

wavfile_name = file_list{file_index};
file_name = extractBefore(wavfile_name, '.WAV');

sentence = struct('filename', file_name);
sentence.time = time;

%% Getting phone to feature mapping.

feature_names = {'vowel', 'open', 'front', 'round', 'dipthong', 'rhoticized', 'place', 'manner', 'voice', 'affricate', 'syllabic'};
feature_dim = length(feature_names);
sentence.feature_names = feature_names;

if tsylb_option

    phone2feature_name = 'tsylbPhones2features.mat';
    phone_suffix = '.tsylbPHN';

else

    phone2feature_name = 'phones2features.mat';
    phone_suffix = '.PHN';

end

phone2feature_data = load(phone2feature_name);

unique_phones = phone2feature_data.unique_phones;
feature_mat = phone2feature_data.feature_mat;
features = phone2feature_data.features;

sentence.phone_list = unique_phones;
sentence.feature_mat = feature_mat;
sentence.features = features;

%% Getting list of phones/words and transition times.

phone_fid = fopen([sentence_dir, file_name, phone_suffix]);

phone_data = textscan(phone_fid, '%d%d%s');
phone_starts = round(100*phone_data{1}/16);
phone_ends = round(100*phone_data{2}/16);
phone_transition_times = time(sort(unique([phone_starts; phone_ends]))+1);
phones = phone_data{3};

num_phones = length(phones);

word_fid = fopen([sentence_dir, file_name, '.WRD']);

word_data = textscan(word_fid, '%d%d%s');
word_starts = round(100*word_data{1}/16);
word_ends = round(100*word_data{2}/16);
word_transition_times = time(sort(unique([word_starts; word_ends]))+1);
words = word_data{3};


%% Constructing time series of feature vectors.

this_phone_list = phones;
this_feature_mat = nan(num_phones, feature_dim);

feature_vec = zeros(length(time), feature_dim);

prev_transition = phone_transition_times(1);
prev_features = this_feature_mat(1, :);

phone_index = 1;

for p = 1:num_phones
    
    this_phone = phones{p};
    these_features = unique(features{strcmpi(unique_phones, this_phone)}, 'rows');
    this_feature_dim = size(these_features, 1);
    prev_feature_dim = size(prev_features, 1);
    
    this_start = phone_transition_times(p);
    this_end = phone_transition_times(p + 1);
    this_duration = this_end - this_start;
    
    max_transition_length = this_duration*transition_factor;
    
    first_transition = this_start + rand*max_transition_length;
    final_transition = this_end - rand*max_transition_length*(p ~= length(this_phone_list));
    
    transition_indicator = time > prev_transition & time <= first_transition;
    sustain_indicator = time > first_transition & time <= final_transition;

    if this_feature_dim == 1 && prev_feature_dim == 1

        if sum(transition_indicator) > 0

            transition_path = ((time(transition_indicator) - prev_transition)'*these_features + (first_transition - time(transition_indicator))'*prev_features)/(first_transition - prev_transition);
            feature_vec(transition_indicator, :) = transition_path;

        end

        feature_vec(sustain_indicator, :) = repmat(these_features, sum(sustain_indicator), 1);

    else

        [transition_f_choice, ~] = find(mnrnd(1, ones(1, this_feature_dim*prev_feature_dim)/(this_feature_dim*prev_feature_dim), sum(transition_indicator))');
        [sustain_f_choice, ~] = find(mnrnd(1, ones(1, this_feature_dim)/this_feature_dim, sum(sustain_indicator))');

        for f = 1:this_feature_dim

            for p = 1:prev_feature_dim

                this_choice = (f - 1)*prev_feature_dim + p;

                transition_path = ((time(transition_indicator) - prev_transition)'*these_features(f, :) + (first_transition - time(transition_indicator))'*prev_features(p, :))/(first_transition - prev_transition);
                
                transition_index = transition_indicator;
                transition_index(transition_indicator) = transition_f_choice == this_choice;
                
                feature_vec(transition_index, :) = transition_path(transition_f_choice == this_choice, :);

            end

            sustain_index = sustain_indicator;
            sustain_index(sustain_indicator) = sustain_f_choice == f;
            
            feature_vec(sustain_index, :) = repmat(these_features(f, :), sum(sustain_index), 1);

        end
        
    end

    prev_transition = final_transition;
    prev_features = these_features;
    
end

sentence.phone_transition_times = phone_transition_times;
sentence.phone_sequence = this_phone_list;
sentence.feature_vec = feature_vec;

%% Adding noise to time series of feature vectors.

noise_cov = 2*cov(zscore(double(feature_mat)))/3;
noise_mean = zeros(feature_dim, 1);
feature_noise = mvnrnd(noise_mean, noise_cov, length(time));

input_vec = feature_vec + feature_noise;

sentence.input_vec = input_vec;

%% Plotting noise covariance & time series of input vectors.

make_fig = 1;

if make_fig == 1
    
    % Noise covariance.
    figure, imagesc(noise_cov), colorbar
    xticks(1:11), xticklabels(feature_names)
    xtickangle(45)
    yticks(1:11), yticklabels(feature_names)
    axis xy
    
    % Input vector.
    figure

    subplot(211)
    
    imagesc(time, 1:feature_dim, feature_vec')%*diag(1./max(feature_mat)))')
    
    hold on
    
    plot(repmat(phone_transition_times, 2, 1), repmat([0; size(input_vec, 1)], 1, length(phone_transition_times)), 'w', 'LineWidth', .5)
    
    xticks(phone_transition_times(1:(end - 1))+diff(phone_transition_times)/2)
    xticklabels(this_phone_list)
    xtickangle(45)
    
    xlim([min(phone_transition_times), max(phone_transition_times)])
    
    yticks(1:feature_dim)
    yticklabels(feature_names)
    
    axis xy
    
    colorbar

    subplot(212)
    
    imagesc(time, 1:feature_dim, input_vec')%*diag(1./max(feature_mat)))')
    
    hold on
    
    plot(repmat(phone_transition_times, 2, 1), repmat([0; size(input_vec, 1)], 1, length(phone_transition_times)), 'w', 'LineWidth', .5)
    
    xticks(phone_transition_times(1:(end - 1))+diff(phone_transition_times)/2)
    xticklabels(this_phone_list)
    xtickangle(45)
    
    xlim([min(phone_transition_times), max(phone_transition_times)])
    
    yticks(1:feature_dim)
    yticklabels(feature_names)
    
    axis xy
    
    colorbar
    
end

% %% Getting syllable boundaries.
% 
% syl_boundaries = round(100*syl_boundaries/16);
% 
% syllable_indicator = zeros(size(cochlear, 1), 1);
% 
% syllable_indicator(onset*100 + syl_boundaries) = 1;
% 
% syllable_indicator = conv(syllable_indicator, ones(1,500)/500, 'same');
% 
% sentence.syllables = syllable_indicator;

% %% Getting wav file & applying cochlear filterbank.
% 
% this_wavfile = audioread([sentence_dir, wavfile_name]);
% wavfile_length = 100*length(this_wavfile)/16;
% 
% loadload, paras(1) = 1;% paras(4) = 0;
% this_auditory = wav2aud(this_wavfile, paras);
% 
% CF = 440 * 2 .^ ((-30:97)/24 - 1);
% CF_included = CF(CF >= lowfreq & CF <= highfreq);
% %CF_limits = [min(CF_included), max(CF_included)];
% this_auditory = this_auditory(:, CF >= lowfreq & CF <= highfreq);
% 
% [~, aud_cols] = size(this_auditory);
% 
% sentence_length = ceil(max(length(time), wavfile_length + onset*100));
% 
% cochlear = zeros(sentence_length, no_channels);
% 
% channels_averaged = floor(aud_cols/no_channels);
% for c = 1:no_channels
%     %this_channel_freqs(c, :) = CF_included([(c - 1)*channels_averaged + 1, c*channels_averaged]);
%     this_channel = nanmean(this_auditory(:,((c - 1)*channels_averaged + 1):c*channels_averaged), 2);
%     this_channel = interp(interp(this_channel, 5), 10);
%     if norm == 1
%         this_channel = this_channel/sum(this_channel);
%     elseif norm == 2
%         this_channel = this_channel/max(this_channel);
%     end
%     cochlear(onset*100 + (1:length(this_channel)), c) = conv(this_channel, ones(1,500)/500, 'same');  
% end
% 
% sentence.cochlear = cochlear;

end


function [format]=make_format(no_entries,entry_format)

format=[];
for i=1:no_entries
    format=[format,'%',entry_format,'\t'];
end
format=[format(1:end-1),'n'];

end

