function sentence = getSentence(time, no_channels, file_index, onset, lowfreq, highfreq, norm)

if nargin < 2, no_channels = []; end
if isempty(no_channels), no_channels = 1; end
if nargin < 3, file_index = []; end
if nargin < 4, onset = []; end
if isempty(onset), onset = 0; end
if nargin < 5, lowfreq = []; end
if isempty(lowfreq), lowfreq = 0; end
if nargin < 5, highfreq = []; end
if isempty(highfreq), highfreq = inf; end
if nargin < 6, norm = []; end
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

%% Getting phone to feature mapping.

phone2feature_fid = fopen([panphon_dir, 'timit_features.txt']);
format = '%s';
for i = 1:24, format = [format, '%d']; end
phone2feature_data = textscan(phone2feature_fid, format);
phone_list = phone2feature_data{1};
feature_mat = [phone2feature_data{2:end}];
feature_dim = size(feature_mat, 2);

unique_phones = unique(phone_list);

for u = 1:length(unique_phones)
    features{u} = feature_mat(strcmpi(phone_list, unique_phones{u}), :);
end

%% Getting list of phones and transition times.

phone_fid = fopen([sentence_dir, file_name, '.PHN']);

phone_data = textscan(phone_fid, '%d%d%s');
phone_starts = round(100*phone_data{1}/16);
phone_ends = round(100*phone_data{2}/16);
phone_transition_times = time(sort(unique([phone_starts; phone_ends]))+1);
phones = phone_data{3};

num_phones = length(phones);

%% Constructing specific feature & transition matrices for this sentence (adding extra transitions for paired feature vectors).

this_phone_list = phones;
this_feature_mat = nan(num_phones, feature_dim);

phone_index = 1;

for p = 1:num_phones
    
    this_phone = phones{p};
    this_feature = features{strcmpi(unique_phones, this_phone)};
    this_feature_dim = size(this_feature, 1);
    
    if this_feature_dim == 1
        
        this_feature_mat(phone_index, :) = this_feature;
        
    elseif this_feature_dim == 2
        
        this_phone_list((phone_index + 2):(end + 1)) = this_phone_list((phone_index + 1):end);
        this_phone_list(phone_index:(phone_index + 1)) = {this_phone(1), this_phone(2)};
        
        phone_transition_times((phone_index + 2):(end + 1)) = phone_transition_times((phone_index + 1):end);
        phone_transition_times(phone_index + 1) = mean(phone_transition_times(phone_index:(phone_index + 1)));
        
        this_feature_mat(phone_index:(phone_index + 1), :) = this_feature;
        
    end
    
    phone_index = phone_index + this_feature_dim;
    
end

%% Constructing time series of feature vectors.

feature_vec = nan(length(time), feature_dim);

prev_transition = phone_transition_times(1);
prev_feature = this_feature_mat(1, :);

for p = 1:length(this_phone_list)
    
    %this_phone = this_phone_list{p};
    this_feature = this_feature_mat(p, :);
    this_start = phone_transition_times(p);
    this_end = phone_transition_times(p + 1);
    this_duration = this_end - this_start;
    
    max_transition_length = this_duration/2;
    
    first_transition = this_start + rand*max_transition_length;
    final_transition = this_end - rand*max_transition_length*(p ~= length(this_phone_list));
    
    transition_indicator = time > prev_transition & time <= first_transition;
    sustain_indicator = time > first_transition & time <= final_transition;
    
    if sum(transition_indicator) > 0
        transition_path = ((time(transition_indicator) - prev_transition)'*this_feature + (first_transition - time(transition_indicator))'*prev_feature)/(first_transition - prev_transition);
        feature_vec(transition_indicator, :) = transition_path;
    end
    
    feature_vec(sustain_indicator, :) = repmat(this_feature, sum(sustain_indicator), 1);
    
    prev_transition = final_transition;
    prev_feature = this_feature;
end

sentence.feature_vec = feature_vec;

make_fig = 1;

if make_fig == 1
    
    figure
    
    imagesc(time, 1:size(feature_vec, 2), feature_vec')
    
    hold on
    
    plot(repmat(phone_transition_times, 2, 1), repmat([0; size(feature_vec, 1)], 1, length(phone_transition_times)), 'w', 'LineWidth', .5)
    
    xticks(phone_transition_times(1:(end - 1))+diff(phone_transition_times)/2)
    xticklabels(this_phone_list)
    xtickangle(45)
    
    xlim([min(phone_transition_times), max(phone_transition_times)])
    
end

%% Getting syllable boundaries.

syl_boundaries = round(100*syl_boundaries/16);

syllable_indicator = zeros(size(cochlear, 1), 1);

syllable_indicator(onset*100 + syl_boundaries) = 1;

syllable_indicator = conv(syllable_indicator, ones(1,500)/500, 'same');

sentence.syllables = syllable_indicator;

%% Getting wav file & applying cochlear filterbank.

this_wavfile = audioread([sentence_dir, wavfile_name]);
wavfile_length = 100*length(this_wavfile)/16;

loadload, paras(1) = 1;% paras(4) = 0;
this_auditory = wav2aud(this_wavfile, paras);

CF = 440 * 2 .^ ((-30:97)/24 - 1);
CF_included = CF(CF >= lowfreq & CF <= highfreq);
%CF_limits = [min(CF_included), max(CF_included)];
this_auditory = this_auditory(:, CF >= lowfreq & CF <= highfreq);

[~, aud_cols] = size(this_auditory);

sentence_length = ceil(max(length(time), wavfile_length + onset*100));

cochlear = zeros(sentence_length, no_channels);

channels_averaged = floor(aud_cols/no_channels);
for c = 1:no_channels
    %this_channel_freqs(c, :) = CF_included([(c - 1)*channels_averaged + 1, c*channels_averaged]);
    this_channel = nanmean(this_auditory(:,((c - 1)*channels_averaged + 1):c*channels_averaged), 2);
    this_channel = interp(interp(this_channel, 5), 10);
    if norm == 1
        this_channel = this_channel/sum(this_channel);
    elseif norm == 2
        this_channel = this_channel/max(this_channel);
    end
    cochlear(onset*100 + (1:length(this_channel)), c) = conv(this_channel, ones(1,500)/500, 'same');  
end

sentence.cochlear = cochlear;

end


function [format]=make_format(no_entries,entry_format)

format=[];
for i=1:no_entries
    format=[format,'%',entry_format,'\t'];
end
format=[format(1:end-1),'n'];

end

