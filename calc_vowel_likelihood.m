function result = calc_vowel_likelihood(sentence, window_length)

if nargin < 2, window_length = []; end
if nargin < 3, method = []; end
if isempty(window_length), window_length = 10000; end % 2501; end

phone2feature_data = load('phones2features.mat');
phone_list = phone2feature_data.phone_list;
phone_features = phone2feature_data.feature_mat;
feature_names = phone2feature_data.feature_names;

result = struct('phones', {phone_list}, 'phone_features', phone_features, 'feature_names', {feature_names});

vowel_feature_index = find(strcmp(feature_names, 'vowel'));

feature_vec = sentence.input_vec;
vowel_vec = feature_vec(:, vowel_feature_index);
time_vec = sentence.time;
num_phones = length(phone_list);
total_obs = length(vowel_vec);
        
sigma_p = 2*cov(zscore(double(phone_features)))/3;

window_start_obs = total_obs - window_length;
window_step_length = window_length/10;
num_windows = floor(window_start_obs/window_step_length);

[mu_v, sigma_v, v_stats] = deal(nan(num_windows, 1)); % (num_phones + 1)*num_phones);

for w = 1:num_windows

    window_start = (w - 1)*window_step_length + 1;
    window_end = (w - 1)*window_step_length + window_length;
    window_indices = window_start:window_end;

    time(w) = mean(time_vec(window_indices));

    this_vowel_vec = vowel_vec(window_indices);
    
    mu_v(w) = nanmean(this_vowel_vec);

    sigma_v(w) = nanmean((this_vowel_vec - mu_v(w)).^2);

    v_tstat = mu_v(w)/sqrt(sigma_v(w));
    v_stats(w) = 1 - tcdf(v_tstat, window_length - 1, 'upper');

end

[v_peak, v_avg, v_std] = ThresholdingAlgo(v_stats, 5, 3, 0.25);
[d2v_peak, d2v_avg, d2v_std] = ThresholdingAlgo(diff(diff([0; v_stats; 0])), 5, 3, 0.25);
v_indicator = d2v_peak == -1 & v_peak == 1; % v_peak(2:end - 1) == 1;
v_indicator = diff(v_indicator) == 1;
v_indicator(end + 1) = 0;

% indicator_time = time(1:end - 1) + diff(time)/2;
% indicator_time = indicator_time(1:end - 1) + diff(indicator_time)/2;

result.time = time;
result.vowel_vec = vowel_vec;
result.vowel_likelihood = v_stats;
% result.indicator_time = indicator_time;
result.v_indicator = v_indicator;
result.v_avg_filter = v_avg;
result.v_std_filter = v_std;
result.mu_v = mu_v;
result.sigma_v = sigma_v;
    
end