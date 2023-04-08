function candidates = get_candidates_from_likelihood(likelihood)

if nargin < 2, window_length = []; end
if nargin < 3, method = []; end
if isempty(window_length), window_length = 1000; end % 2501; end
if isempty(method), method = 'ttest'; end

% phone_data = load('normPhoneData.mat');
% phones = phone_data.id;
% phone_prob = phone_data.prob;

phone2feature_data = load('phones2features.mat');
phone_list = phone2feature_data.phone_list;
phone_features = phone2feature_data.feature_mat;

result = struct('phones', {phone_list}, 'phone_features', phone_features);

feature_vec = sentence.input_vec;
time_vec = sentence.time;
num_features = size(feature_vec, 2);
num_phones = length(phone_list);
total_obs = size(feature_vec, 1);
        
sigma_p = 2*cov(zscore(double(phone_features)))/3;

num_windows = floor(total_obs/window_length);

mu_f = nan(num_features, num_windows); % (num_phones + 1)*num_phones);
sigma_f = nan(num_features, num_features, num_windows);
phone_stats = nan(num_phones, num_windows);
[trans_stats, inv_entropy] = deal(nan(num_windows, 1));

for w = 1:num_windows

    if trans_stats(w)

end