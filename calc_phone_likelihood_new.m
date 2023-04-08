function result = calc_phone_likelihood_new(sentence, window_length, method, tsylb_option)

if nargin < 2, window_length = []; end
if nargin < 3, method = []; end
if nargin < 4, tsylb_option = []; end
if isempty(window_length), window_length = 1000; end % 2501; end
if isempty(method), method = 'ttest'; end
if isempty(tsylb_option), tsylb_option = 1; end

% phone_data = load('normPhoneData.mat');
% phones = phone_data.id;
% phone_prob = phone_data.prob;
if tsylb_option
    phone2feature_data = load('tsylbPhones2features.mat');
else
    phone2feature_data = load('phones2features.mat');
end

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

[mu_f, mu_diff, mu_delta] = deal(nan(num_features, num_windows)); % (num_phones + 1)*num_phones);
sigma_f = nan(num_features, num_features, num_windows);
phone_stats = nan(num_phones, num_windows);
[inv_entropy, trans_stats, delta_stats] = deal(nan(num_windows, 1));

for w = 1:num_windows

    window_indices = ((w - 1)*window_length + 1):w*window_length;

    time(w) = mean(time_vec(window_indices));

    this_feature_vec = feature_vec(window_indices, :);
    
    mu_p_pages = repmat(permute(phone_features, [3, 2, 1]), [window_length, 1, 1]);
    feature_pages = repmat(this_feature_vec, [1, 1, num_phones]);
    feature_mean_diff = feature_pages - mu_p_pages;

    sigma_f = pagemtimes(feature_mean_diff, 'transpose', feature_mean_diff, 'none')/(window_length - 1);

    mu_diff_pages = repmat(mean(this_feature_vec), [1, 1, num_phones]) - permute(phone_features, [3, 2, 1]);

    phone_tstat = squeeze(pagemtimes(mu_diff_pages, pagemldivide(sigma_f, pagetranspose(mu_diff_pages))));
    phone_stats(:, w) = tcdf(phone_tstat, window_length - 1, 'upper');

    norm_phone_stats = nanunitsum(phone_stats(:, w));
    inv_entropy(w) = (log(num_phones) + nansum(norm_phone_stats.*log(norm_phone_stats)))/num_phones;

    mu_f(:, w) = mean(this_feature_vec);
%     mu_f_column = repmat(mu_f(:, w)', num_phones, 1);
% 
    sigma_f(:, :, w) = cov(this_feature_vec);
% 
%     %phone_tstat = (mu_f_column - phone_features)*(sigma_f(:, :, w)\(mu_f_column - phone_features)');
%     phone_tstat = (mu_f_column - phone_features)*(sigma_p\(mu_f_column - phone_features)');


    if w > 1

        feature_trans = this_feature_vec - repmat(mu_f(:, w - 1)', window_length, 1);
        sigma_t = cov(feature_trans);
        mu_diff(:, w) = diff(mu_f(:, (w - 1):w), [], 2);

        trans_tstat = (mu_diff(:, w)')*(sigma_t\mu_diff(:, w));

        trans_stats(w) = tcdf(trans_tstat, window_length - 1);

        feature_delta = diff(this_feature_vec);
        mu_delta(:, w) = mean(feature_delta);
        sigma_delta = cov(feature_delta);

        delta_tstat = (mu_delta(:, w)')*(sigma_delta\mu_delta(:, w));

        delta_stats(w) = tcdf(delta_tstat, window_length - 1);

    end

end

result.time = time;
result.phone_likelihood = phone_stats;
result.trans_likelihood = trans_stats;
result.delta_likelihood = delta_stats;
result.inv_entropy = inv_entropy;
    
end