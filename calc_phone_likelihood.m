function result = calc_phone_likelihood(sentence, window_length, method)

if nargin < 2, window_length = []; end
if nargin < 3, method = []; end
if isempty(window_length), window_length = 2501; end
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
num_obs = size(feature_vec, 1);

if num_obs > window_length
    
    num_windows = floor(num_obs/window_length);
    
    likelihood = ones(num_windows, num_phones); % (num_phones + 1)*num_phones);
    
    for w = 1:num_windows
        
        window_indices = ((w - 1)*window_length + 1):w*window_length;
        
        this_feature_vec = feature_vec(window_indices, :);
        
        this_struct = struct('input_vec', this_feature_vec, 'time', time_vec(window_indices));
        
        this_result = calc_phone_likelihood(this_struct, window_length);
        
        likelihood(w, :) = this_result.likelihood;
        
        time(w) = this_result.time;
        
    end
    
    result.likelihood = likelihood;
    
    result.time = time;
    
else
            
    stats = nan(1, num_phones);
    
    for p = 1:num_phones
        
        mu_p = phone_features(p, :);
        
        sigma_p = 2*cov(zscore(double(phone_features)))/3;
        
        switch method
            
            case 'ttest'
                
                mu_f = sum(feature_vec)/num_obs;
                
                feature_mean_diff = feature_vec - repmat(mu_p, num_obs, 1);
                
                sigma_f = ((feature_mean_diff')*feature_mean_diff)/(num_obs - 1);
                
                tstat = (mu_f - mu_p)*(sigma_f\(mu_f - mu_p)');
                
                stats(p) = tcdf(tstat, num_obs - 1, 'upper');
                
            case 'meanpdf'
                
                mu_f = sum(feature_vec)/num_obs;
                
                stats(p) = mvnpdf(mu_f, mu_p, sigma_p);
                
            case 'prodpdf'
                
                stats(p) = prod(mvnpdf(feature_vec, mu_p, sigma_p));
                
        end
        
    end
    
    result.likelihood = stats; % result.likelihood = nanunitsum(stats, 2);
    
    result.time = nanmean(time_vec);
    
end