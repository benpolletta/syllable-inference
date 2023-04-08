function [pvalues, phone_dist, times] = calc_phone_likelihood_ttest(window_length, feature_vec, phones, phone_features, time_vec)

num_features = size(feature_vec, 2);
num_phones = length(phones);
num_obs = size(feature_vec, 1);

if num_obs > window_length
    
    num_windows = floor(num_obs/window_length);
    
    tstats = ones(num_windows, num_phones);
    
    for w = 1:num_windows
        
        window_indices = ((w - 1)*window_length + 1):w*window_length;
        
        this_feature_vec = feature_vec(window_indices, :);
        
        [this_pvalues, this_phone_dist] = calc_phone_likelihood_ttest(window_length, this_feature_vec, phones, phone_features);
        
        pvalues(w, :) = this_pvalues;
        phone_dist(w, :) = this_phone_dist;
        
        this_time = time_vec(window_indices);
        
        times(w) = mean(this_time);
        
    end
    
else
    
    tstats = nan(1, num_phones);
    
    for p = 1:num_phones
        
        mu_f = sum(feature_vec)/num_obs;
        
        mu_p = phone_features(p, :);
        
        feature_mean_diff = feature_vec - repmat(mu_p, num_obs, 1);
        
        sigma_f = ((feature_mean_diff')*feature_mean_diff)/(num_obs - 1);
        
        tstats(p) = (mu_f - mu_p)*inv(sigma_f)*(mu_f - mu_p)';
        
    end
    
    pvalues = tcdf(tstats, num_obs - 1, 'upper');
    
    phone_dist = nanunitsum(pvalues, 2);
    
end