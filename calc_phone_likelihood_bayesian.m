function result = calc_phone_likelihood_bayesian(window_length, sentence)

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
        
        this_result = calc_phone_likelihood_bayesian(window_length, this_struct);
        
        likelihood(w, :) = this_result.likelihood;
        
        time(w) = this_result.time;
        
    end
    
    result.likelihood = likelihood;
    result.time = time;
    
else
    
    likelihood = nan(1, num_phones);
    
    for p = 1:num_phones
        
        mu_p = phone_features(p, :);
        
        sigma_p = 2*cov(zscore(double(phone_features)))/3;
        
        likelihood(p) = prod(mvnpdf(feature_vec, mu_p, sigma_p));
        
%         for q = 1:num_phones
%             
%             mu_pq = 
%            
%             likelihood((p + 1)*num_phones + q) = 
%             
%         end
        
    end
    
    result.likelihood = nanunitsum(likelihood, 2);
    
    result.time = nanmean(time_vec);
    
end