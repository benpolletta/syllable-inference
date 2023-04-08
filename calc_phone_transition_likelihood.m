function result = calc_phone_transition_likelihood(sentence, window_length, method)

if nargin < 2, window_length = []; end
if nargin < 3, method = []; end
if isempty(window_length), window_length = 1000; end
if isempty(method), method = 'ttest'; end

% phone_data = load('normPhoneData.mat');
% phones = phone_data.id;
% phone_prob = phone_data.prob;

phone2feature_data = load('phones2features.mat');
feature_names = phone2feature_data.feature_names;
phone_list = phone2feature_data.phone_list;
phone_features = phone2feature_data.feature_mat;
transition_list = phone2feature_data.phone_transition_list;
transition_vectors = phone2feature_data.phone_transition_vectors;

bad_nx_index = find(strcmpi(phone_list, 'nx'), 1, 'last');
transition_vectors(bad_nx_index, :) = [];
transition_vectors(:, bad_nx_index) = [];
transition_list(bad_nx_index, :) = [];
transition_list(:, bad_nx_index) = [];

trans_vector_column = cell2mat(transition_vectors(:));
zero_trans = sum(abs(trans_vector_column), 2) < eps*100;
trans_vector_column(zero_trans, :) = [];
trans_list_column = transition_list(:);
trans_list_column(zero_trans) = [];

transition_data = load('phoneTransitions.mat');
prob = transition_data.prob;

t_count = transition_data.count;
all_transitions = zscore(trans_vector_column); % zscore(diag(t_count(:))*trans_vector_column);
sigma_t = 2*cov(all_transitions)/3; % (all_transitions'*all_transitions)/(num_trans - 1);

% figure, imagesc(sigma_t)
% set(gca, 'XTick', 1:length(feature_names), 'XTickLabel', feature_names)
% set(gca, 'YTick', 1:length(feature_names), 'YTickLabel', feature_names)

result = struct('phones', {phone_list}, 'phone_features', {phone_features},...
    'transition_vectors', {transition_vectors}, 'transition_list', {trans_list_column}, 'sigma_t', sigma_t, 'prob', prob);

feature_vec = diff(sentence.input_vec);
time_vec = sentence.time(1:(end - 1));
num_features = size(feature_vec, 2);
num_trans = size(trans_vector_column, 1);
num_obs = size(feature_vec, 1);

if num_obs > window_length
    
    num_windows = floor(num_obs/window_length);
    
    likelihood = ones(num_windows, num_trans); % (num_phones + 1)*num_phones);
    
    for w = 1:num_windows
        
        window_indices = ((w - 1)*window_length + 1):w*window_length;
        
        this_feature_vec = feature_vec(window_indices, :);
        
        this_struct = struct('input_vec', this_feature_vec, 'time', time_vec(window_indices));
        
        this_result = calc_phone_transition_likelihood(this_struct, window_length);
        
        likelihood(w, :) = this_result.likelihood;
        
        time(w) = this_result.time;
        
    end
    
    result.likelihood = likelihood;
    
    result.time = time;
    
else
            
    stats = nan(1, num_trans);
    
    switch method
        
        case 'ttest'
            
            mu_t_pages = repmat(permute(trans_vector_column, [3, 2, 1]), [num_obs, 1, 1]);
            feature_pages = repmat(feature_vec, [1, 1, num_trans]);
            feature_mean_diff = feature_pages - mu_t_pages;
            
            sigma_f = pagemtimes(feature_mean_diff, 'transpose', feature_mean_diff, 'none')/(num_obs - 1);
            
            mu_diff_pages = repmat(sum(feature_vec)/num_obs, [1, 1, num_trans]) - permute(trans_vector_column, [3, 2, 1]);
            
            tstat = squeeze(pagemtimes(mu_diff_pages, pagemldivide(sigma_f, pagetranspose(mu_diff_pages))));
            
            stats = tcdf(tstat, num_obs - 1, 'upper');
            
        case 'meanpdf'
            
            mu_f = sum(feature_vec)/num_obs;
            
            stats = mvnpdf(mu_f, trans_vector_column, sigma_t);
            
        case 'prodpdf'
            
            stats = prod(mvnpdf(feature_vec, trans_vector_column, sigma_t));
        
    end
    
    result.likelihood = stats; % result.likelihood = nanunitsum(stats, 2);
    
    result.time = nanmean(time_vec);
    
end