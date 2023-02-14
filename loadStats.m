function stats = loadStats

%% Getting phoneme probabilities.
phone_data = load('phoneData.mat');
stats = phone_data;

%% Getting phone to feature mapping.
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

zero_trans_mat = cellfun(@(x) sum(abs(x), 2) < eps*100, transition_vectors);

trans_vector_column = cell2mat(transition_vectors(:));
zero_trans = sum(abs(trans_vector_column), 2) < eps*100;
trans_vector_column(zero_trans, :) = [];
trans_list_column = transition_list(:);
trans_list_column(zero_trans) = [];

%% Getting transition probabilities.
transition_data = load('phoneTransitions.mat');
trans_prob = transition_data.prob;
trans_prob(zero_trans_mat) = 0;

trans_count = transition_data.count;
trans_count(zero_trans_mat) = 0;
all_transitions = zscore(trans_vector_column); % zscore(diag(t_count(:))*trans_vector_column);
sigma_t = 2*cov(all_transitions)/3; % (all_transitions'*all_transitions)/(num_trans - 1);

stats.trans_phones = transition_data.timit_phonemes;
stats.trans_prob = trans_prob;
stats.trans_count = trans_count;
stats.sigma_t = sigma_t;

end