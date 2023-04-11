function stats = loadStats(tsylb_option)

if nargin < 1, tsylb_option = []; end
if isempty(tsylb_option); tsylb_option = 0; end

if tsylb_option
    phone_prefix = 'tsylbPhone';
else
    phone_prefix = 'phone'; 
end

%% Getting phoneme probabilities.
phone_data = load([phone_prefix, 'Data.mat']);
stats.phones = phone_data;

%% Getting phone to feature mapping.
phone2feature_data = load([phone_prefix, 's2features.mat']);
feature_names = phone2feature_data.feature_names;
phone_list = phone2feature_data.phone_list;
phone_features = phone2feature_data.feature_mat;
transition_list = phone2feature_data.unique_transition_list;
transition_vectors = phone2feature_data.unique_transition_vectors;

% bad_nx_index = find(strcmpi(phone_list, 'nx'), 1, 'last');
% transition_vectors(bad_nx_index, :) = [];
% transition_vectors(:, bad_nx_index) = [];
% transition_list(bad_nx_index, :) = [];
% transition_list(:, bad_nx_index) = [];

zero_trans_mat = cellfun(@(x) all(sum(abs(x), 2) < eps*100), transition_vectors);

trans_vector_column = transition_vectors(:);
zero_trans = cellfun(@(x) all(sum(abs(x), 2) < eps*100), trans_vector_column);
trans_vector_column(zero_trans, :) = [];
trans_list_column = transition_list(:);
trans_list_column(zero_trans) = [];

%% Getting transition probabilities.
transition_data = load([phone_prefix, 'Transitions.mat']);
trans_prob = transition_data.prob;
trans_prob(zero_trans_mat) = 0;

trans_count = transition_data.count;
trans_count(zero_trans_mat) = 0;
all_transitions = zscore(cell2mat(trans_vector_column)); % zscore(diag(t_count(:))*trans_vector_column);
sigma_t = 2*cov(all_transitions, 'partialrows')/3; % (all_transitions'*all_transitions)/(num_trans - 1);

stats.phone_trans.phones = transition_data.phonemes;
stats.phone_trans.prob = trans_prob;
stats.phone_trans.count = trans_count;
stats.phone_trans.sigma = sigma_t;

%% Loading syllable & word data.

datafiles = {'normSylbData.mat', 'sylbTransitions.mat', 'normWordCanonData.mat', 'wordTransitions.mat'};
struct_names = {'sylbs', 'sylb_trans', 'words', 'word_trans'};

for d = 1:length(datafiles)
    data = load(datafiles{d});
    fields = fieldnames(data);

    for f = 1:length(fields)
        eval(sprintf('stats.%s.%s = data.%s;', struct_names{d}, fields{f}, fields{f}));
    end

end
