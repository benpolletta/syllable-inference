function [phone_list, feature_mat, unique_phones, features] = getWord2PhoneMap(transition_option)

if nargin == 0, transition_option = []; end
if isempty(transition_option), transition_option = 1; end

panphon_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/panphon/panphon/data/';
addpath(genpath(panphon_dir))

%% Getting phone to feature mapping.

feature_names = {'vowel', 'open', 'front', 'round', 'dipthong', 'rhoticized', 'place', 'manner', 'voice', 'affricate', 'syllabic'};
feature_dim = length(feature_names);

format = '%s';
for i = 1:feature_dim, format = [format, '%f']; end

phone2feature_fid = fopen([panphon_dir, 'timit_to_mfa.txt']); %'timit_features.txt']);
phone2feature_data = textscan(phone2feature_fid, format);
fclose(phone2feature_fid)

phone_list = phone2feature_data{1};
num_phones = length(phone_list);
feature_mat = [phone2feature_data{2:end}];

unique_phones = unique(phone_list);
    
for u = 1:length(unique_phones)
    features{u} = feature_mat(strcmpi(phone_list, unique_phones{u}), :);
end

if transition_option
    
    [phone_transition_list, phone_transition_vectors] = deal(cell(num_phones));
    
    for p = 1:num_phones
        
        for q = 1:num_phones
            
            phone_transition_list{p, q} = sprintf('%s2%s', phone_list{[p, q]});
            
            phone_transition_vectors{p, q} = feature_mat(q, :) - feature_mat(p, :);
            
        end
        
    end
    
end

phone_transition_list = phone_transition_list(:);
phone_transition_vectors = cell2mat(phone_transition_vectors(:));

save('phones2features.mat', 'feature_names', 'phone_list', 'feature_mat', 'unique_phones', 'features', 'phone_transition_list', 'phone_transition_vectors')