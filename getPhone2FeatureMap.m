function [phone_list, feature_mat, unique_phones, features] = getPhone2FeatureMap(transition_option, tsylb_flag)

if nargin == 0, transition_option = []; end
if isempty(transition_option), transition_option = 0; end
if nargin < 2, tsylb_flag = []; end
if isempty(tsylb_flag), tsylb_flag = 0; end

panphon_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/panphon/panphon/data/';
addpath(genpath(panphon_dir))

%% Getting phone to feature mapping.

feature_names = {'vowel', 'open', 'front', 'round', 'dipthong', 'rhoticized', 'place', 'manner', 'voice', 'affricate', 'syllabic'};
feature_dim = length(feature_names);

format = '%s';
for i = 1:feature_dim, format = [format, '%f']; end

if tsylb_flag
    save_name = 'tsylbPhones2features.mat';

    phone2feature_fid = fopen([panphon_dir, 'tsylb_to_mfa.txt']); %'timit_features.txt']);
    phone2feature_data = textscan(phone2feature_fid, format);
    fclose(phone2feature_fid)

    phone_list = phone2feature_data{2};
    num_phones = length(phone_list);
    feature_mat = [phone2feature_data{3:end}];
else
    save_name = 'phones2features.mat';

    phone2feature_fid = fopen([panphon_dir, 'timit_to_mfa.txt']);
    phone2feature_data = textscan(phone2feature_fid, format);
    fclose(phone2feature_fid)

    phone_list = phone2feature_data{1};
    num_phones = length(phone_list);
    feature_mat = [phone2feature_data{2:end}];
end

unique_phones = unique(phone_list);

if tsylb_flag, unique_phones(strcmpi(unique_phones, '-')) = []; end
    
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
    
    save(save_name, 'feature_names', 'phone_list', 'feature_mat', 'unique_phones', 'features', 'phone_transition_list', 'phone_transition_vectors')

else

    save(save_name, 'feature_names', 'phone_list', 'feature_mat', 'unique_phones', 'features')

end

% phone_transition_list = phone_transition_list(:);
% phone_transition_vectors = cell2mat(phone_transition_vectors(:));
