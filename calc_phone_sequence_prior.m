function prior = calc_phone_sequence_prior(phone_sequence)

phone_data = load('normPhoneData.mat');
phones = phone_data.id;
phone_prob = phone_data.prob;
cdf = phone_data.cdf;


trans_data = load('phoneTransitions.mat');
trans_prob = trans_data.prob;
trans_phone_id = trans_data.phonemes;

prob_phones = ones(size(phone_sequence));

num_phones = length(phone_sequence);

timit_indices = cellfun(@(x) find(strcmpi(trans_phone_id, x), 1), phone_sequence);

prob_phones(1) = phone_prob(timit_indices(1));

for p = 2:num_phones

    prob_phones(p) = trans_prob(timit_indices(p), timit_indices(p - 1));

end

prior = prod(prob_phones);
    
end