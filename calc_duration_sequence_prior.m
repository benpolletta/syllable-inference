function result = calc_duration_sequence_prior(phone_sequence, duration_sequence)

phone_data = load('normPhoneData.mat');
phones = phone_data.id;
phone_prob = phone_data.prob;
cdf = phone_data.cdf;


trans_data = load('phoneTransitions.mat');
trans_prob = trans_data.prob;
trans_phone_id = trans_data.timit_phonemes;

[prob_phones, prob_durations] = deal(ones(size(phone_sequence)));

num_phones = length(phone_sequence);

timit_indices = cellfun(@(x) find(strcmpi(trans_phone_id, x), 1), phone_sequence);

prob_phones(1) = phone_prob(timit_indices(1));

for p = 2:num_phones

    prob_phones(p) = trans_prob(timit_indices(p), timit_indices(p - 1));

    this_cdf = cdf{timit_indices(p)};

    [~, cdf_index] = min(abs(this_cdf(:, 1) - duration_sequence(p));

    prob_durations(p - 1) = this_cdf(cdf_index, 2);
    prob_durations(p - 1) = min(1 - prob_durations, prob_durations);

end

prior = prod(prob_phones)*prod(prob_durations);
    
end