function result = calc_phone_sequence_posterior_alt(likelihood, phone_sequence, transition_sequence, tsylb_option)

global stats

stats = loadStats(tsylb_option);

fields = fieldnames(likelihood);

for f = 1:length(fields)
    eval(sprintf('%s = likelihood.%s;', fields{f}, fields{f}))
end

num_phones = length(phone_sequence);

phones_prior = calc_phone_sequence_prior(phone_sequence);

[lh_phones, prob_durations] = deal(ones(size(phone_sequence)));
duration_sequence = diff(transition_sequence);

num_phones = length(phone_sequence);

prob_indices = cellfun(@(x) find(strcmpi(stats.id, x), 1), phone_sequence);
trans_indices = cellfun(@(x) find(strcmpi(stats.trans_phones, x), 1), phone_sequence);

for p = 1:num_phones

    phone_time_indicator = time >= transition_sequence(p) & time <= transition_sequence(p + 1);
    lh_phones(p) = sum(likelihood(prob_indices(p), phone_time_indicator));

    prior_phones(p) = trans_prob(trans_indices(p), trans_indices(p - 1));

    this_cdf = stats.cdf{prob_indices(p)};

    [~, cdf_index] = min(abs(this_cdf(:, 1) - duration_sequence(p)));

    prob_durations(p - 1) = this_cdf(cdf_index, 2);
    prob_durations(p - 1) = min(1 - prob_durations, prob_durations);

end

prior = phones_prior*prod(lh_phones)*prod(prob_durations);
    
end