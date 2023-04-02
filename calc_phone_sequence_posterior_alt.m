function result = calc_phone_sequence_posterior_alt(likelihood, phone_sequence, transition_sequence, tsylb_option)

global stats

stats = loadStats(tsylb_option);

fields = fieldnames(likelihood);

for f = 1:length(fields)
    eval(sprintf('%s = likelihood.%s;', fields{f}, fields{f}))
end

num_phones = length(phone_sequence);

[lh_phones, prob_durations] = deal(ones(size(phone_sequence)));
duration_sequence = diff(transition_sequence);

num_phones = length(phone_sequence);

prob_indices = cellfun(@(x) find(strcmpi(stats.id, x), 1), phone_sequence);
trans_indices = cellfun(@(x) find(strcmpi(stats.trans_phones, x), 1), phone_sequence);

phones_prior = calc_phone_sequence_prior(phone_sequence);
cdf_sequence = stats.cdfs(prob_indices);
[~, expected_duration_indices] = cellfun(@(x) min(abs(x(:, 2) - .5)), cdf_sequence, 'UniformOutput', false);
expected_durations = cellfun(@(x, y) x(y), cdf_sequence, expected_duration_indices);

expected_transitions = cumsum([transition_sequence(1), expected_durations]);
TElength = length(expected_transitions);
TOlength = length(transition_sequence);

G = zeros(TElength + 1, TOlength + 1);

G(1, :) = 0:TElength;
G(:, 1) = 0:TOlength;

for i = 2:(TElength + 1)
    
    for j = 2:(TOlength + 1)

        % Inserting missed transition after ith expected duration.
        insertTEi_indicator = time >= transitions(i - 1) & time <= transitions(i - 1) + expected_durations(i);
        insertTEi = G(i - 1, j) - sum(log(likelihood(prob_indices(i), phone_time_indicator))) - log(.5);

        % Assuming jth observed transition is spurious.
        insertTOj_indicator = time >= transitions(i - 1) & time <= transitions(i - 1) + sum(duration_sequence(j:(j + 1)));
        [~, prob_index] = min(abs(cdf_sequence{i}(:, 1) - sum(duration_sequence(j:(j + 1)))));
        duration_prob = cdf_sequence{i}(prob_index, 2);
        duration_prob = min(1 - duration_prob, duration_prob);
        insertTOj = G(i, j - 1) - sum(log(likelihood(prob_indices(i), phone_time_indicator))) - log(duration_prob);
        
        % Assuming next observed transition is the right duration.
        insertTOj_indicator = time >= transitions(i - 1) & time <= transitions(i - 1) + duration_sequence(j);
        [~, prob_index] = min(abs(cdf_sequence{i}(:, 1) - sum(duration_sequence(j))));
        duration_prob = cdf_sequence{i}(prob_index, 2);
        duration_prob = min(1 - duration_prob, duration_prob);
        moveTEitoTOj = G(i - 1, j - 1) - sum(log(likelihood(prob_indices(i), phone_time_indicator))) - log(duration_prob);

        [G(i, j), index] = min([insertTEi, insertTOj, moveTEitoTOj]);
        
        switch index
            
            case 1
                
                durations(i) = expected_durations(i);
                duration_sequence((j + 2):(end + 1)) = duration_sequence((j + 1):end);
                duration_sequence(j) = duration_sequence(j) - expected_durations(i);
                
            case 2
                
                durations(i) = sum(duration_sequence(j:(j + 1)));
                duration_sequence((j + 1):(end - 1)) = duration_sequence((j + 2):end);
                
            case 3
                
                durations(i) = duration_sequence(j);
            
        end

        transitions(i) = transitions(i - 1) + durations(i);
        
    end
    
end

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


function hazard = phone_hazard(phone, duration) % (dist, duration)

global stats

[~, cdf_index] = cellfun(@(x) min(abs(x(:, 1) - duration)), stats.cdf, 'unif', 0);
hazard = cellfun(@(x, y) x(y, 2), stats.cdf, cdf_index); 
% duration_hazard = cellfun(@(x, y) x(y, 2), stats.cdf, cdf_index);

% hazard = sum(duration_hazard.*nanunitsum(dist));

end