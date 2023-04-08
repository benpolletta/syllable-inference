function results = calc_phone_sequence_posterior_alt(likelihood, phone_sequence, transition_sequence, tsylb_option, method)

if nargin < 4, tsylb_option = []; end
if nargin < 5, method = []; end
if isempty(tsylb_option), tsylb_option = 1; end
if isempty(method), method = 'linear'; end

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
cdf_sequence = stats.cdf(prob_indices);
[~, expected_duration_indices] = cellfun(@(x) min(abs(x(:, 2) - .5)), cdf_sequence, 'UniformOutput', false);
expected_durations = cellfun(@(x, y) x(y), cdf_sequence, expected_duration_indices);

switch method

    case 'linear'

        total_duration = range(transition_sequence);
        expected_durations = nanunitsum(expected_durations)*total_duration;
        expected_transitions = cumsum([transition_sequence(1); expected_durations]);

        evidence_intervals = arrayfun(@(i) time >= expected_transitions(i) & time < expected_transitions(i + 1), (1:num_phones)', 'UniformOutput', false);
        evidence_intervals = cell2mat(evidence_intervals);
        phoneme_evidence = arrayfun(@(i) mean(likelihood.phone_likelihood(prob_indices(i), evidence_intervals(i, :))), 1:num_phones);

        duration_prob = arrayfun(@(i) prob_duration(phone_sequence(i), expected_durations(i)), 1:num_phones);

    case 'viterbi'

        expected_transitions = cumsum([transition_sequence(1), expected_durations]);

        TElength = length(expected_transitions);
        TOlength = length(transition_sequence);

        [P, D, T] = deal(zeros(TElength + 1, TOlength + 1));

        P(1, :) = 0:TElength;
        P(:, 1) = 0:TOlength;

        for i = 2:(TElength + 1)

            insertTEi_indicator = time >= expected_transitions(i - 1) & time <= expected_transitions(i - 1) + expected_durations(i);
            P(1, i) = P(1, i - 1) - sum(log(likelihood(prob_indices(i), insertTEi_indicator))) - log(.5);
            D(1, i) = expected_durations(i);
            T(1, i) = expected_transitions(i);

            for j = 2:(TOlength + 1)

                % Assuming next observed transition (jth) is same (ith) phoneme.
                absorbTOj_indicator = time >= T(i, j - 1) & time <= T(i, j - 1) + transition_sequence(i);
                absortbTOj = P(i, j - 1) - sum(log(likelihood(prob_indices(i), insertTEi_indicator))) - log(.5);

                % Assuming jth observed transition is spurious.
                insertTOj_indicator = time >= transitions(i - 1) & time <= transitions(i - 1) + sum(duration_sequence(j:(j + 1)));
                [~, prob_index] = min(abs(cdf_sequence{i}(:, 1) - sum(duration_sequence(j:(j + 1)))));
                duration_prob = cdf_sequence{i}(prob_index, 2);
                duration_prob = min(1 - duration_prob, duration_prob);
                insertTOj = P(i, j - 1) - sum(log(likelihood(prob_indices(i), phone_time_indicator))) - log(duration_prob);

                % Assumingr next observed transition is the right duration.
                insertTOj_indicator = time >= transitions(i - 1) & time <= transitions(i - 1) + duration_sequence(j);
                [~, prob_index] = min(abs(cdf_sequence{i}(:, 1) - sum(duration_sequence(j))));
                duration_prob = cdf_sequence{i}(prob_index, 2);
                duration_prob = min(1 - duration_prob, duration_prob);
                moveTEitoTOj = P(i - 1, j - 1) - sum(log(likelihood(prob_indices(i), phone_time_indicator))) - log(duration_prob);

                [P(i, j), index] = min([insertTEi, insertTOj, moveTEitoTOj]);

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

end

prior = calc_phone_sequence_prior(phone_sequence);

posterior = prior*prod(phoneme_evidence)*prod(duration_prob);

results = struct('prior', prior, 'posterior', posterior, 'phoneme_evidence', phoneme_evidence, 'duration_prob', duration_prob)
    
end

function prob = prob_duration(phone, duration) % (dist, duration)

global stats

cdf = stats.cdf{strcmpi(stats.id, phone)};

[~, cdf_index] = min(abs(cdf(:, 1) - duration));
prob = cdf(cdf_index, 2); 
prob = min(prob, 1 - prob);
% duration_hazard = cellfun(@(x, y) x(y, 2), stats.cdf, cdf_index);

% hazard = sum(duration_hazard.*nanunitsum(dist));

end