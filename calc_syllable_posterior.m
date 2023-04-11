function [sylb_posterior, vocalic_nuclei] = calc_syllable_posterior(sentence, vowel_likelihood, likelihood, method, tsylb_option, cutoff)

if nargin < 4, method = ''; end
if nargin < 5, tsylb_option = []; end
if nargin < 6, cutoff = []; end
if isempty(method), method = 'syllable-specific'; end
if isempty(tsylb_option), tsylb_option = 1; end
if isempty(cutoff), cutoff = 0; end

global stats

stats = loadStats(tsylb_option);
sylbs = stats.sylbs.id;
trans_matrix = stats.sylb_trans.prob_cols;

[tsylb_phonemes, class_indicator, class_names] = getPhones(1);
vowels = tsylb_phonemes(class_indicator(:, strcmpi(class_names, 'vowels')));
vowels = {vowels{:}, 'el', 'em', 'en', 'enx'};

% fields = fieldnames(likelihood);
% 
% for f = 1:length(fields)
%     eval(sprintf('%s = likelihood.%s;', fields{f}, fields{f})) 
% end

time = likelihood.time;
phone_likelihood = likelihood.phone_likelihood;
trans_likelihood = likelihood.trans_likelihood;

nucleus_onset_times = vowel_likelihood.time(vowel_likelihood.v_indicator);
nucleus_onset = zeros(size(time));
for t = 1:length(nucleus_onset_times)
    [~, nuc_index] = min(abs(time - nucleus_onset_times(t)));
    nucleus_onset(nuc_index) = 1;
end

%% Get estimated syllable duration.

%% Calculate posterior for each time.

vocalic_nuclei = [];
[sylb_posterior, sylb_likelihood, sylb_prior] = deal(zeros(length(sylbs), 50));
vn_i = 1;
% [last_transition_index, last_posterior_index] = deal(1);

for w = 1:length(time)

    %% Determining onset.
    if trans_likelihood(w) > .9 && isempty(vocalic_nuclei)

        vocalic_nuclei(1) = time(w);
        sylb_prior(:, vn_i) = stats.sylbs.prob;
        continue

    end

    %% Defining duration.
    if ~isempty(vocalic_nuclei)

        duration = time(w) - vocalic_nuclei(end);

    else

        duration = 0;

    end

    %% Determining presence of vocalic nucleus.
    bu_nucleus_flag = nucleus_onset(w) == 1 && nucleus_onset(w - 1) == 0;

    [hazards, haz_exp, has_ent] = calc_hazard(sylb_prior(:, vn_i), duration, method);
    td_nucleus_flag = haz_exp > .6;

    if bu_nucleus_flag || td_nucleus_flag && ~isempty(vocalic_nuclei)

        %% Recording new vocalic nucleus.
        vn_i = vn_i + 1;
        vocalic_nuclei(vn_i) = time(w);

        sylb_end = vocalic_nuclei(vn_i);

        sylb_start = vocalic_nuclei(vn_i - 1);

        sylb_index = time >= sylb_start & time <= sylb_end;
        trans_sequence = [0; diff(trans_likelihood(sylb_index) > .9) == 1];
        trans_sequence([1, end]) = 1;
        this_sylb_likelihood = truncate_likelihood(likelihood, sylb_index);

        sylb_prior(:, vn_i) = trans_matrix*sylb_posterior(:, vn_i - 1);
        prev_candidate_sylbs = sylbs(sylb_posterior(:, vn_i - 1) > cutoff);

        % Getting sequences of syllables.
        chunks = generate_sequences(prev_candidate_sylbs, sylbs, trans_matrix, 2, cutoff);
        chunk_likelihood = nan(size(chunks));

        for c = 1:length(chunks)
            % Turning syllables into sequences of phonemes.
            chunk_phones = cellfun(@(x) split(x, '/'), chunks{c}, 'unif', 0);
            chunk_phones = cat(1, chunk_phones{:});
            chunk_phones(cellfun(@isempty, chunk_phones)) = [];

            % Finding vocalic nuclei & truncating appropriately.
            vowel_indicator = cellfun(@(x) any(strcmp(vowels, x)), chunk_phones);
            vowel_bounds = find(diff([0; vowel_indicator]) == 1);
            if vn_i == 2
                phone_sequence = chunk_phones(1:vowel_bounds(1));
            else
                phone_sequence = chunk_phones(vowel_bounds(1):(vowel_bounds(2) - 1));
            end

            % Calculating phone sequence likelihood.
            phone_sequence_likelihood = calc_phone_sequence_likelihood(this_sylb_likelihood, phone_sequence, trans_sequence);
            chunk_likelihood(c) = phone_sequence_likelihood.likelihood;
        end

        for s = 1:length(sylbs)
            % Finding chunks ending with this syllable.
            chunk_indicator = cellfun(@(x) strcmp(x{end}, sylbs{s}), chunks);
            % Adding up those chunks.
            sylb_likelihood(s, vn_i) = sum(chunk_likelihood(chunk_indicator));
        end

        %% Calculating sylb posterior distribution.
        sylb_posterior(:, vn_i) = sylb_likelihood(:, vn_i).*sylb_prior(:, vn_i);
        [~, i] = max(sylb_posterior, 'first', 10);
        top_candidates{vn_i} = sylbs(i)

        %% Calculating rate posterior distribution.

    end

end

end

function likelihood_out = truncate_likelihood(likelihood_in, index)

    index_length = length(index);

    these_fields = fieldnames(likelihood_in);

    for f = 1:length(these_fields)

        this_field = likelihood_in.(these_fields{f});

        [r, c] = size(this_field);

        if r == index_length

            this_field = this_field(index, :);

        elseif c == index_length

            this_field = this_field(:, index);

        end

        likelihood_out.(these_fields{f}) = this_field;

    end

end

function [hazards, haz_exp, haz_ent] = calc_hazard(dist, duration, method) % (dist, duration)

global stats

switch method

    case 'syllable-specific'

        [~, cdf_index] = cellfun(@(x) min(abs(x(:, 1) - duration)), stats.sylbs.cdf, 'unif', 0);
        hazards = cellfun(@(x, y) x(y, 2), stats.sylbs.cdf, cdf_index);
        haz_exp = sum(hazards.*nanunitsum(dist));
        haz_ent = -nansum(nanunitsum(hazards).*log(nanunitsum(hazards)))/log(length(dist)); % (log(num_phones) + nansum(hazard(:, w).*log(hazard(:, w))))/num_phones;

    case 'non-syllable-specific'

        [~, cdf_index] = min(abs(stats.sylbs.cdf(:, 1) - duration));
        haz_exp = stats.sylbs_cdf(cdf_index, 2);
        haz_ent = 0;
        hazards = haz_est*ones(size(this_dist));

    case 'zero_hazard'

        haz_exp(w) = 0;
        haz_ent(w) = 0;
        hazard(w) = zeros(size(dist));

end

end

function intensity = phone_intensity(dist, duration)

global stats

pdf = cellfun(@(x) [cumsum(x(:, 1) + diff([x(:, 1); x(end, 1)])/2), diff]);

[~, cdf_index] = cellfun(@(x) min(abs(x(:, 1) - duration)), stats.cdf, 'unif', 0);
duration_hazard = cellfun(@(x, y) x(y, 2), stats.cdf, cdf_index);

hazard = sum(duration_hazard.*nanunitsum(dist));

end