function [sylb_posterior, vocalic_nuclei, top_candidates] = calc_syllable_posterior_parallel(sentence, vowel_likelihood, likelihood, method, tsylb_option, cutoff)

if nargin < 4, method = ''; end
if nargin < 5, tsylb_option = []; end
if nargin < 6, cutoff = []; end
if isempty(method), method = 'syllable-specific'; end
if isempty(tsylb_option), tsylb_option = 1; end
if isempty(cutoff), cutoff = 0; end

global stats

stats = loadStats(tsylb_option);
sylbs = stats.sylbs.id;
probs = stats.sylbs.prob;
trans_matrix = stats.sylb_trans.prob;
trans_matrix = trans_matrix*diag(1./probs);

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

[onset, vocalic_nuclei] = deal([]);
sylb_prior = nanunitsum(ones(size(sylbs)));
hazards = nan(length(sylbs), length(likelihood.time));
[haz_exp, haz_ent] = deal(nan(size(likelihood.time)));
vn_i = 0;
% [last_transition_index, last_posterior_index] = deal(1);

for w = 1:length(time)

    %% Determining onset.
    if trans_likelihood(w) > .9 && isempty(onset)

        onset = time(w);
        continue

    end

    %% Defining duration.
    if vn_i == 0
        duration = 0;      
    elseif vn_i == 1
        duration = time(w) - onset;
    else
        duration = time(w) - vocalic_nuclei(vn_i);
    end

    %% Determining presence of vocalic nucleus.
    bu_nucleus_flag = nucleus_onset(w) == 1 && nucleus_onset(w - 1) == 0;

    [hazards(:, w), haz_exp(w), haz_ent(w)] = calc_hazard(sylb_prior(:, vn_i + 1), duration, method);
    td_nucleus_flag = haz_exp(w) > .6;

    if bu_nucleus_flag || td_nucleus_flag && vn_i > 0

        %% Recording new vocalic nucleus.
        vn_i = vn_i + 1;
        [vocalic_nuclei(vn_i), sylb_end] = deal(time(w));

        if vn_i == 1
            sylb_start = onset;
        else
            sylb_start = vocalic_nuclei(vn_i - 1);
        end

        sylb_index = time >= sylb_start & time <= sylb_end;
        trans_indicator = diff([0; trans_likelihood(sylb_index) > .9]) == 1;
        trans_indicator([1, end]) = 1;
        trans_time = time(sylb_index == 1);
        trans_times = trans_time(trans_indicator);
        this_sylb_likelihood = truncate_likelihood(likelihood, sylb_index);

        %% Getting candidate sequences of syllables.
        if vn_i == 1
            prev_candidate_sylbs = {};
        else
            prev_candidate_sylbs.sylbs = sylbs(sylb_posterior(:, vn_i - 1) > cutoff);
        end

        chunk_length = min(vn_i, 2);
        % length_function = @(x) get_vowel_onsets(x, vowels);
        tic
        chunk_structs = generate_sylb_sequences(prev_candidate_sylbs, sylbs, probs, trans_matrix, chunk_length, cutoff); %, @(x) get_vowel_onsets(x, vowels));
        toc

        % Retrieving fields & truncating chunks appropriately.
        chunks = cellfun(@(x) x.sylbs, chunk_structs, 'unif', 0);
        chunk_phones = cellfun(@(x) x.phones, chunk_structs, 'unif', 0);
        vowel_indicators = cellfun(@(x) x.vowel_indicator, chunk_structs, 'unif', 0);
        vowel_onsets = cellfun(@(x) find(diff(x) == 1), vowel_indicators, 'unif', 0);
        if vn_i == 1
            phone_sequences = cellfun(@(x,y) x(1:max(y(1) - 1, 1)), chunk_phones, vowel_onsets, 'unif', 0);
        else
            phone_sequences = cellfun(@(x,y) x(max(y(1), 1):max(y(2) - 1, 1)), chunk_phones, vowel_onsets, 'unif', 0);
        end

        %% Finding unique initial phone sequences.
        %joined_sequences = cellfun(@(x) strjoin(x, '/'), phone_sequences, 'UniformOutput', 0);
        [phone_seq_index, unique_phone_sequences] = findgroups(cellfun(@(x) strjoin(x, '/'), phone_sequences, 'UniformOutput', 0));
        unique_phone_sequences = cellfun(@(x) strsplit(x, '/'), unique_phone_sequences, 'UniformOutput', 0);
        ups_likelihood = zeros(size(unique_phone_sequences));

        %% Calculating likelihood of unique initial phone sequences & translating to chunk-coordinates.
        tic
        parfor ps = 1:length(unique_phone_sequences)
            phone_sequence_likelihood = calc_phone_sequence_likelihood(this_sylb_likelihood, unique_phone_sequences{ps}, trans_times);
            ups_likelihood(ps) = phone_sequence_likelihood.likelihood;
        end
        toc
        chunk_likelihood = ups_likelihood(phone_seq_index);

        %% Calculating syllable likelihood from chunk likelihood.
        tic
        parfor s = 1:length(sylbs)
            for n = 1:chunk_length % This is a kludge to get likelihood over syllable sequences including intermediate syllables...
                % Finding chunks ending with this syllable.
                chunk_indicator = cellfun(@(x) strcmp(x{n}, sylbs{s}), chunks);
                % Adding up those chunks.
                sl_from_cl(s, n) = sum(chunk_likelihood(chunk_indicator));
            end
        end
        toc
        sylb_likelihood(:, (vn_i - chunk_length + 1):vn_i) = flip(sl_from_cl, 2);

        %% Calculating sylb posterior distribution.
        sylb_posterior(:, vn_i) = nanunitsum(sylb_likelihood(:, vn_i).*sylb_prior(:, vn_i));
        [sorted_prob, sort_order] = sort(sylb_posterior(:, vn_i), 'descend');
        top_candidates{vn_i} = sylbs(sort_order(1:10));

        %% Calculating next prior from this posterior.
        sylb_prior(:, vn_i + 1) = trans_matrix*sylb_posterior(:, vn_i);

        %% Calculating rate posterior distribution.

    end

end

save(sprintf('sylbPosterior_%s.mat', datetime('now', 'Format', 'yy-MM-dd_HH-mm-ss')),...
    'sylb_posterior', 'top_candidates', 'sylb_likelihood', 'sylb_prior', 'hazards',...
    'haz_exp', 'haz_ent')

end

% function [vowel_length, vowel_onsets, vowel_phones] = get_vowel_onsets(chunk, vowels)
% 
% % Turning syllables into sequences of phonemes.
% phones = cellfun(@(y) split(y, '/'), chunk, 'unif', 0);
% phones = cat(1, phones{:});
% empty_phones = cellfun(@isempty, phones);
% phones = phones(~empty_phones);
% 
% % Finding vocalic nuclei & getting length (in vocalic nuclei).
% 
% vowel_length = 0;
% 
% if length(phones) > 0
% 
%     vowel_indicators = strcmp(vowels, phones);
%     vowel_bounds = find(diff([0; vowel_indicators]) == 1);
% 
%     vowel_length = length(vowel_bounds);
% 
% end
% 
% end

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