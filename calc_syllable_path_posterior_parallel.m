function [path_posterior, vocalic_nuclei, top_candidates] = calc_syllable_path_posterior_parallel(sentence, vowel_likelihood, likelihood, method, tsylb_option, cutoff, debug_flag)

if nargin < 4, method = ''; end
if nargin < 5, tsylb_option = []; end
if nargin < 6, cutoff = []; end
if nargin < 7, debug_flag = []; end
if isempty(method), method = 'non-syllable-specific'; end
if isempty(tsylb_option), tsylb_option = 1; end
if isempty(cutoff), cutoff = 0; end
if isempty(debug_flag), debug_flag = false; end

name_prefix = 'sylbPosterior';
if debug_flag, name_prefix = 'sylbNOposterior'; end
sim_name = sprintf('%s_%s_%s.mat', name_prefix, replace(sentence.filename, '/', '-'), datetime('now', 'Format', 'yy-MM-dd_HH-mm-ss'));

global stats vowels

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
    nucleus_onset(nuc_index(1)) = 1;
end

%% Get estimated syllable duration.

mean_dur = stats.sylb_dur.stats(1); %.mean;
var_dur = stats.sylb_dur.stats(2); %.var;

%% Calculate posterior for each time.

[onset, vocalic_nuclei] = deal([]);
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

    [hazards(:, w), haz_exp(w), haz_ent(w)] = calc_hazard(probs, duration, method);
    % td_nucleus_flag = haz_exp(w) > .6;

    if bu_nucleus_flag % || td_nucleus_flag && vn_i > 0

        %% Recording new vocalic nucleus.
        vn_i = vn_i + 1;
        [vocalic_nuclei(vn_i), sylb_end] = deal(time(w));

        if vn_i == 1
            sylb_start = onset;
        else
            sylb_start = vocalic_nuclei(vn_i - 1);
        end

        % Estimating number of syllables.
        abs_dur = sylb_end - sylb_start;
        est_sylb_dur = ceil(abs_dur/(mean_dur - var_dur));

        sylb_index = time >= sylb_start & time <= sylb_end;
        trans_indicator = diff([0; trans_likelihood(sylb_index) > .9]) == 1;
        trans_indicator([1, end]) = 1;
        trans_time = time(sylb_index == 1);
        trans_times = trans_time(trans_indicator);
        this_phone_likelihood = truncate_likelihood(likelihood, sylb_index);

        %% Getting candidate sequences of syllables.
        if vn_i == 1 % Initializing empty paths.
            prev_paths = generate_sylb_sequences({}, sylbs, probs, trans_matrix, vowels, 0);
            % prev_posterior = cellfun(@(x) x.prob, prev_paths);
        end
        % likely_paths = find(prev_posterior >= cutoff);

        min_path_vowel_onsets = vn_i;
        max_chunk_sylbs = max(3, est_sylb_dur);
        % length_function = @(x) get_vowel_onsets(x, vowels);
        [path_structs, path_priors] = deal(cell(size(prev_paths)));
        fprintf('\n Calculating continuations for %d paths.\n', length(prev_paths))
        tic
        parfor p = 1:length(prev_paths) % min(100, length(prev_paths))
            this_path = setfield(setfield(prev_paths{p}, 'recur', 0), 'prob', 1);
            path_structs{p} = generate_sylb_sequences(this_path, sylbs, probs, trans_matrix, vowels, min_path_vowel_onsets, cutoff, max_chunk_sylbs); %, @(x) get_vowel_onsets(x, vowels));
            path_priors{p} = prev_paths{p}.prob*ones(size(path_structs{p}));
        end
        toc
        path_structs = cat(1, path_structs{:});
        candidate_paths{vn_i} = path_structs;
        chunk_prior{vn_i} = cellfun(@(x) x.prob^(1/x.recur), path_structs);
        path_priors = chunk_prior{vn_i}.*cat(1, path_priors{:});

        %% Calculating likelihood of chunk.
        % Retrieving fields & truncating chunks appropriately.
        % path_sylbs = cellfun(@(x) x.sylbs, path_structs, 'unif', 0);
        path_phones = cellfun(@(x) x.phones, path_structs, 'unif', 0);
        vowel_indicators = cellfun(@(x) x.vowel_indicator, path_structs, 'unif', 0);
        vowel_onsets = cellfun(@(x) find(diff(x) == 1), vowel_indicators, 'unif', 0);
        if vn_i == 1
            chunk_phones = cellfun(@(x,y) x(1:max(y(1) - 1, 1)), path_phones, vowel_onsets, 'unif', 0);
        else
            chunk_phones = cellfun(@(x,y) x(max(y(vn_i - 1), 1):max(y(vn_i) - 1, 1)), path_phones, vowel_onsets, 'unif', 0);
        end

        % Finding unique initial phone sequences.
        %joined_sequences = cellfun(@(x) strjoin(x, '/'), phone_sequences, 'UniformOutput', 0);
        [phone_seq_index, unique_phone_sequences] = findgroups(cellfun(@(x) strjoin(x, '/'), chunk_phones, 'UniformOutput', 0));
        unique_phone_sequences{vn_i} = cellfun(@(x) strsplit(x, '/'), unique_phone_sequences, 'UniformOutput', 0);
        ups_likelihood = zeros(size(unique_phone_sequences{vn_i}));

        % Calculating likelihood of unique initial phone sequences & translating to chunk-coordinates.
        fprintf('\n Calculating likelihood for %d phone sequences.\n', length(ups_likelihood))
        tic
        if debug_flag
            ups_likelihood = nanunitsum(exprnd(5*ones(size(ups_likelihood))));
        else
            parfor ps = 1:length(unique_phone_sequences{vn_i})
                phone_sequence_likelihood = calc_phone_sequence_likelihood(this_phone_likelihood, unique_phone_sequences{vn_i}{ps}, trans_times);
                % ups_likelihood_struct(ps) = phone_sequence_likelihood;
                ups_likelihood(ps) = phone_sequence_likelihood.likelihood;
            end
        end
        toc
        chunk_likelihood{vn_i} = ups_likelihood(phone_seq_index);

        %% Calculating path posterior distribution.
        path_posterior = nanunitsum(chunk_likelihood{vn_i}.*path_priors);
        posterior_cell = mat2cell(path_posterior, ones(size(path_posterior, 1), 1), ones(size(path_posterior, 2), 1));
        path_structs = cellfun(@(x, y) setfield(x, 'prob', y), path_structs, posterior_cell, 'unif', 0);
        [sorted_prob, sort_order] = sort(path_posterior, 'descend');
        top_index = find(sorted_prob == sorted_prob(1), 1, 'last');
        cutoff_index = find(sorted_prob >= cutoff, 1, 'last');
        if isempty(cutoff_index), cutoff_index = 0; end
        if vn_i == 1
            keep_index = max([100, top_index, cutoff_index]);
        else
            keep_index = max([10, top_index, cutoff_index]);
        end
        top_paths{vn_i} = path_structs(sort_order(1:keep_index));
%         top_paths{vn_i} = cellfun(@(x) setfield(x, 'recur', 0), top_paths, 'unif', 0);
%         top_paths{vn_i} = cellfun(@(x) setfield(x, 'prob', 1), top_paths, 'UniformOutput', 0);
        prev_paths = top_paths{vn_i};

        if exist(sim_name, "file")
            if debug_flag
                save(sim_name, '-v7.3', '-append', 'top_paths')
            else
                save(sim_name, '-v7.3', '-append', 'top_paths', 'candidate_paths', 'chunk_likelihood', 'chunk_prior', 'hazards', 'haz_exp', 'haz_ent')
            end

        else
            if debug_flag
                save(sim_name, '-v7.3', 'top_paths')
            else
                save(sim_name, '-v7.3', 'top_paths', 'candidate_paths', 'chunk_likelihood', 'chunk_prior', 'hazards', 'haz_exp', 'haz_ent')
            end

        end

    end

end

% save(sim_name, '-v7.3', 'path_posterior', 'top_paths', 'candidate_paths', 'chunk_likelihood', 'chunk_prior', 'hazards',...
%     'haz_exp', 'haz_ent')

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

if nargin < 3, method = ''; end
if isempty(method), method = 'non-syllable-specific'; end

global stats

switch method

    case 'syllable-specific'

        [~, cdf_index] = cellfun(@(x) min(abs(x(:, 1) - duration)), stats.sylbs.cdf, 'unif', 0);
        hazards = cellfun(@(x, y) x(y, 2), stats.sylbs.cdf, cdf_index);
        haz_exp = sum(hazards.*nanunitsum(dist));
        haz_ent = -nansum(nanunitsum(hazards).*log(nanunitsum(hazards)))/log(length(dist)); % (log(num_phones) + nansum(hazard(:, w).*log(hazard(:, w))))/num_phones;

    case 'non-syllable-specific'

        this_cdf = stats.sylb_dur.cdf{1};
        [~, cdf_index] = min(abs(this_cdf(:, 1) - duration));
        haz_exp = this_cdf(cdf_index, 2);
        haz_ent = 0;
        hazards = haz_exp*ones(size(dist));

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