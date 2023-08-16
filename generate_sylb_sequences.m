function seqs = generate_sylb_sequences(current_seq, symbols, probs, transition_matrix, n, cutoff, max_recur, save_opt)
% Generate all possible sequences of length n
% starting from the given current sequence
% using the given transition matrix.

if isempty(current_seq), current_seq = struct('sylbs', {{}}, 'phones', {{}}, 'recur', 0, 'prob', 1, 'vowel_indicator', [], 'num_onsets', 0); end
if nargin < 6, cutoff = []; end
if isempty(cutoff), cutoff = 0; end
if nargin < 7, max_recur = []; end
if isempty(max_recur), max_recur = inf; end
if nargin < 8, save_opt = []; end
if isempty(save_opt), save_opt = 0; end

[tsylb_phonemes, class_indicator, class_names] = getPhones(1);
vowels = tsylb_phonemes(class_indicator(:, strcmpi(class_names, 'vowels')));
vowels = {vowels{:}, 'el', 'em', 'en', 'enx'};

% Check if current sequence is of length n
if current_seq.num_onsets >= n
    seqs = {current_seq}; % Return as a single sequence
    return
end

% Loop over possible symbols
seqs = {};
if isempty(current_seq.sylbs)
    for i = 1:length(symbols)
        % Get current symbol
        this_symbol = symbols(i);
        [this_indicator, this_phones] = get_vowel_onsets(this_symbol, vowels);
        % Append current symbol to current sequence
        next_seq.sylbs = {current_seq.sylbs{:}, this_symbol{:}};
        next_seq.recur = 1;
        next_seq.prob = 1;
        next_seq.phones = {current_seq.phones{:}, this_phones{:}};
        next_seq.vowel_indicator = [current_seq.vowel_indicator, this_indicator];
        vowel_onsets = find(diff(next_seq.vowel_indicator) == 1);
        next_seq.num_onsets = length(vowel_onsets);
        % Recursive call to generate sequences with next sequence
        subseqs = generate_sylb_sequences(next_seq, symbols, probs, transition_matrix, n, cutoff);
        % Add generated sequences to current sequences
        seqs = [seqs; subseqs]; %#ok<AGROW>
    end
else
    % Check which symbols can be appended to current sequence
    current_index = find(strcmp(symbols, current_seq.sylbs(end)));
    possible_indices = find(transition_matrix(:, current_index) > cutoff);
    % fprintf('Possible indices: %d\n', length(possible_indices))
    for i = 1:length(possible_indices)
        % Get current symbol & length
        this_index = possible_indices(i);
        this_symbol = symbols(this_index);
        [this_indicator, this_phones] = get_vowel_onsets(this_symbol, vowels);
        % Append current symbol to current sequence
        next_seq.sylbs = {current_seq.sylbs{:}, this_symbol{:}};
        next_seq.recur = current_seq.recur + 1;
        next_seq.prob = current_seq.prob*transition_matrix(current_index, this_index);
        next_seq.phones = {current_seq.phones{:}, this_phones{:}};
        next_seq.vowel_indicator = [current_seq.vowel_indicator; this_indicator];
        vowel_onsets = find(diff(next_seq.vowel_indicator) == 1);
        next_seq.num_onsets = length(vowel_onsets);
        % Recursive call to generate sequences with next sequence
        if next_seq.recur < max_recur % next_seq.prob^(1/next_seq.recur) > 10^-3
            subseqs = generate_sylb_sequences(next_seq, symbols, probs, transition_matrix, n, cutoff);
            % Add generated sequences to current sequences
            seqs = [seqs; subseqs]; %#ok<AGROW>
        else
            return
        end
    end
end

if save_opt

    datestr = datetime('now', 'Format', 'yy-MM-dd_HH-mm-ss');
    name = sprintf('sylbSequences_%d%g_%s.mat', n, cutoff, datestr);

    save(name)

end

end


function [vowel_indicator, phones, vowel_length] = get_vowel_onsets(chunk, vowels)

% Turning syllables into sequences of phonemes.
phones = cellfun(@(y) split(y, '/'), chunk, 'unif', 0);
phones = cat(1, phones{:});

% Finding vocalic nuclei & getting length (in vocalic nuclei).
vowel_length = 0;
vowel_onsets = [];

if length(phones) > 0

    empty_phones = cellfun(@isempty, phones);
    phones = phones(~empty_phones);

    vowel_indicator = cellfun(@(x) any(strcmp(vowels, x)), phones);
    vowel_onsets = find(diff(vowel_indicator) == 1);

    vowel_length = length(vowel_onsets);

end

end