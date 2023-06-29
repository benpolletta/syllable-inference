function seqs = generate_sequences(current_seq, symbols, transition_matrix, n, cutoff, length_function)
% Generate all possible sequences of length n
% starting from the given current sequence
% using the given transition matrix.

if nargin < 5, cutoff = []; end
if isempty(cutoff), cutoff = 0; end
if nargin < 6, length_function = []; end
if isempty(length_function)
    [tsylb_phonemes, class_indicator, class_names] = getPhones(1);
    vowels = tsylb_phonemes(class_indicator(:, strcmpi(class_names, 'vowels')));
    vowels = {vowels{:}, 'el', 'em', 'en', 'enx'};
    length_function = @(x) get_vowel_onsets(x, vowels); 
end

% Check if current sequence is of length n
if feval(length_function, current_seq) >= n
    seqs = {current_seq}; % Return as a single sequence
    return
end

% Loop over possible symbols
seqs = {};
if isempty(current_seq)
    for i = 1:length(symbols)
        % Get current symbol
        symbol = symbols{i};
        % Append current symbol to current sequence
        next_seq = {current_seq{:}, symbol};
        % Recursive call to generate sequences with next sequence
        subseqs = generate_sequences(next_seq, symbols, transition_matrix, n, cutoff, length_function);
        % Add generated sequences to current sequences
        seqs = [seqs subseqs]; %#ok<AGROW>
    end
else
    % Check which symbols can be appended to current sequence
    possible_indices = find(transition_matrix(:, strcmp(symbols, current_seq(end))) > cutoff);
    for i = 1:length(possible_indices)
        % Get current symbol
        symbol = symbols{possible_indices(i)};
        % Append current symbol to current sequence
        next_seq = {current_seq{:}, symbol};
        % Recursive call to generate sequences with next sequence
        subseqs = generate_sequences(next_seq, symbols, transition_matrix, n, cutoff, length_function);
        % Add generated sequences to current sequences
        seqs = [seqs subseqs]; %#ok<AGROW>
    end
end

end

function [vowel_length, vowel_onsets, phones] = get_vowel_onsets(chunk, vowels)

% Turning syllables into sequences of phonemes.
phones = cellfun(@(y) split(y, '/'), chunk, 'unif', 0);
phones = cat(1, phones{:});

% Finding vocalic nuclei & getting length (in vocalic nuclei).
vowel_length = 0;
vowel_onsets = [];

if length(phones) > 0

    empty_phones = cellfun(@isempty, phones);
    phones = phones(~empty_phones);

    vowel_indicators = cellfun(@(x) any(strcmp(vowels, x)), phones);
    vowel_onsets = find(diff([0; vowel_indicators]) == 1);

    vowel_length = length(vowel_onsets);

end

end