function seqs = generate_word_sequences(current_seq, symbols, lengths, transition_matrix, n, cutoff)
% Generate all possible sequences of length n
% starting from the given current sequence
% using the given transition matrix.

if nargin < 5, cutoff = []; end
if isempty(current_seq), current_seq.seq = {}; current_seq.length = 0; end
if isempty(cutoff), cutoff = 0; end

% Check if current sequence is of length n
if current_seq.length >= n
    seqs = {current_seq}; % Return as a single sequence
    return
end

% Loop over possible symbols
seqs = {};
if isempty(current_seq.seq)
    for i = 1:length(symbols)
        % Get current symbol
        this_symbol = symbols{i};
        this_length = lengths(i);
        % Append current symbol to current sequence
        next_seq.seq = [current_seq.seq this_symbol];
        next_seq.length = current_seq.length + this_length;
        % Recursive call to generate sequences with next sequence
        subseqs = generate_word_sequences(next_seq, symbols, lengths, transition_matrix, n, cutoff);
        % Add generated sequences to current sequences
        seqs = [seqs subseqs]; %#ok<AGROW>
    end
else
    % Check which symbols can be appended to current sequence
    possible_indices = find(transition_matrix(:, strcmp(symbols, current_seq.seq(end))) > cutoff);
    for i = 1:length(possible_indices)
        % Get current symbol & length
        this_symbol = symbols{i};
        this_length = lengths(i);
        % Append current symbol to current sequence
        next_seq.seq = [current_seq.seq this_symbol];
        next_seq.length = current_seq.length + this_length;
        % Recursive call to generate sequences with next sequence
        subseqs = generate_word_sequences(next_seq, symbols, lengths, transition_matrix, n, cutoff);
        % Add generated sequences to current sequences
        seqs = [seqs subseqs]; %#ok<AGROW>
    end
end

end