function seqs = generate_sequences(current_seq, symbols, transition_matrix, n, cutoff)
% Generate all possible sequences of length n
% starting from the given current sequence
% using the given transition matrix.

if nargin < 5, cutoff = []; end
if isempty(cutoff), cutoff = 0; end

% Check if current sequence is of length n
if length(current_seq) == n
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
        next_seq = [current_seq symbol];
        % Recursive call to generate sequences with next sequence
        subseqs = generate_sequences(next_seq, symbols, transition_matrix, n, cutoff);
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
        next_seq = [current_seq symbol];
        % Recursive call to generate sequences with next sequence
        subseqs = generate_sequences(next_seq, symbols, transition_matrix, n, cutoff);
        % Add generated sequences to current sequences
        seqs = [seqs subseqs]; %#ok<AGROW>
    end
end

end