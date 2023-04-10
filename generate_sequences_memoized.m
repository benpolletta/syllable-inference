function seqs = generate_sequences_memoized(current_seq, transition_matrix, n, cache)
% Generate all possible sequences of length n
% starting from the given current sequence
% using the given transition matrix and memoization cache.

% Check if current sequence is already in the cache
if isKey(cache, current_seq)
    seqs = cache(current_seq);
    return
end

% Check if current sequence is of length n
if length(current_seq) == n
    seqs = {current_seq}; % Return as a single sequence
    cache(current_seq) = seqs; % Store in the cache
    return
end

% Loop over possible symbols
seqs = {}; % Initialize sequences for this level
for i = 1:size(transition_matrix, 1)
    % Get current symbol
    symbol = char('A' + i - 1);
    % Check if current symbol can be appended to current sequence
    if isempty(current_seq) || transition_matrix(current_seq(end)-'A'+1, i) > 0
        % Append current symbol to current sequence
        next_seq = [current_seq symbol];
        % Recursive call to generate sequences with next sequence
        subseqs = generate_sequences_memoized(next_seq, transition_matrix, n, cache);
        % Add generated sequences to current sequences
        seqs = [seqs subseqs]; %#ok<AGROW>
    end
end

% Store result in the cache
cache(current_seq) = seqs;

end