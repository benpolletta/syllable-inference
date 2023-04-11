function [seqs, cache] = generate_sequences_memoized(current_seq, symbols, transition_matrix, n, cache)
% Generate all possible sequences of length n
% starting from the given current sequence
% using the given transition matrix and memoization cache.
if nargin < 5, cache = []; end
if isempty(cache), cache = containers.Map(); end

% Check if current sequence is already in the cache
if isKey(cache, strjoin(current_seq, '#'))
    seqs = cache(strjoin(current_seq, '#'));
    return
end

% Check if current sequence is of length n
if length(current_seq) == n
    seqs = {current_seq}; % Return as a single sequence
    cache(strjoin(current_seq, '#')) = seqs; % Store in the cache
    return
end

% Loop over possible symbols
seqs = {};
for i = 1:length(symbols)
    % Get current symbol
    symbol = symbols{i};
    % Check if current symbol can be appended to current sequence
    if isempty(current_seq) || transition_matrix(strcmp(symbols, current_seq(end)), i) > 0
        % Append current symbol to current sequence
        next_seq = [current_seq symbol];
        % Recursive call to generate sequences with next sequence
        subseqs = generate_sequences_memoized(next_seq, symbols, transition_matrix, n, cache);
        % Add generated sequences to current sequences
        seqs = [seqs subseqs]; %#ok<AGROW>
    end
end
% Add sequence to memoization cache
cache(strjoin(current_seq, '#')) = seqs;

end