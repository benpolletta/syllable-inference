function sequences = enumerate_sequences_recursive_pruning_memo(symbols, transition_matrix, n)
% Enumerate all possible sequences of symbols of length n
% given a list of possible symbols and a transition matrix
% using recursion and pruning to minimize computation time,
% and memoization to avoid redundant calculations.

% Initialize memoization cache
cache = containers.Map();

% Recursive function to generate sequences
function seq = generate_sequences(current_seq)
    % Check if sequence is already of length n
    if length(current_seq) == n
        seq = {current_seq}; % Return as cell array
        return
    end
    % Check if sequence has already been computed
    if isKey(cache, current_seq)
        seq = cache(current_seq); % Return from cache
        return
    end
    % Loop over possible symbols
    seq = {};
    for i = 1:length(symbols)
        % Get current symbol
        symbol = symbols(i);
        % Check if current symbol can be appended to current sequence
        if isempty(current_seq) || transition_matrix(current_seq(end), symbol) > 0
            % Append current symbol to current sequence
            next_seq = [current_seq symbol];
            % Recursive call to generate sequences with next sequence
            next_seq_cell = generate_sequences(next_seq);
            % Add new sequences to list of sequences
            seq = [seq, next_seq_cell];
        end
    end
    % Add sequence to memoization cache
    cache(current_seq) = seq;
end

% Generate sequences starting from empty sequence
sequences = generate_sequences('');

end
