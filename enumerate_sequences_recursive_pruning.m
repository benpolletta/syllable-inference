function sequences = enumerate_sequences(symbols, transition_matrix, n)
% Enumerate all possible sequences of symbols of length n
% given a list of possible symbols and a transition matrix
% using recursion and pruning to minimize computation time.

% Initialize sequences with empty sequence
sequences = {''};

% Recursive function to generate sequences
function generate_sequences(current_seq)
    % Check if sequence is already of length n
    if length(current_seq) == n
        sequences{end+1} = current_seq; % Add to list of sequences
        return
    end
    % Loop over possible symbols
    for i = 1:length(symbols)
        % Get current symbol
        symbol = symbols(i);
        % Check if current symbol can be appended to current sequence
        if isempty(current_seq) || transition_matrix(current_seq(end), symbol) > 0
            % Append current symbol to current sequence
            next_seq = [current_seq symbol];
            % Recursive call to generate sequences with next sequence
            generate_sequences(next_seq);
        end
    end
end

% Generate sequences starting from empty sequence
generate_sequences('');

end
