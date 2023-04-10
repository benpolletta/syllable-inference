function sequences = enumerate_sequences_recursive(symbols, transition_matrix, n)
% ENUMERATE_SEQUENCES_RECURSIVE Enumerate all possible sequences of symbols of
%   length n given a list of possible symbols and a transition matrix using
%   recursive backtracking.
%   sequences = ENUMERATE_SEQUENCES_RECURSIVE(symbols, transition_matrix, n)
%   returns a cell array of all possible sequences of symbols of length n,
%   where each sequence is represented as a row vector of indices into the
%   symbols list.

% Check input arguments
if ~iscell(symbols) || ~isnumeric(transition_matrix) || ~isscalar(n)
    error('Invalid input arguments');
end

% Initialize sequences array
sequences = {};

% Recursively generate sequences
generate_sequence([], 0);

    function generate_sequence(sequence, sequence_length)
        % Base case: sequence is of length n
        if sequence_length == n
            sequences{end+1} = sequence;  % Add sequence to list
            return;
        end
        
        % Recursive case: generate next symbol and continue generating sequence
        num_symbols = numel(symbols);
        prev_symbol_idx = sequence(end);
        transition_probs = transition_matrix(prev_symbol_idx, :);
        
        for i = 1:num_symbols
            next_symbol_idx = i;
            next_sequence = [sequence next_symbol_idx];
            generate_sequence(next_sequence, sequence_length+1);
        end
    end

end
