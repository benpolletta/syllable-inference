function sequences = enumerate_sequences(symbols, transition_matrix, n)
% ENUMERATE_SEQUENCES Enumerate all possible sequences of symbols of length n
%   given a list of possible symbols and a transition matrix.
%   sequences = ENUMERATE_SEQUENCES(symbols, transition_matrix, n) returns a
%   cell array of all possible sequences of symbols of length n, where each
%   sequence is represented as a row vector of indices into the symbols list.
%   The transition matrix specifies the probabilities of transitioning from
%   one symbol to another in the sequence.

% Check input arguments
if ~iscell(symbols) || ~isnumeric(transition_matrix) || ~isscalar(n)
    error('Invalid input arguments');
end

% Initialize sequences array
num_symbols = numel(symbols);
sequences = cell(num_symbols^n, 1);

% Enumerate all possible sequences
sequence_idx = 1;
for i1 = 1:num_symbols
    for i2 = 1:num_symbols
        % Initialize sequence
        sequence = [i1];
        
        % Generate remaining symbols
        for j = 2:n
            % Get transition probabilities for previous symbol
            prev_symbol_idx = sequence(j-1);
            transition_probs = transition_matrix(prev_symbol_idx, :);
            
            % Choose next symbol based on transition probabilities
            next_symbol_idx = randsample(num_symbols, 1, true, transition_probs);
            sequence = [sequence next_symbol_idx];
        end
        
        % Add sequence to list
        sequences{sequence_idx} = sequence;
        sequence_idx = sequence_idx + 1;
    end
end

% Remove unused cells from sequences array
sequences(sequence_idx:end) = [];
