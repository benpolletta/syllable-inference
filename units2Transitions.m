function transition_dict = units2Transitions%(Tsylb_indicator)

wsp_map = load('word2sylb2phone_bysentence.mat');
wsp_map = wsp_map.results;

all_units = {cat(1, wsp_map.canonical_word_cell),...
    cat(2, wsp_map.sylb_cell), cat(2, wsp_map.phone_cell)};

unit_vocabs = cellfun(@unique, all_units, 'UniformOutput', false);
for v = 1:length(unit_vocabs)
    this_vocab = unit_vocabs{v};
    nonunits = cellfun(@(x) ismember(x, {'/','','#'}), this_vocab);
    this_vocab(nonunits) = [];
    this_vocab = [{'#'}, this_vocab(:)'];
    unit_vocabs{v} = this_vocab;
end
% unit_vocabs = cellfun(@(x) {'START', x{:}, 'END'}, unit_vocabs, 'UniformOutput', false);

subunits = {cat(2, wsp_map.word_sylb_cell),...
    cat(2, wsp_map.sylb_phone_cell)};

transition_dicts = cell(size(subunits));

for u = 1:length(subunits) %length(subunits) %
    curr_units = all_units{u};
    curr_subunits = subunits{u};
    curr_vocab = unit_vocabs{u+1};
    vocab_length = length(curr_vocab);
    %category = categorical(units{u});           
    %onehot{u} = onehotencode(category, 2);
    %cats = categories(category);
    curr_groups = findgroups(curr_units);
    num_groups = max(curr_groups);
    transition_dict = cell(1, num_groups);
    for g = 1:num_groups
        group_indices = find(curr_groups == g);
        group_subunits = curr_subunits(group_indices);
        group_size = length(group_indices);
        index = 1;
        % unit_transition = nan(vocab_length); %, vocab_length, group_size]);
        [I, J, Z] = deal([]);
        for j = 1:length(group_size)
            % transition_matrix = zeros(length(curr_vocab));
            these_subunits = group_subunits{j};
            these_subunits(cellfun(@(x) ismember(x, {'/','','#'}), these_subunits)) = [];
            if ~isempty(these_subunits)
                these_subunits = [{'#'}, these_subunits(:)'];
                num_subunits = length(these_subunits);
                [~, su_indices] = ismember(these_subunits, curr_vocab);
                % this_onehot = zeros(num_subunits, vocab_length);
                % idx = sub2ind([num_subunits, vocab_length], 1:num_subunits, su_indices);
                % this_onehot(idx) = 1;
                for k = 1:(num_subunits - 1)
                    I(index) = su_indices(k + 1);
                    J(index) = su_indices(k);
                    Z(index) = 1/group_size;
                    index = index + 1;
                    % transition_matrix =  transition_matrix + (this_onehot(k+1, :)')*this_onehot(k, :);
                end
                %imagesc(transition_matrix)
                %unit_transition = unit_transition + transition_matrix;
            end
        end
        transition_dict{g} = sparse(I, J, Z, vocab_length, vocab_length);% sparse(unit_transition/group_size);
    end
    transition_dicts{u} = dictionary(1:max(curr_groups), transition_dict);
    %transition_dict{u} = splitapply(@(x)mean(x{x}, 3), unit_transitions, curr_groups);
end

save('units2Transitions.mat', 'unit_vocabs', 'transition_dicts')