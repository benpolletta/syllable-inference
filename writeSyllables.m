function results = writeSyllables

timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';
tsylb_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/tsylb2-1.1/';

timit_filenames_file = [timit_dir, 'DOC/allphonelist_filenames.txt'];
filename_fid = fopen(timit_filenames_file, 'r');
timit_filenames = textscan(filename_fid, '%s', 'Delimiter', '\n');
fclose(filename_fid);
timit_filenames = timit_filenames{1};

tsylb_output_file = [timit_dir, 'DOC/allphonelist_tsylb_words_out_l43897fix.txt'];
tsylb_fid = fopen(tsylb_output_file, 'r');
tsylb_output = textscan(tsylb_fid, '%s', 'headerLines', 1, 'Delimiter', '\n');
fclose(tsylb_fid);
tsylb_output = tsylb_output{1};

output_line_indicator = contains(tsylb_output, 'Enter');
tsylb_output(~output_line_indicator) = [];
tsylb_output(end) = [];
%tsylb_output(cellfun(@isempty, tsylb_output)) = [];
tsylb_output = cellfun(@(x) extractBetween(x, 'Basic pron is /', '/'), tsylb_output, 'unif', 0); % extractBetween(x, '/#', '#/'), 

for s = 1:length(tsylb_output)
    
    sentence_filename = [timit_dir, timit_filenames{s}, '.tsylbPHN'];
    fid = fopen(sentence_filename, 'r');
    sentence_phones = textscan(fid, '%s');
    fclose(fid);
    sentence_phones = sentence_phones{1};
    sentence_phones = reshape(sentence_phones, 3, length(sentence_phones)/3);
    phone_cell = sentence_phones(3, :);
    phone_cell(cellfun(@(x) contains(x, '#'), phone_cell)) = [];
    num_phones = length(phone_cell);
    phone_indices = cellfun(@str2num, sentence_phones(1:2, :))';
    phone_indices([1 end], :) = [];
   
    sentence_tsylb = tsylb_output{s}{1};

    %% Splitting into words.
    word_cell = strsplit(sentence_tsylb, '#');
    zero_length_indicator = @(x) length(x) == 0;
    empty_word_indicator = cellfun(@(x) length(x) == 0, word_cell);
    word_cell(empty_word_indicator) = [];
    word_cell = cellfun(@strip, word_cell, 'unif', 0);
    
    %% Splitting words into syllables.
    word_sylb_cell = cellfun(@(x) extractBetween(x, '[', ']'), word_cell, 'unif', false); 
    word_sylb_cell = cellfun(@strip, word_sylb_cell, 'unif', 0);

    %% Splitting syllables into phones.
    word_sylb_phone_cell = cellfun(@(x) cellfun(@(y) strsplit(y, {' ','''0'}, 'CollapseDelimiters', true), x, 'unif', 0), word_sylb_cell, 'unif', 0);
    word_sylb_phone_nonzero_length = cellfun(@(x) cellfun(@(y) cellfun(@(z) length(z) ~= 0, y), x, 'unif', 0), word_sylb_phone_cell, 'unif', 0);
    word_sylb_phone_cell = cellfun(@(x, y) cellfun(@(z, w) z(w), x, y, 'unif', 0), word_sylb_phone_cell, word_sylb_phone_nonzero_length, 'unif', 0);
    
    %% Rejoining phones.
    word_sylb_cell = cellfun(@(x) cellfun(@(y) strjoin(y, '/'), x, 'unif', 0), word_sylb_phone_cell, 'unif', 0);
    word_cell = cellfun(@(x) strjoin(x, '*'), word_sylb_cell, 'unif', 0);

    %% Counting syllables & phones.
    word_sylb_phone_num = cellfun(@(x) cellfun(@length, x), word_sylb_phone_cell, 'unif', 0);
    word_sylb_num = cellfun(@length, word_sylb_cell);

    %% Aggregating into syllables.
    sylb_phone_cell = cat(1, word_sylb_phone_cell{:});
    sylb_cell = cat(1, word_sylb_cell{:});
    % %%Splitting into syllables.
    %     sylb_cell = strsplit(sentence_tsylb, {'# [', '[', '] ', '#'}, 'CollapseDelimiters', true);
    %
    %     % Splitting syllables into phones.
    %     % sylb_split = cellfun(@(x) textscan(x, '%s', 'delimiter', {' ', '''0'}), sentence_split, 'unif', false); %
    %     sylb_phone_cell = cellfun(@(x) strsplit(x, {' ','''0'}, 'CollapseDelimiters', true), sylb_cell, 'unif', 0);
    %     sylb_phone_nonzero_length = cellfun(@(x) cellfun(@(y) length(y) ~= 0, x), sylb_phone_cell, 'unif', 0);
    %     sylb_phone_cell = cellfun(@(x, y) x(y), sylb_phone_cell, sylb_phone_nonzero_length, 'unif', 0);
    %     sylb_phone_cell = sylb_phone_cell(~cellfun(@isempty, sylb_phone_cell));
    
    sylb_phone_num = cellfun(@length, sylb_phone_cell);

    %% Retrieving start and end indices of syllables.
    
    
    sylb_end_phone = min(cumsum(sylb_phone_num), num_phones);
    sylb_start_phone = min([0; sylb_end_phone((1:(end - 1)))] + 1, num_phones);
    %sylb_phone_bounds = [sylb_start_phone, sylb_end_phone];
    
    sylb_start_indices = phone_indices(sylb_start_phone, 1);
    sylb_end_indices = phone_indices(sylb_end_phone, 2);

    for_print = mat2cell([sylb_start_indices, sylb_end_indices], ones(length(sylb_start_indices), 1), [1 1]);
    for_print(:, end + 1) = sylb_cell(:);
    for_print = for_print';
    
    boundary_filename = [timit_dir, timit_filenames{s}, '.SYLB'];
    fid = fopen(boundary_filename, 'w');
    fprintf(fid, '%d %d %s\n', for_print{:});
    fclose(fid);

    %% Saving results.
    results(s) = struct('word_cell', {word_cell}, 'word_sylb_cell', {word_sylb_cell}, 'word_sylb_phone_cell', {word_sylb_phone_cell},...
        'sylb_phone_cell', {sylb_phone_cell}, 'sylb_cell', {sylb_cell}, 'phone_cell', {phone_cell});
    
end

save('word2sylb2phone_bysentence.mat', 'results')