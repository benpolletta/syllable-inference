function results = writeSyllables

wt_data = load('/projectnb/crc-nak/brpp/SI_Model/wordTranscriptions.mat');
wt_data = wt_data.transcription_data;
wt_data = wt_data.results;

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

cwords = load('canonicalWords.mat');
cwords = cwords.results;
cwords = {cwords(:).canonical_words};

dict = load('word2sylb2phone_DICT.mat');
dict = dict.results;
dict_words = {dict(:).word};
%dict_words = cellfun(@(x) erase(x, '-'), dict_words, 'unif', 0);
dict_sylbs = {dict(:).sylb_cell};

pound_fid = fopen([timit_dir, 'DOC/replaced_pauses.txt'], 'w');
ncsn_fid = fopen([timit_dir, 'DOC/non_canonical_sylb_num.txt'], 'w');
rncsn_fid = fopen([timit_dir, 'DOC/resolved_non_canonical_sylb_num.txt'], 'w');

for s = 1:length(tsylb_output)
    
    file_name = timit_filenames{s};
    phone_filename = [timit_dir, file_name, '.tsylbPHN'];
    [phone_cell, phone_indices, phone_durations] = getUnits(phone_filename); % sentence_phones(3, :);
    phone_cell(cellfun(@(x) contains(x, '$'), phone_cell)) = [];
    num_phones = length(phone_cell);
    %     phone_indices = cellfun(@str2num, sentence_phones(1:2, :))';
    phone_indices([1 end], :) = [];

    word_filename = [timit_dir, file_name, '.WRD'];
    [word_cell, word_indices, word_durations] = getUnits(word_filename);
   
    sentence_tsylb = tsylb_output{s}{1};

    %% Splitting into words.
    word_cell = strsplit(sentence_tsylb, '#');
    zero_length_indicator = @(x) length(x) == 0;
    empty_word_indicator = cellfun(@(x) length(x) == 0, word_cell);
    word_cell(empty_word_indicator) = [];
    word_cell = cellfun(@strip, word_cell, 'unif', 0);
    word_num = length(word_cell);
    
    %% Splitting words into syllables.
    word_sylb_cell = cellfun(@(x) strsplit(x, {'[', '] ', ']'}, 'CollapseDelimiters', true), word_cell, 'unif', false);
    word_sylb_nonzero_indicator = cellfun(@(x) cellfun(@(y) length(y) ~= 0, x), word_sylb_cell, 'unif', false);
    word_sylb_cell = cellfun(@(x, y) x(y), word_sylb_cell, word_sylb_nonzero_indicator, 'unif', 0);
    
    zero_sylb_words = cellfun(@(x) ~contains(x, {'[', ']'}), word_cell);
    % zero_sylb_words = cellfun(@(x) length(x) == 0, word_sylb_extract);
    %     word_sylb_cell = cellfun(@strip, word_sylb_cell, 'unif', 0);

    %% Splitting syllables into phones.
    word_sylb_phone_cell = cellfun(@(x) cellfun(@(y) strsplit(y, {' ','''0'}, 'CollapseDelimiters', true), x, 'unif', 0), word_sylb_cell, 'unif', 0);
    word_sylb_phone_cell = word_sylb_phone_cell(~cellfun(@isempty, word_sylb_phone_cell));
    word_sylb_phone_nonzero_indicator = cellfun(@(x) cellfun(@(y) cellfun(@(z) length(z) ~= 0, y), x, 'unif', 0), word_sylb_phone_cell, 'unif', 0);
    word_sylb_phone_cell = cellfun(@(x, y) cellfun(@(z, w) z(w), x, y, 'unif', 0), word_sylb_phone_cell, word_sylb_phone_nonzero_indicator, 'unif', 0);

    %% Adding in words with zero durations.
    zero_duration_indices = [];
    if any(word_durations == 0)
        zero_duration_indices = find(word_durations == 0);
        for i = 1:length(zero_duration_indices)
%             word_sylb_cell((zero_duration_indices(i):end) + 1) = word_sylb_cell(zero_duration_indices(i):end);
%             word_sylb_cell{(zero_duration_indices(i))} = {''};
            word_sylb_phone_cell((zero_duration_indices(i):end) + 1) = word_sylb_phone_cell(zero_duration_indices(i):end);
            word_sylb_phone_cell{(zero_duration_indices(i))} = {{''}};
        end
    end

    %% Counting syllables & phones.
    word_sylb_phone_num = cellfun(@(x) cellfun(@length, x), word_sylb_phone_cell, 'unif', 0);
    word_phone_num = cellfun(@sum, word_sylb_phone_num);
    word_sylb_num = cellfun(@length, word_sylb_cell);
    word_sylb_num(zero_sylb_words) = 0;
    word_sylb_num(zero_duration_indices) = 0;
    word_sylb_num(word_phone_num == 0) = 0;

    %% Fixing words containing pauses transcribed as #.

    transcript = wt_data(s).pronunciations;
    transcript = transcript(~cellfun(@(x) cellfun(@(y) isempty(y), x), transcript));
    t_words = cat(2, transcript{:});
    t_pound_count = cellfun(@(x) any(strcmp(x, '#')), t_words);

    if any(t_pound_count)

        t_phone_count = cellfun(@length, t_words);
        t_comp_count = t_phone_count - t_pound_count;
        cum_t_count = cumsum([0 t_comp_count]);

        wsp_words = cellfun(@(x) cat(2, x{:}), word_sylb_phone_cell, 'unif', 0);
        wsp_phones = cat(2, wsp_words{:});
        cum_w_count = cumsum(word_phone_num);

        if length(t_words) == word_num - 1

            %% Assuming that a word is split in two by the #.

            for i = 1:length(transcript)

                this_word_indicator = cum_w_count <= cum_t_count(i + 1) & cum_w_count > cum_t_count(i);
                wsp2t_map(this_word_indicator) = i;
                zero_sylb_words(i) = all(zero_sylb_words(this_word_indicator));

            end

            pound_adjacent = wsp2t_map == find(t_pound_count);
            pa_indices = find(pound_adjacent);

            % Recombining split words;
            split_word = cat(2, word_sylb_phone_cell(pound_adjacent));
            new_syl = {split_word{1}{end}, {'#'}, split_word{2}{1}};
            new_syl = cat(2, new_syl{:});
            new_word = {split_word{1}{1:(end - 1)}, new_syl, split_word{2}{2:end}};

            % Replacing split words with new word.
            word_sylb_phone_cell{pa_indices(1)} = new_word;
            word_sylb_phone_cell(pa_indices(2)) = [];

        elseif length(t_words) ~= word_num

            display(sprintf('Error w/ %s: transcript has %d words & word_sylb_phone_cell has %d words.', file_name, length(t_words), word_num))

        elseif any(t_phone_count - word_phone_num == 1)

            % Finding word that's affected.
            pound_adjacent = t_phone_count - word_phone_num == 1;
            pa_index = find(pound_adjacent);
            poundless_word = word_sylb_phone_cell(pound_adjacent);

            % Finding phone indices of that word.
            word_phone_indices = mat2cell(1:length(wsp_phones), 1, word_phone_num);
            word_sylb_phone_indices = cellfun(@(x, y) mat2cell(x, 1, y), word_phone_indices, word_sylb_phone_num, 'unif', 0);
            poundless_indices = word_sylb_phone_indices{pound_adjacent};

            t_phones = cat(2, t_words{:});
            pound_index = find(strcmp(t_phones, '#'), 1);
            
            % Adding phone to beginning or end of word.
            new_word = poundless_word;
            if pound_index < min(poundless_indices{:})
                new_word = cat(2, {'#'}, new_word{1}{:});
            elseif pound_index > max(poundless_indices{:})
                new_word = cat(2, new_word{1}{:}, {'#'});
            end

            % Replacing poundless word with new word.
            word_sylb_phone_cell{pa_index} = {new_word};

        elseif any(t_phone_count ~= word_phone_num)

            display(sprintf('Error w/ %s: transcript phones per word:', file_name))
            t_phone_count
            display('word_sylb_phone_cell phones per word:')
            word_phone_num

        end

        % Recounting phones per syllable and word.
        word_sylb_phone_num = cellfun(@(x) cellfun(@length, x), word_sylb_phone_cell, 'unif', 0);
        word_phone_num = cellfun(@sum, word_sylb_phone_num);


%         % Adjusting pound_index.
%         pound_phone_index = pound_phone_index + (cumsum(ones(size(pound_phone_index))) > p);

        % Printing out fix.
        for_print = cellfun(@(x) cellfun(@(y) strjoin(y, '/'), x, 'unif', 0), word_sylb_phone_cell, 'unif', 0);
        for_print = cellfun(@(x) strjoin(x, '*'), for_print, 'unif', 0);
        for_print = strjoin(for_print,'|');

        fprintf(pound_fid, '%s\n%s\n%s\n', file_name, sentence_tsylb, for_print);

    end

    %% Rejoining phones.
    word_sylb_cell = cellfun(@(x) cellfun(@(y) strjoin(y, '/'), x, 'unif', 0), word_sylb_phone_cell, 'unif', 0);
    word_cell = cellfun(@(x) strjoin(x, '*'), word_sylb_cell, 'unif', 0);

    %% Comparing to canonical words to fix syllabification errors.
    these_cwords = cwords{s};
    
    for w = 1:length(word_cell)

        % Finding canonical pronunciations.

        dict_index = find(strcmp(these_cwords{w}, dict_words));
        these_dict_sylbs = dict_sylbs(dict_index);
        these_dict_phones = cellfun(@(x) strsplit(x, '/'), these_dict_sylbs{:}, 'unif', 0);
        these_dict_phones = cat(2, these_dict_phones{:});
        these_dict_prons = cellfun(@(x) strjoin(x, '*'), these_dict_sylbs, 'unif', 0);
        initial_pron = word_cell{w};
        initial_phones = cat(2, word_sylb_phone_cell{w}{:});
        resolution = '';

        if length(these_dict_sylbs{:}) == length(word_sylb_cell{w})

            canonical_sylbs{w} = these_dict_sylbs;

        else
            
            if all(strcmp(these_dict_sylbs{:}, 'y/ih/r'))

                canonical_sylbs{w} = word_sylb_cell{w};

                resolution = strjoin(word_sylb_cell{w}, '*');

            elseif length(these_dict_phones) == length(initial_phones)

                sound_families = {{'q', 't', 'd', 'dh', 'dx'}, {'s', 'z', 'zh', 'jh'}};
                delimiters = {'/','*'};
                strings = {these_dict_prons{:}, word_cell{w}};
                [patterns, locations] = deal(cell(length(delimiters), length(sound_families)));

                for d = 1:length(delimiters)
                    patterns(d, :) = cellfun(@(y)['(', strjoin(cellfun(@(x) [delimiters{d},x], y, 'unif', 0), '|'), ')'], sound_families, 'unif', 0);
                    locations(d, :) = cellfun(@(x) regexp(strings{d}, x), patterns(d, :), 'unif', 0);
                end
                
                family = find(all(~cellfun(@isempty, locations)));

                if ~isempty(family)

                    fam_home = locations(:, family);

                    if fam_home{1} == fam_home{2}

                        word_cell{w}(fam_home{2}) = '/';

                        word_sylb_cell{w} = strsplit(word_cell{w}, '*');
                        word_sylb_phone_cell{w} = cellfun(@(x) strsplit(x, '/'), word_sylb_cell{w}, 'unif', 0);

                        canonical_sylbs{w} = these_dict_sylbs;

                        resolution = strjoin(word_sylb_cell{w}, '*');

                    end

                end

            end

            if ~isempty(resolution)

                fprintf(rncsn_fid, '%s %s %s %s\n', file_name, strjoin(these_dict_sylbs{:}, '*'), initial_pron, resolution);

            else

                fprintf(ncsn_fid, '%s %s %s\n', file_name, strjoin(these_dict_sylbs{:}, '*'), initial_pron);

            end

        end


    end


    %% Recounting syllables per word.
    word_sylb_num = cellfun(@length, word_sylb_cell);
    word_sylb_num(zero_sylb_words) = 0;
    word_sylb_num(zero_duration_indices) = 0;
    word_sylb_num(word_phone_num == 0) = 0;

    %% Aggregating into syllables.
    sylb_phone_cell = cat(2, word_sylb_phone_cell{:});
    sylb_cell = cat(2, word_sylb_cell{:});
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
    sylb_start_phone = min([0, sylb_end_phone((1:(end - 1)))] + 1, num_phones);
    %sylb_phone_bounds = [sylb_start_phone, sylb_end_phone];
    
    sylb_start_indices = phone_indices(sylb_start_phone, 1);
    sylb_end_indices = phone_indices(sylb_end_phone, 2);

    for_print = mat2cell([sylb_start_indices, sylb_end_indices], ones(length(sylb_start_indices), 1), [1 1]);
    for_print(:, end + 1) = sylb_cell(:);
    for_print = for_print';

    %     %% Removing words that have no syllables.
    %     word_sylb_sum = cumsum(word_sylb_num);
    %     word_zero_sylbs = word_sylb_sum(word_sylb_num == 0) + 1;
    %     for_print(:, word_zero_sylbs) = [];


    %% Writing syllables.
    
    sylb_filename = [timit_dir, timit_filenames{s}, '.SYLB'];
    fid = fopen(sylb_filename, 'w');
    fprintf(fid, '%d %d %s\n', for_print{:});
    fclose(fid);

    %% Saving results.
    results(s) = struct('word_cell', {word_cell}, 'word_sylb_cell', {word_sylb_cell}, 'word_sylb_phone_cell', {word_sylb_phone_cell},...
        'word_sylb_phone_num', {word_sylb_phone_num}, 'word_phone_num', word_phone_num, 'word_sylb_num', word_sylb_num,...
        'sylb_phone_cell', {sylb_phone_cell}, 'sylb_cell', {sylb_cell}, 'phone_cell', {phone_cell});
    
end

fclose(pound_fid);

save('word2sylb2phone_bysentence.mat', 'results')

end


function [units, indices, durations] = getUnits(filename)

onset_time = 1000;

fid = fopen(filename, 'r');
data = textscan(fid, '%s');
fclose(fid);
data = data{1};
data = reshape(data, 3, length(data)/3);

units = data(3, :);

indices = cellfun(@str2num, data(1:2, :))';
% times = (indices/16 + onset_time);

durations = diff(indices, [], 2);

end
