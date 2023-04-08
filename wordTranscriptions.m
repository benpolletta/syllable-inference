function transcription_data = wordTranscriptions
% Transcribes sentences in the TIMIT corpus with word boundaries marked by
% hashtags (#) - i.e., in a format readable by the syllabification program
% tsylb2.
% Writes word transcriptions (as sequences of phonemes) in files w/ extension .WRDphn.
% Produces wordTranscriptions.mat, which collects pronunciations by (unique) word.

SI = (1:6300)/6300;

time = 0:.1:10000;

onset_time = 1000;

name = 'Data';
if length(SI) ~= 6300
    name = sprintf('%s_%dsentences', name, length(SI));
end 

%results(length(SI)) = struct();
    
timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';
    
file_list_id = fopen([timit_dir, 'DOC/allphonelist_filenames.txt'], 'r');
file_list = textscan(file_list_id, '%s');
fclose(file_list_id);

file_list = file_list{1};

outfile_name = [timit_dir, 'DOC/allphonelist_tsylb_words.txt'];
outfile_fid = fopen(outfile_name, 'w');
fprintf(outfile_fid, '%s\n', '');

nxfile_name = [timit_dir, 'DOC/nx_words.txt'];
nxfile_fid = fopen(nxfile_name, 'w');

cwords = load('canonicalWords.mat');
cwords = cwords.results;
cwords = {cwords(:).canonical_words};

dict = load('word2sylb2phone_DICT.mat');
dict = dict.results;
dict_words = {dict(:).word};
dict_sylbs = {dict(:).sylb_cell};

[nx_prons, nx_dict_prons, resolutions] = deal({});

for s = 1:length(SI)
    
    file_index = SI(s);
    
    if isfloat(file_index) && file_index < 1 && file_index > 0
        file_index = round(file_index*length(file_list));
    end
    
    file_name = file_list{file_index};

    % if filename == 'TEST/DR8/MRES0/SI587', dbstop, end
    
    %% Retrieving words and their start and end times.
    
    word_filename = [timit_dir, file_name, '.WRD'];
    fid = fopen(word_filename, 'r');
    word_data = textscan(fid, '%s');
    fclose(fid);
    word_data = word_data{1};
    word_data = reshape(word_data, 3, length(word_data)/3);
    
    words = word_data(3, :);
    
    word_indices = cellfun(@str2num, word_data(1:2, :))';
    word_times = (word_indices/16 + onset_time);
    
    word_lengths = diff(word_times, [], 2); 
    
    %% Retrieving phonemes and their start and end times.
    
    phone_filename = [timit_dir, file_name, '.tsylbPHN'];
    fid = fopen(phone_filename, 'r');
    phone_data = textscan(fid, '%s');
    fclose(fid);
    phone_data = phone_data{1};
    phone_data = reshape(phone_data, 3, length(phone_data)/3);
    
    phones = phone_data(3, :);
    
    phone_indices = cellfun(@str2num, phone_data(1:2, :))';
    %mid_phone_indices = mean(phone_indices);
    
    phone_times = (phone_indices/16)';
    phone_lengths = diff(phone_times, [], 2);
    
    %% Retrieving & writing phonemic transcriptions.

    pronunciations = cell(size(word_lengths));

    pron_filename = [timit_dir, file_name, '.WRDphn'];
    fid = fopen(pron_filename, 'w');

    fprintf(outfile_fid, '%s ', '$##');

    these_cwords = cwords{s};
    
    for w = 1:length(word_lengths)
        
        this_word_indicator = all(phone_indices >= word_indices(w, 1) & phone_indices <= word_indices(w, 2), 2);
        
        this_pronunciation = phones(this_word_indicator);

        nx_index = find(cellfun(@(x) strcmp(x, 'nx'), this_pronunciation));

        if ~isempty(nx_index)

            dict_index = find(strcmp(these_cwords{w}, dict_words));
            these_dict_sylbs = dict_sylbs(dict_index);
            these_dict_phones = cellfun(@(x) strsplit(x, '/'), these_dict_sylbs{:}, 'unif', 0);
            these_dict_phones = cat(2, these_dict_phones{:});

            if ~isequal(this_pronunciation, these_dict_phones)

                for n = 1:length(nx_index)
                    
                    this_nxi = nx_index(n);

                    this_dict_pron = cellfun(@(x) strjoin(x, '*'), these_dict_sylbs, 'unif', 0);

                    this_pron = strjoin(this_pronunciation, '/');

                    nx_dict_phone = '';
                    if length(these_dict_phones) >= this_nxi
                        nx_dict_phone = these_dict_phones{this_nxi};
                    end

                    if strcmp(nx_dict_phone, 'n')
                        this_pronunciation{this_nxi} = 'n';
                        
                        nx_prons(end + 1) = {this_pron};
                        nx_dict_prons(end + 1) = this_dict_pron;
                        resolutions(end + 1) = {this_pronunciation};
                    elseif strcmp(nx_dict_phone, 'nx')
                        
                        nx_prons(end + 1) = {this_pron};
                        nx_dict_prons(end + 1) = this_dict_pron;
                        resolutions(end + 1) = {this_pronunciation};
                    end

                    nx_reverse_index = length(this_pronunciation) - this_nxi;
                    
                    nx_reverse_dict_phone = '';
                    if length(these_dict_phones) > nx_reverse_index
                        nx_reverse_dict_phone = these_dict_phones{end - nx_reverse_index};
                    end

                    if strcmp(nx_reverse_dict_phone, 'n')
                        this_pronunciation{end - nx_reverse_index} = 'n';

                        nx_prons(end + 1) = {this_pron};
                        nx_dict_prons(end + 1) = this_dict_pron;
                        resolutions(end + 1) = {this_pronunciation};
                    elseif strcmp(nx_reverse_dict_phone, 'nx') || (nx_reverse_index == 0 && strcmp(nx_reverse_dict_phone, 'm'))
                        
                        nx_prons(end + 1) = {this_pron};
                        nx_dict_prons(end + 1) = this_dict_pron;
                        resolutions(end + 1) = {this_pronunciation};
                    end

                    if ~strcmp(nx_dict_phone, 'nx') && ~strcmp(nx_reverse_dict_phone, 'nx') && ~strcmp(nx_dict_phone, 'n') && ~strcmp(nx_reverse_dict_phone, 'n')

                        prev_nx_index = find(strcmp(this_dict_pron{:}, nx_dict_prons) & strcmp(this_pron, nx_prons));

                        if ~isempty(prev_nx_index)

                            this_pronunciation = resolutions{prev_nx_index};

                        else

                            nx_prons(end + 1) = {this_pron};
                            nx_dict_prons(end + 1) = this_dict_pron;

                            prompt = sprintf('Sentence %d \n Canonical pronunciation: %s \n This pronunciation: %s \n Phoneme %d should be: \n 1) nx \n 2) n \n',...
                                s, this_dict_pron{:}, this_pron, this_nxi);

                            index = input(prompt);

                            if logical(index - 1)

                                this_pronunciation{this_nxi} = 'n';

                            end

                            resolutions(end + 1) = {this_pronunciation};

                        end

                    end

                    fprintf(nxfile_fid, '%s %s phone_num=%d %s\n', this_dict_pron{:}, this_pron, this_nxi, strjoin(this_pronunciation, '/'));

                end

            end

        end

        pronunciations{w} = {this_pronunciation};
        
        fprintf(fid, '%d %d %s\n', word_indices(w, :), strjoin(this_pronunciation(:), '/'));

        format = join(repmat({'%s '}, 1, length(this_pronunciation)), '');
        format = format{1}; % , '\n'];
        %fprintf(fid, [format, '\n'], this_pronunciation{:});

        fprintf(outfile_fid, format, this_pronunciation{:});

        if w < length(word_lengths)
            fprintf(outfile_fid, '%s ', '#');
        end
        
    end

    fprintf(outfile_fid, '%s\n', '##$');

    fclose(fid);
    
    %% Saving results.
    
    results(s) = struct('word_lengths', word_lengths, 'words', {words}, 'pronunciations', {pronunciations});
    
end

save('wordTranscriptions_results.m', 'results', 'nx_prons', 'nx_dict_prons', 'resolutions')

fclose(outfile_fid);

fclose(nxfile_fid);

word_vec = cat(2, results.words)';

pronunciation_vec_cell = cat(1, results.pronunciations);
pronunciation_vec_string = cellfun(@(x) strjoin(x{1}, '/'), pronunciation_vec_cell, 'unif', 0);


%% Collecting word lengths across sentences.

%%% Grouping by words.

[word_index, word_id] = findgroups(word_vec);

% Finding number of pronunciations.

pronunciation_map_cell = splitapply(@(x) {x}, pronunciation_vec_cell, word_index);
pronunciation_map_string = splitapply(@(x) {x}, pronunciation_vec_string, word_index);
[pronunciation_index, unique_pronunciations] = cellfun(@findgroups, pronunciation_map_string, 'unif', 0);
num_pronunciations = cellfun(@length, unique_pronunciations);

transcription_data = struct('results', {results}, 'word_id', {word_id}, 'unique_pronunciations', {unique_pronunciations}, 'num_pronunciations', {num_pronunciations});

save('wordTranscriptions.mat', 'transcription_data')

end

