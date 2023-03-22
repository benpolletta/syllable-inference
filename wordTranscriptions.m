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
    
sentence_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';
    
file_list_id = fopen([sentence_dir, 'DOC/allphonelist_filenames.txt'], 'r');
file_list = textscan(file_list_id, '%s');
fclose(file_list_id);

file_list = file_list{1};

outfile_name = [sentence_dir, 'DOC/allphonelist_tsylb_words.txt'];
outfile_fid = fopen(outfile_name, 'w');

for s = 1:length(SI)
    
    file_index = SI(s);
    
    if isfloat(file_index) && file_index < 1 && file_index > 0
        file_index = round(file_index*length(file_list));
    end
    
    file_name = file_list{file_index};

    % if filename == 'TEST/DR8/MRES0/SI587', dbstop, end
    
    %% Retrieving words and their start and end times.
    
    word_filename = [sentence_dir, file_name, '.WRD'];
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
    
    phone_filename = [sentence_dir, file_name, '.tsylbPHN'];
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

    pron_filename = [sentence_dir, file_name, '.WRDphn'];
    fid = fopen(pron_filename, 'w');

    fprintf(outfile_fid, '%s ', '$##');
    
    for w = 1:length(word_lengths)
        
        this_word_indicator = all(phone_indices >= word_indices(w, 1) & phone_indices <= word_indices(w, 2), 2);
        
        this_pronunciation = phones(this_word_indicator);

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

fclose(outfile_fid);

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

