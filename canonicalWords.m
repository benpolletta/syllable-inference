function results = canonicalWords(SI)
% Function used to disambiguate words found in the corpus but not in the
% TIMIT lexicon via user input. Stores and reuses particular ambiguous
% sentence/word pairs to reduce amount of user input necessary to
% disambiguate all words in the corpus.

global onset_time

if nargin < 1, SI = []; end
if isempty(SI), SI = (1:6300)/6300; end
if nargin < 2, reuse_results = []; end
if isempty(reuse_results), reuse_results = true; end

name = 'Data';
if length(SI) ~= 6300
    name = sprintf('%s_%dsentences', name, length(SI));
end 

timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

timit_filenames_file = [timit_dir, 'DOC/allphonelist_filenames.txt'];
filename_fid = fopen(timit_filenames_file, 'r');
file_list = textscan(filename_fid, '%s', 'Delimiter', '\n');
fclose(filename_fid);
file_list = file_list{1};

wsp_map = load('word2sylb2phone_bysentence.mat');
wsp_map = wsp_map.results;

dict = load('word2sylb2phone_DICT.mat');
dict = dict.results;
dict_words = {dict(:).word};
% dict_words = cellfun(@(x) erase(x, '-'), dict_words, 'unif', 0);

[ambiguous_words, ambiguous_sentences, resolutions] = deal({});

for s = 1:length(SI)

    file_index = SI(s);

    if isfloat(file_index) && file_index < 1 && file_index > 0
        file_index = round(file_index*length(file_list));
    end

    file_name = file_list{file_index};

    %% Retrieving words, phonemes, and syllables and their start and end times.

    word_filename = [timit_dir, file_name, '.WRD'];
    words = getUnits(word_filename);

    %% Finding canonical words.

    canonical_words = cell(size(words));

    for w = 1:length(words)

        if any(strcmp(words{w}, dict_words))

            canonical_words{w} = dict_words(strcmp(words{w}, dict_words));

        else

            possible_words = dict_words(contains(dict_words, words{w}));

            if length(possible_words) == 1

                canonical_words(w) = possible_words(:);

            else

                sentence = strjoin(words, ' ');

                amb_index = find(strcmp(words{w}, ambiguous_words) & strcmp(sentence, ambiguous_sentences));

                if ~isempty(amb_index)

                    canonical_words(w) = resolutions(amb_index);

                else

                    ambiguous_words(end + 1) = words(w);

                    ambiguous_sentences(end + 1) = {sentence};

                    prompt = sprintf('No unique match for %s in sentence\n ''%s'' \n is found in TIMITDIC. Possible matches are:\n %s.\n Please enter the index of the correct word:\n',words{w}, strjoin(words, ' '), strjoin(possible_words, ', '));

                    index = input(prompt);

                    [canonical_words(w), resolutions(end + 1)] = deal(possible_words(index));

                end

            end

        end

    end

    %% Saving results.

    results(s) = struct('words', {words}, 'canonical_words', {canonical_words}); % 'syllabifications', {syllabifications});

end

save('canonicalWords.mat', 'results', 'ambiguous_words', 'ambiguous_sentences', 'resolutions') %, 'wsp_map')

end

function [units, times, lengths] = getUnits(filename)

global onset_time

fid = fopen(filename, 'r');
data = textscan(fid, '%d%d%s'); % '%s');
fclose(fid);
% data = data{1};
% data = reshape(data, 3, length(data)/3);

units = data{3}; % (3, :);

indices = cell2mat(data(1:2)); % cellfun(@str2num, data(1:2, :))';
times = (indices/16 + onset_time);

lengths = diff(times, [], 2);

end