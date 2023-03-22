function results = writeDICsyllables

timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

dic_file = [timit_dir, 'DOC/TIMITDIC_noheader.txt'];
dic_fid = fopen(dic_file, 'r');
dic_words = textscan(dic_fid, '%s', 'headerLines', 1, 'Delimiter', '\n');
fclose(dic_fid);
dic_words = dic_words{1};
dic_words = extractBefore(dic_words, '/');

tsylb_output_file = [timit_dir, 'DOC/TIMITDIC_tsylb_out.txt'];
tsylb_fid = fopen(tsylb_output_file, 'r');
tsylb_output = textscan(tsylb_fid, '%s', 'Delimiter', '\n');
fclose(tsylb_fid);
tsylb_output = tsylb_output{1};

output_line_indicator = contains(tsylb_output, 'Enter');
tsylb_output(~output_line_indicator) = [];
tsylb_output(end) = [];
%tsylb_output(cellfun(@isempty, tsylb_output)) = [];
tsylb_output = cellfun(@(x) extractBetween(x, 'Basic pron is /#', '#/'), tsylb_output, 'unif', 0);

for w = 1:length(tsylb_output)
    
    %% Selecting word & pronunciation.
    word = strip(dic_words{w});
    tsylb = strip(tsylb_output{w}{1});
    
    %% Splitting pronunciation into syllables.
    sylb_cell = strsplit(tsylb, {' [','[', '] ', ']'}, 'CollapseDelimiters', true);
    sylb_nonzero_indicator = cellfun(@(x) length(x) ~= 0, sylb_cell);
    sylb_cell = sylb_cell(sylb_nonzero_indicator);

    %% Splitting syllables into phones.
    sylb_phone_cell = cellfun(@(x) strsplit(x, {' ','''0'}, 'CollapseDelimiters', true), sylb_cell, 'unif', 0);
    sylb_phone_nonzero_indicator = cellfun(@(x) cellfun(@(y) length(y) ~= 0, x), sylb_phone_cell, 'unif', 0);
    sylb_phone_cell = cellfun(@(x, y) x(y), sylb_phone_cell, sylb_phone_nonzero_indicator, 'unif', 0);

    %% Rejoining phones.
    sylb_cell = cellfun(@(x) strjoin(x, '/'), sylb_phone_cell, 'unif', 0);
    pron = strjoin(sylb_cell, '*');

    %% Recounting syllables & phones per word.
    sylb_phone_num = cellfun(@length, sylb_phone_cell, 'unif', 0);
    phone_num = cellfun(@sum, sylb_phone_num);
    sylb_num = length(sylb_cell);

    %% Saving results.
    results(w) = struct('word', word, 'pron', pron, 'sylb_cell', {sylb_cell}, 'sylb_phone_cell', {sylb_phone_cell},...
        'sylb_phone_num', {sylb_phone_num}, 'phone_num', phone_num, 'sylb_num', sylb_num);
    
end

save('word2sylb2phone_DICT.mat', 'results')

end
