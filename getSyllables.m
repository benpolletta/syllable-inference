timit_dir = 'timit/TIMIT/';
tsylb_dir = 'tsylb2-1.1/';

timit_filenames_file = [timit_dir, 'DOC/allphonelist_filenames.txt'];
filename_fid = fopen(timit_filenames_file, 'r');
timit_filenames = textscan(filename_fid, '%s', 'Delimiter', '\n');
fclose(filename_fid);
timit_filenames = timit_filenames{1};

tsylb_output_file = [timit_dir, 'DOC/allphonelist_tsylb_out.txt'];
tsylb_fid = fopen(tsylb_output_file, 'r');
tsylb_output = textscan(tsylb_fid, '%s', 'headerLines', 1, 'Delimiter', '\n');
fclose(tsylb_fid);
tsylb_output = tsylb_output{1};

output_line_indicator = contains(tsylb_output, 'Enter');
tsylb_output(~output_line_indicator) = [];
tsylb_output(end) = [];
%tsylb_output(cellfun(@isempty, tsylb_output)) = [];
tsylb_output = cellfun(@(x) extractAfter(x, 'Basic pron is'), tsylb_output, 'unif', 0); % extractBetween(x, '/#', '#/'), 

for s = 1:length(tsylb_output)
    
    sentence_filename = [timit_dir, timit_filenames{s}, '.tsylbPHN'];
    fid = fopen(sentence_filename, 'r');
    sentence_phones = textscan(fid, '%s');
    fclose(fid);
    sentence_phones = sentence_phones{1};
    sentence_phones = reshape(sentence_phones, 3, length(sentence_phones)/3);
   
    sentence_tsylb = tsylb_output{s};
    
    sentence_split = strsplit(sentence_tsylb, {'[','] '}, 'CollapseDelimiters', true);
    
    sylb_split = cellfun(@(x) textscan(x, '%s', 'delimiter', {' ', '''0'}), sentence_split, 'unif', false); % strsplit(x, {' ','''0'}, 'CollapseDelimiters', true), sentence_split, 'unif', 0);
    
    sylb_split = cellfun(@(x) x{:}, sylb_split, 'unif', 0);
    
    phones_per_sylb = cellfun(@length, sylb_split) - cellfun(@(x) sum(cellfun(@isempty, x)), sylb_split);
    
    boundary_phone_numbers = cumsum(phones_per_sylb);
    
    boundary_indices = cellfun(@str2num, sentence_phones(2, boundary_phone_numbers(1:(end - 1))));
    
    boundary_filename = [timit_dir, timit_filenames{s}, '.TSYLB'];
    fid = fopen(boundary_filename, 'w');
    fprintf(fid, '%d\n', boundary_indices');
    fclose(fid);
    
end