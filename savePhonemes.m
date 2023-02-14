timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';
tsylb_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/tsylb2-1.1/';

timit_phonemes_file = [timit_dir, 'DOC/allphonelist.txt'];
fid = fopen(timit_phonemes_file, 'r');
timit_phonemes = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid)
timit_phonemes = timit_phonemes{1};

tsylb_phonemes_file = [timit_dir, 'DOC/allphonelist_tsylb.txt'];
sentence_fid = fopen(tsylb_phonemes_file, 'w');

timit_filenames_file = [timit_dir, 'DOC/allphonelist_filenames.txt'];
filename_fid = fopen(timit_filenames_file, 'w');


for s = 1:length(timit_phonemes)
   
    sentence = timit_phonemes{s};
    
    sentence_split = strsplit(sentence);
    
    filename = upper(sentence_split{1});
    
    fprintf(filename_fid, '%s\n', filename);
    
    % sentence_phones = sentence_split(2:(end - 1));
    
    sentence_phones = timit_transcription_brpp( filename ); % sentence_phones = cellfun(@tsylb_substitutions, sentence_phones, 'unif', 0);
    
    sentence = join(sentence_phones, ' ');
    
    fprintf(sentence_fid, '%s\n', sentence{:});
    
end

fclose(sentence_fid)
fclose(filename_fid)

system(sprintf('%s./tsylb2 -n phon1ax.pcd < %s > %sDOC/allphonelist_tsylb_out.txt', tsylb_dir, tsylb_phonemes_file, timit_dir))

function phone = tsylb_substitutions(phone)

phones_to_sub = {'cl', 'h#', 'hv', 'epi', 'pau', 'ng', 'ax-h';...
    'q', '#', 'hh', '#', '#', 'nx', 'ax'};

for p = 1:length(phones_to_sub)
    
    if contains(phone, phones_to_sub{1, p})
        
        phone = phones_to_sub{2, p};
        
    end
    
end
    
end