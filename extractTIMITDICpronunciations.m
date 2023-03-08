timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

dict_file = [timit_dir, 'DOC/TIMITDIC_noheader.txt'];
dict_fid = fopen(dict_file, 'r');
dict_lines = textscan(dict_fid, '%s', 'headerLines', 1, 'Delimiter', '\n');
fclose(dict_fid);
dict_lines = dict_lines{1};

pron_file = [timit_dir, 'DOC/TIMITDIC_pronunciations.txt'];
pron_fid = fopen(pron_file, 'w');

for w = 1:length(dict_lines)

    word_pron = extractBetween(dict_lines{w}, ' /','/');
    word_pron = erase(word_pron, {'1', '2'});
    phonemes = strsplit(word_pron{:}, ' ', 'CollapseDelimiters', true);

    phonemes = tsylb_transcribe(phonemes);

    fprintf(pron_fid, '%s\n', strjoin(phonemes, ' '));

end

fclose(pron_fid);

function phonemes_out = tsylb_transcribe(phonemes);

% Closure labels suppression:
indbcl=strmatch('bcl', phonemes);
TF = strcmp('b', phonemes(indbcl+1)); % if not followed by matching plosive
phonemes(indbcl(~TF)) = cellstr('b'); % rewrite as corresponding phoneme
phonemes(indbcl(TF))=[];              % else delete.

inddcl=strmatch('dcl', phonemes);
TF = logical(strcmp('d', phonemes(inddcl+1))+strcmp('jh', phonemes(inddcl+1)));
phonemes(inddcl(~TF)) = cellstr('d');
phonemes(inddcl(TF))=[];

indgcl=strmatch('gcl', phonemes);
TF = strcmp('g', phonemes(indgcl+1));
phonemes(indgcl(~TF)) = cellstr('g');
phonemes(indgcl(TF))=[];

indpcl=strmatch('pcl', phonemes);
TF = strcmp('p', phonemes(indpcl+1));
phonemes(indpcl(~TF)) = cellstr('p');
phonemes(indpcl(TF))=[];

indtcl=strmatch('tcl', phonemes);
TF = logical(strcmp('t', phonemes(indtcl+1))+strcmp('ch', phonemes(indtcl+1)));
phonemes(indtcl(~TF)) = cellstr('t');
phonemes(indtcl(TF))=[];

indkcl=strmatch('kcl', phonemes);
TF = strcmp('k', phonemes(indkcl+1));
phonemes(indkcl(~TF)) = cellstr('k');
phonemes(indkcl(TF))=[];

% Rewrite sequence /hv w/ as /wh/
indhv=strmatch('hv', phonemes);
TF = strcmp('w', phonemes(indhv+1));
phonemes(indhv(TF)+1)=cellstr('wh');
phonemes(indhv(TF))=[];

% Replace /pau/ by #:
indpau=strmatch('pau', phonemes);
phonemes(indpau)=cellstr('#');

% % Replace /h#/ by $## or ##$:
% indh=strmatch('h#', phonemes);
% phonemes(indh(1))=cellstr('$##');
% phonemes(indh(2))=cellstr('##$');

% Replace /eng/ by /enx/
indeng=strmatch('eng', phonemes);
phonemes(indeng)=cellstr('enx');

%--------------------------------------------------------------------------
% Autres modifications dont je suis moins sure:

% Replace /ng/ by /nx/
indng=strmatch('ng', phonemes);
phonemes(indng)=cellstr('nx');

%Replace /epi/ w/ /#/. % Suppress /epi/:
indepi=strmatch('epi', phonemes);
phonemes(indepi)=[];
% phonemes(indepi)=cellstr('#');

% Replace /ax-h/ by /ah/
indaxh=strmatch('ax-h', phonemes);
phonemes(indaxh)=cellstr('ah');

% Replace /hv/ by /hh/
indhv=strmatch('hv', phonemes);
phonemes(indhv)=cellstr('hh');

phonemes_out = phonemes;

end