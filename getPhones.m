function [phonemes, class_indicator, class_names] = getPhones(tsylb_option)

if nargin == 0, tsylb_option = []; end
if isempty(tsylb_option), tsylb_option = 0; end

phone_prefix = 'TIMIT';

%% Getting lists of phonemes (and phoneme classes) used in TIMIT.

timit_dir = '/projectnb/crc-nak/brpp/Speech_Stimuli/timit/TIMIT/';

timit_phonemes_file = [timit_dir, 'DOC/TIMITphones.txt'];
fid = fopen(timit_phonemes_file, 'r');
phonemes = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
phonemes = phonemes{1};

class_borders = [find(contains(phonemes, '%%')); length(phonemes) + 1];

class_indicator = false(length(phonemes), length(class_borders) - 1);

for b = 1:(length(class_borders) - 1)
    
    class_indicator((class_borders(b) + 1):(class_borders(b + 1) - 1), b) = true;
    
end

class_names = phonemes(class_borders(1:(end - 1)));

class_names = cellfun(@(x) x(4:end), class_names, 'unif', 0);

class_indicator(class_borders(1:(end - 1)), :) = [];

phonemes(class_borders(1:(end - 1))) = [];

phonemes = cellfun(@(x) textscan(x, '%s'), phonemes, 'unif', 0);
phonemes = cellfun(@(x) x{:}, phonemes, 'unif', 0);
phonemes = cellfun(@(x) x{:}, phonemes, 'unif', 0);

%% Subsitutions for tsylb.

if tsylb_option

    phone_prefix = 'TSYLB';

    % Closure labels suppression:
    stops = {'b', 'd', 'g', 'p', 't', 'k'};

    for s = 1:length(stops)

        indcl=strmatch([stops{s}, 'cl'], phonemes);
        phonemes(indcl) = [];
        class_indicator(indcl, :) = [];

    end

    % Replace /pau/ by #:
    indpau=strmatch('pau', phonemes);
    phonemes(indpau)=cellstr('#');

    % Replace /h#/ by $## or ##$:
    indh=strmatch('h#', phonemes);
    phonemes(indh)=cellstr('$##');
    phonemes((indh + 2):(end + 1)) = phonemes((indh + 1):end);
    class_indicator((indh + 2):(end + 1), :) = class_indicator((indh + 1):end, :);
    phonemes(indh + 1)=cellstr('##$');
    class_indicator(indh + 1, :) = class_indicator(indh, :);

    % Replace /eng/ by /enx/
    indeng=strmatch('eng', phonemes);
    phonemes(indeng)=cellstr('enx');

    % Replace /ng/ by /nx/, which already exists in th elist of phonemes.
    indng=strmatch('ng', phonemes);
    phonemes(indng)=[];
    class_indicator(indng, :) = [];
    % timit_phonemes(indng)=cellstr('nx');

    % Suppress /epi/:
    indepi=strmatch('epi', phonemes);
    phonemes(indepi)=[];
    class_indicator(indepi, :) = [];
    
    %     % Replace /epi/ by #:
    %     indepi=strmatch('epi', timit_phonemes);
    %     timit_phonemes(indepi)=cellstr('#');

    % Replace /ax-h/ by /ah/, which already exists in the list of phonemes.
    indaxh=strmatch('ax-h', phonemes);
    phonemes(indaxh)=[];
    class_indicator(indaxh, :) = [];

    % Replace /hv/ by /hh/, which already exists in the list of phonemes.
    indhv=strmatch('hv', phonemes);
    phonemes(indhv)=[];
    class_indicator(indhv, :) = [];

    %     % Add /wh/
    %     indw = strmatch('w', timit_phonemes);
    %     timit_phonemes((indw + 2):(end + 1)) = timit_phonemes((indw + 1):end);
    %     class_indicator((indw + 2):(end + 1), :) = class_indicator((indw + 1):end, :);
    %     timit_phonemes(indw + 1)=cellstr('wh');
    %     class_indicator(indw + 1, :) = class_indicator(indw, :);
    
end

save([phone_prefix, 'phones.mat'], 'phonemes', 'class_names', 'class_indicator')
