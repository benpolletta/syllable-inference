
global stats vowels

tsylb_option = 1;
stats = loadStats(tsylb_option);
prons = stats.prons.id;
prons{end + 1} = '#';
counts = stats.prons.count;
counts(end + 1) = 2*6300;
probs = counts/sum(counts); % stats.prons.prob;
trans_matrix = stats.pron_trans.prob;
trans_matrix = trans_matrix*diag(1./probs);

[tsylb_phonemes, class_indicator, class_names] = getPhones(1);
vowels = tsylb_phonemes(class_indicator(:, strcmpi(class_names, 'vowels')));
vowels = {vowels{:}, 'el', 'em', 'en', 'enx'};

five_prons = prons(round(rand(5,1)*length(prons)));