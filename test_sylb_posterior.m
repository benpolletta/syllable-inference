cd /projectnb/crc-nak/brpp/
startupSImodel
[sentence, vowel_lh, likelihood, posterior] = plotSLPV;
% vowel_lh_trunc = truncate_struct(vowel_lh, vowel_lh.time < 2);
% sentence_trunc = truncate_struct(sentence, sentence.time < 2);
% likelihood_trunc = truncate_struct(likelihood, likelihood.time < 2);
[sylb_posterior, vocalic_nuclei] = calc_syllable_path_posterior_parallel(sentence, vowel_lh, likelihood, [], [], 10^-2);