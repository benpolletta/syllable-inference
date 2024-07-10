function sentence = generateSentence(num_words)

stats = loadStats(1);

meanrate = stats.sylb_dur.stats(1);
stdrate = stats.sylb_dur.stats(2);
rate = normrnd(meanrate, stdrate);

word_mc = dtmc(stats.word_trans.t_prob);

time = 0;
w0 = randsample(total_words, 1, true, word_trans(1, :));
s0 = randsample(total_sylbs, 1, true, stats.word2trans(w0));
p0 = randsample(total_sylbs, 1, true, stats.sylb2trans(s0));
f0 = zeros(1,num_features);

i = 1;

while n_word <= num_words

    time = time + dt;

    rate_prod(i) = rate; %*word_rate*sylb_rate*phone_rate;

    phoneme_time(i + 1) = phoneme_time(i) + rate_prod*dt; % phoneme_dur = (1/rate_prod - time + onset)/rate_prod;

    phone_transition = 0;
    if phone_time(i+1) >= 1
        phone_transition = 1;
        phone_time(i+1) = 0;
    end

    phoneme_state(i + 1) = (1 - phoneme_transition)*phoneme_state(i) + phoneme_transition*sylb2trans(sylb_state(i))*phoneme_state(i);

    spect(i+1, :) = phone2spect(phoneme_time(i+1), phoneme_state(i+1));

    sylb_time(i + 1) = sylb_time(i) + phoneme_transition/num_phones(sylb_state(i));

    sylb_transition = 0;
    if sylb_time(i+1) >= 1
        sylb_transition = 1;
        sylb_time(i+1) = 0;
    end

    sylb_state(i + 1) = (1 - sylb_transition)*sylb_state(i) + sylb_transition*word2trans(word_state(i))*sylb_state(i);

    word_time(i + 1) = word_time(i) + sylb_transition/num_phones(word_state(i));

    word_transition = 0;
    if word_time(i+1) >= 1
        word_transition = 1;
        word_time(i+1) = 0;
    end

    word_state(i + 1) = (1 - word_transition)*word_state(i) + word_transition*word2trans(word_state(i))*word_state(i);

end

sentence = {word_state, sylb_state, phoneme_state};

% words = simulate(word_mc, num_words);
% 
% for w = 1:num_words
% 
%     word_dur_factor = normrnd(0, 1); %word_duration = word2numsylbs(words(w))/rate;
% 
%     sylb_mc = dtmc(stats.word2trans(words(w)));
%     sylbs = simulate(sylb_mc, 10, 'X0', 1);
%     num_sylbs_remaining = length(sylbs);
% 
%     for s = 1:num_sylbs
% 
%         sylb_dur_factor = normrnd(0, 1);
% 
%         phone_mc = dtmc(stats.sylb2trans(sylbs(s)));
%         phones = simulate(phone_mc, 10, 'X0', 1);
%         num_phones_remaining = length(phones);
% 
%         for p = 1:num_phones
% 
%             phone_dur_factor = normrnd(0, 1);
% 
%             phone_tau = (time - phone_onset)/(rate*word_dur_factor*sylb_dur_factor*phone_dur_factor);
% 
% 
% 
%         end
% 
%     end
% 
% end

end