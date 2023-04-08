function posterior = calc_word_posterior(word_seq, landmark)
% word_seq = (cell) sequence of words used to construct template to compare
%               to time series of observed (prior) distributions over phonemes 
%               (i.e., observed phoneme sequence).
% landmark = (struct w/ fields phase, measure) indicates current temporal
%               landmark used to align word sequence template w/ observed
%               phoneme sequence.

if nargin < 2, method = ''; end

global stats

stats = loadStats;

fields = fieldnames(likelihood);

for f = 1:length(fields)
%     if ~strcmpi(fields{f}, 'likelihood')
    eval(sprintf('%s = likelihood.%s;', fields{f}, fields{f}))
%     end    
end
% likelihood_data = likelihood;
% likelihood = likelihood.likelihood;
num_phones = length(phones);

phones2timit = cellfun(@(x) strcmpi(x, phones), stats.id, 'unif', 0);
phones2timit_mat = cat(2, phones2timit{:});

trans2timit = cellfun(@(x) strcmpi(x, phones), stats.trans_phones, 'unif', 0);
trans2timit_mat = cat(2, trans2timit{:});
trans_prob = nanunitsum(trans2timit_mat*stats.trans_prob*trans2timit_mat');

[hazard, transition, trans_posterior, inv_entropy] = deal(zeros(size(trans_likelihood)));
phone_posterior = nan(size(phone_likelihood));

phone_posterior(:, 1) = phone_likelihood(:, 1).*(phones2timit_mat*stats.prob);
trans_posterior(1) = eps;

trans_lookback = 10; % In timesteps.

last_transition = 0;
[last_transition_index, last_posterior_index] = deal(1);

for w = 2:length(time)

    duration = time(w) - last_transition;
    
    time_since_last_prior = time(w) - time(last_posterior_index);% last_transition;

    if transition(w) % time_since_certainty <= mean(diff(time))

        phone_prior = trans_prob*phone_posterior(:, last_posterior_index);

    else
        
        phone_prior = phone_posterior(:, w - 1);

    end

    switch method

        case ''

            hazard(w) = phone_hazard((phones2timit_mat')*phone_prior, duration);


        case 'exponential'

            hazard(w) = 1 - exp(-time_since_last_prior/mean(diff(time)));

        case 'zero_hazard'

            hazard(w) = 0;

    end
    
    phone_posterior(:, w) = nanunitsum(phone_likelihood(:, w).*(phone_prior*(1 - hazard(w)) + trans_prob*phone_prior*hazard(w)));

    % tp_derivative(w) = trans_likelihood(w)*
    trans_posterior(w) = trans_likelihood(w); %*hazard(w); % trans_posterior(w - 1)*hazard(w);

    inv_entropy(w) = (log(num_phones) + nansum(phone_posterior(:, w).*log(phone_posterior(:, w))))/num_phones;

    pt_dist(:, w) = nanunitsum([phone_posterior(:, w); trans_posterior(w)]);

    if trans_posterior(w) > 0.9
        
        last_transition = time(w);
        last_transition_index = w;
        transition(w) = 1;

%     else

%         last_posterior_index = w;

        if ~transition(w - 1)

            last_posterior_index = w - 1;

        end
    
    end

end

posterior = struct('likelihood', likelihood, 'phone_posterior', phone_posterior, 'trans_posterior', trans_posterior,...
    'hazard', hazard, 'transition', transition, 'pt_dist', pt_dist, 'inv_entropy', inv_entropy);

end

function hazard = phone_hazard(dist, duration)

global stats

[~, cdf_index] = cellfun(@(x) min(abs(x(:, 1) - duration)), stats.cdf, 'unif', 0);
duration_hazard = cellfun(@(x, y) x(y, 2), stats.cdf, cdf_index);

hazard = sum(duration_hazard.*nanunitsum(dist));

end

function intensity = phone_intensity(dist, duration)

global stats

pdf = cellfun(@(x) [cumsum(x(:, 1) + diff([x(:, 1); x(end, 1)])/2), stats.cdf, 'unif', 0);

[~, cdf_index] = cellfun(@(x) min(abs(x(:, 1) - duration)), stats.cdf, 'unif', 0);
duration_hazard = cellfun(@(x, y) x(y, 2), stats.cdf, cdf_index);

hazard = sum(duration_hazard.*nanunitsum(dist));

end

    

% for c = 1:num_candidates
% 
%     this_candidate = candidates{c};
%     phone_indices = cellfun(@(x) find(strcmpi(phone_list, this_candidate)));
% 
%     num_units = length(this_candidate);
% 
%     unit_index = 1;
% 
% 
%             if transition(w - 1) = 1
% 
%                 new_phone = this_candidate{unit_index + 1};
%                 new_phone_index = find(strcmpi(phone_list, new_phone));
% 
%                 phoneme_posterior(:, w) = phone_likelihood(:, w)*calc_phone_sequence_prior(this_candidate(1:(unit_index + 1));
% 
%                 phoneme_posterior(w) = phone_likelihood(phone_indices(unit_index + 1), w)*calc_phone_sequence_prior(this_candidate(1:(unit_index + 1));
% 
%                 transition_posterior(w) = phone_cdf(duration)*
% 
%             else
% 
%                 phoneme_posterior(:, w) = phone_likelihood(:, w)*phone_posterior(:, w - 1);
% 
%                 transition_posterior(w) = phone_cdf(duration)*trans_likelihood(w - 1);
% 
% 
% 
%             end
% 
%         end
% 
%     end
% 
% end
% 
% end
% 