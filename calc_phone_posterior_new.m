function posterior = calc_phone_posterior_new(likelihood)

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

phones2timit = cellfun(@(x) strcmpi(x, phones), stats.id, 'unif', 0);
phones2timit_mat = cat(2, phones2timit{:});

trans2timit = cellfun(@(x) strcmpi(x, phones), stats.trans_phones, 'unif', 0);
trans2timit_mat = cat(2, trans2timit{:});
trans_prob = trans2timit_mat*stats.trans_prob*trans2timit_mat';

[hazard, transition, trans_posterior] = deal(zeros(size(trans_likelihood)));
phone_posterior = nan(size(phone_likelihood));

phone_posterior(:, 1) = phone_likelihood(:, 1).*(phones2timit_mat*stats.prob);
trans_posterior(1) = eps;

trans_lookback = 10; % In timesteps.

last_transition = 0;
[last_transition_index, last_posterior_index] = deal(1);

for w = 2:length(time)

    duration = time(w) - last_transition;

    if transition(w) % duration <= mean(diff(time))

        phone_prior = trans_prob*phone_posterior(:, last_posterior_index);

    else
        
        phone_prior = phone_posterior(:, w - 1);

    end

    hazard(w) = phone_hazard((phones2timit_mat')*phone_prior, duration);

    phone_posterior(:, w) = nanunitsum(phone_likelihood(:, w).*(phone_prior*(1 - hazard(w)) + trans_prob*phone_prior*hazard(w)));
    % tp_derivative(w) = trans_likelihood(w)*
    trans_posterior(w) = trans_likelihood(w); %*hazard(w); % trans_posterior(w - 1)*hazard(w);

    pt_dist(:, w) = nanunitsum([phone_posterior(:, w); trans_posterior(w)]);

    if trans_posterior(w) > 0.9
        
        last_transition = time(w);
        last_transition_index = w;
        transition(w) = 1;

        if transition(w - 1) = 0;

            last_posterior_index = w - 1;

        end
    
    end

end

posterior = struct('phone_posterior', phone_posterior, 'trans_posterior', trans_posterior,...
    'hazard', hazard, 'transition', transition, 'pt_dist', pt_dist);

end

function hazard = phone_hazard(dist, duration)

global stats

[~, cdf_index] = cellfun(@(x) min(abs(x(:, 1) - duration)), stats.cdf, 'unif', 0);
duration_hazard = cellfun(@(x, y) x(y, 2), stats.cdf, cdf_index);

hazard = sum(duration_hazard.*nanunitsum(dist));

end

function intensity = phone_intensity(dist, duration)

global stats

pdf = cellfun(@(x) [cumsum(x(:, 1) + diff([x(:, 1); x(end, 1)])/2), diff])

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