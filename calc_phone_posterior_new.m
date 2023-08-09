function posterior = calc_phone_posterior_new(likelihood, method, tsylb_option)

if nargin < 2, method = ''; end
if nargin < 3, tsylb_option = []; end
if isempty(method), method = 'phone-specific'; end
if isempty(tsylb_option), tsylb_option = 1; end

global stats

stats = loadStats(tsylb_option);

fields = fieldnames(likelihood);

for f = 1:length(fields)
%     if ~strcmpi(fields{f}, 'likelihood')
    eval(sprintf('%s = likelihood.%s;', fields{f}, fields{f}))
%     end    
end
% likelihood_data = likelihood;
% likelihood = likelihood.likelihood;
num_phones = length(phones);

phones2timit = cellfun(@(x) strcmpi(x, phones), stats.phones.id, 'unif', 0);
phones2timit_mat = cat(2, phones2timit{:});

trans2timit = cellfun(@(x) strcmpi(x, phones), stats.phone_trans.phones, 'unif', 0);
trans2timit_mat = cat(2, trans2timit{:});
trans_prob = nanunitsum(trans2timit_mat*stats.phone_trans.prob*trans2timit_mat');

[haz_est, haz_ent, transition, trans_posterior, inv_entropy] = deal(zeros(size(trans_likelihood)));
phone_posterior = nan(size(phone_likelihood));
hazard = nan(size(phones2timit_mat, 2), size(phone_likelihood, 2));

phone_posterior(:, 1) = phone_likelihood(:, 1).*(phones2timit_mat*stats.phones.prob);
trans_posterior(1) = eps;

trans_lookback = 10; % In timesteps.
% decay_time = 1000; % In timesteps.

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

%     unif_phone_dist = nanunitsum(ones(size(phone_prior)));
%     phone_prior = phone_prior + (unif_phone_dist - phone_prior)/decay_time;

    switch method

        case 'phone-specific'

            this_dist = nanunitsum((phones2timit_mat')*phone_prior);
            hazard(:, w) = phone_hazard(duration);
            haz_est(w) = sum(hazard(:, w).*this_dist);
            % haz_var(w) = sum((hazard(:, w).^2).*this_dist) - haz_est(w)^2;
            haz_ent(w) = -nansum(nanunitsum(hazard(:, w)).*log(nanunitsum(hazard(:, w))))/log(num_phones); % (log(num_phones) + nansum(hazard(:, w).*log(hazard(:, w))))/num_phones;

        case 'non-phone-specific'

            this_dist = (phones2timit_mat')*phone_prior;
            haz_est(w) = sum(phone_hazard(duration).*nanunitsum(this_dist));
            haz_ent(w) = 0;
            hazard(:, w) = haz_est(w)*ones(size(this_dist));


        case 'exponential'

            haz_est(w) = 1 - exp(-time_since_last_prior/mean(diff(time)));
            haz_ent(w) = 0;
            hazard(:, w) = haz_est(w)*ones(size(phone_prior));

        case 'zero_hazard'

            haz_est(w) = 0;
            haz_ent(w) = 0;
            hazard(w) = zeros(size(phone_prior));

    end
    
    phone_posterior(:, w) = nanunitsum(phone_likelihood(:, w).*(phone_prior.*(1 - phones2timit_mat*hazard(:, w)) + trans_prob*phone_prior.*phones2timit_mat*hazard(:, w)));

    % tp_derivative(w) = trans_likelihood(w)*
    trans_posterior(w) = trans_likelihood(w) + haz_ent(w)*haz_est(w); %*hazard(w); % trans_posterior(w - 1)*hazard(w);

    inv_entropy(w) = (log(num_phones) + nansum(phone_posterior(:, w).*log(phone_posterior(:, w))))/num_phones;

    pt_dist(:, w) = nanunitsum([phone_posterior(:, w); trans_posterior(w)], [], 'uniform');

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
    'hazard', hazard, 'hazard_est', haz_est, 'hazard_var', haz_ent, 'transition', transition, 'pt_dist', pt_dist, 'inv_entropy', inv_entropy);

end

function hazard = phone_hazard(duration) % (dist, duration)

global stats

[~, cdf_index] = cellfun(@(x) min(abs(x(:, 1) - duration)), stats.phones.cdf, 'unif', 0);
hazard = cellfun(@(x, y) x(y, 2), stats.phones.cdf, cdf_index); 
% duration_hazard = cellfun(@(x, y) x(y, 2), stats.cdf, cdf_index);

% hazard = sum(duration_hazard.*nanunitsum(dist));

end

function intensity = phone_intensity(dist, duration)

global stats

pdf = cellfun(@(x) [cumsum(x(:, 1) + diff([x(:, 1); x(end, 1)])/2), diff]);

[~, cdf_index] = cellfun(@(x) min(abs(x(:, 1) - duration)), stats.phones.cdf, 'unif', 0);
duration_hazard = cellfun(@(x, y) x(y, 2), stats.phones.cdf, cdf_index);

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