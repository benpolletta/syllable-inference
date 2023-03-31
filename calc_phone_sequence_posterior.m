function posterior = calc_phone_sequence_posterior(likelihood, phone_sequence, method, tsylb_option)

if nargin < 3, method = ''; end
if nargin < 4, tsylb_option = []; end
if isempty(method), method = 'phone-specific'; end
if isempty(tsylb_option), tsylb_option = 1; end

global stats

stats = loadStats(tsylb_option);

fields = fieldnames(likelihood);

for f = 1:length(fields)
    eval(sprintf('%s = likelihood.%s;', fields{f}, fields{f}))
end
num_phones = length(phones);
dt = mean(diff(time));

phones2timit = cellfun(@(x) strcmpi(x, phones), stats.id, 'unif', 0);
phones2timit_mat = cat(2, phones2timit{:});

trans2timit = cellfun(@(x) strcmpi(x, phones), stats.trans_phones, 'unif', 0);
trans2timit_mat = cat(2, trans2timit{:});
trans_prob = nanunitsum(trans2timit_mat*stats.trans_prob*trans2timit_mat');

evidence_weights = zeros(length(phone_sequence), length(time));
previous_weight_shape = 1 - cdfs{1}(:, 2);
evidence_weights(1, 1:length(previous_weight_shape)) = previous_weight_shape;

for p = 2:length(phone_sequence)
    
    [durations(p), cdf] = phone_duration(phone_sequence{p});
    cdfs{p} = resample_cdf(cdf, dt);

    weight_shape = conv(1-cdfs{p}(:, 2), previous_duration_probability, 'same');
    weight_start = find(evidence_weights(p - 1, :) <= .5, 1);
    weight_loc = 1:length(time) >= weight_start & 1:length(time) <= weight_start + length(weight_shape) - 1;

    max_index = min(weight_start + length(weight_shape) - 1, size(evidence_weights, 2));
    evidence_weights(p, weight_loc) = weight_shape(1:max_index);
end

[haz_est, transition, trans_posterior, inv_entropy] = deal(zeros(size(trans_likelihood)));
phone_posterior = nan(size(phone_likelihood));
hazard = nan(size(phones2timit_mat, 2), size(phone_likelihood, 2));

phone_posterior(:, 1) = phone_likelihood(:, 1).*(phones2timit_mat*stats.prob);
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
            haz_var(w) = sum((hazard(:, w).^2).*this_dist) - haz_est(w)^2;

        case 'non-phone-specific'

            this_dist = (phones2timit_mat')*phone_prior;
            haz_est(w) = sum(phone_hazard(duration).*nanunitsum(this_dist));
            haz_var(w) = 0;
            hazard(:, w) = haz_est(w)*ones(size(this_dist));


        case 'exponential'

            haz_est(w) = 1 - exp(-time_since_last_prior/mean(diff(time)));
            haz_var(w) = 0;
            hazard(:, w) = haz_est(w)*ones(size(phone_prior));

        case 'zero_hazard'

            haz_est(w) = 0;
            haz_var(w) = 0;
            hazard(w) = zeros(size(phone_prior));

    end
    
    phone_posterior(:, w) = nanunitsum(phone_likelihood(:, w).*(phone_prior.*(1 - phones2timit_mat*hazard(:, w)) + trans_prob*phone_prior.*phones2timit_mat*hazard(:, w)));

    % tp_derivative(w) = trans_likelihood(w)*
    trans_posterior(w) = trans_likelihood(w) + haz_est(w); %*hazard(w); % trans_posterior(w - 1)*hazard(w);

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
    'hazard', hazard, 'hazard_est', haz_est, 'transition', transition, 'pt_dist', pt_dist, 'inv_entropy', inv_entropy);

end

function hazard = phone_hazard(duration) % (dist, duration)

global stats

[~, cdf_index] = cellfun(@(x) min(abs(x(:, 1) - duration)), stats.cdf, 'unif', 0);
hazard = cellfun(@(x, y) x(y, 2), stats.cdf, cdf_index); 
% duration_hazard = cellfun(@(x, y) x(y, 2), stats.cdf, cdf_index);

% hazard = sum(duration_hazard.*nanunitsum(dist));

end

function [duration, cdf] = phone_duration(phone)

global stats

phone_index = find(strcmpi(phone, stats.id));
cdf = stats.cdf{phone_index};

[~, duration_index] = min(abs(cdf(:, 2) - .5));
duration = cdf(duration_index, 1);

end

function cdf_out = resample_cdf(cdf_in, dt)

[dt_val, dt_time] = resample(cdf_in(:, 2), cdf_in(:, 1), 1/dt);

first_over_1 = find(dt_val > 1, 1);
bad_end_indicator = 1:length(dt_time) >= first_over_1;
y_end = linspace(dt_val(first_over_1 - 1), 1, sum(bad_end_indicator) + 1);

dt_val(bad_end_indicator) = y_end(2:end);

cdf_out = [dt_time, dt_val];

end

function intensity = phone_intensity(dist, duration)

global stats

pdf = cellfun(@(x) [cumsum(x(:, 1) + diff([x(:, 1); x(end, 1)])/2), diff]);

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