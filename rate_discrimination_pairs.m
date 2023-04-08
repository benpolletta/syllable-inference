function [tone_times, tones, time] = rate_discrimination_pairs(standard_rate, rel_diff)

if nargin < 2, rel_diff = []; end
if isempty(rel_diff), rel_diff = .4*rand - .2; end
if nargin < 1, standard_rate = []; end
if isempty(standard_rate)

    standard_rates = [4, 5.57, 7.14, 8.71, 10.29, 11.86, 13.43, 15];

    standard_rate = standard_rates(ceil(rand*length(standard_rates)));

end

this_comparison = standard_rate*(1 + rel_diff);

n_rand = rand;
n1 = 5*(n_rand < .5) + 7*(n_rand >= .5);
n2 = 5*(n_rand >= .5) + 7*(n_rand < .5);

standard_interval = 1/standard_rate;
comparison_interval = 1/this_comparison;

inter_sequence_interval = (rand + 4.5)*standard_interval;
duration = ceil(n1*(standard_interval)+2*inter_sequence_interval+n2*comparison_interval);

time = 0:.001:duration;

ITIs = [inter_sequence_interval/2, ones(1, n1)*standard_interval, inter_sequence_interval, ones(1,n2)*comparison_interval];
tone_times = cumsum(ITIs);
tones = zeros(size(time));

for t = 1:length(tone_times)
    
    tones(abs(time - tone_times(t)) < 7.5*10^-3) = 1;
end

end