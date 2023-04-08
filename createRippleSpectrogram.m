function [rippleSpectrogram, t, f] = createRippleSpectrogram(varargin)
% createRippleSpectrogram - Create the spectogram of ripple stimuli.
%
% [RIPPLESPECTROGRAM, T, F] = createRippleSpectrogram();
%   Create the spectogram of a ripple stimuli RIPPLESPECTROGRAM
%   with time axis T (across columns) and frequency axis F (across
%   rows). A complete description of such spectrogram (and its stimuli)
%   is given in:
%
%       Chi, Taishih, et al. "Spectro-temporal modulation transfer
%           functions and speech intelligibility." The Journal of the
%           Acoustical Society of America 106.5 (1999): 2719-2732.
%
%       Elhilali, Mounya, Taishih Chi, and Shihab A. Shamma. "A
%           spectro-temporal modulation index (STMI) for assessment of 
%           speech intelligibility." Speech communication 41.2 (2003):
%           331-348.
%
%   The spectogram generated using the default values corresponds to
%   Fig. 4A (right panel) of elhilali2003spectro
%
%   Full credit goes to the authors of the papers.
%
%   Inputs
%       All inputs are optional. If an input is not given, it will use its
%       default value (denoted in parenthesis).
%
%       FS      Sampling rate [Hz] (44100).
%       DUR     Duration [s] (1).
%       L       Level (1).
%       DELTAA  Modulation depth [0 - 1] (1).
%       NTONES  Number of tones (280).
%       FLOW    Lowest frequency [Hz] (250).
%       FHIGH   Highest frequency [Hz] (8000).
%       PHI     Phase [rad] (0).
%       OMEGA   Ripple velocity [cyc/s] (2)
%               The number of ripple cycles per second
%               sweeping past the low-frequency edge of the spectrum.
%       OMEGA2  Ripple density [cyc/oct] (0.6).
%
%
%   Outputs
%
%       RIPPLESPECTROGRAM   Matrix. Ripple spectrogram.
%       T                   Array. Time axis [s].
%       F                   Array. Frequency axis [Hz].
%
%
% Created March 2, 2016.
% Arturo Moncada-Torres
%   arturo.moncadatorres@med.kuleuven.be
%   http://www.arturomoncadatorres.com


%% Manage (optional) inputs.

% Number of optional inputs.
nVarargin = length(varargin);

if nVarargin<0 || nVarargin>10
    error([mfilename,':inputs'],'Invalid number of input.');    
end

% Default values.
optionalArguments = {44100, ...     % Fs [Hz]
                        1, ...      % Duration [s]
                        1, ...      % Level
                        1, ...      % Modulation depth [0 - 1].
                        280, ...    % Number of tones.
                        250, ...    % Bottom frequency [Hz]
                        8000, ...   % Top frequency [Hz].
                        0, ...      % Phase [rad].
                        2, ...      % Ripple velocity [cyc/s]
                        0.6};       % Ripple density [cyc/oct].

% Overwrite the values defined in varargin.
optionalArguments(1:nVarargin) = varargin;

% Empty into appropriate variables.
[fs, dur, L, deltaA, nTones, fLow, fHigh, phi, omega, Omega2] = optionalArguments{:};


%% Important parameters.

% Number of samples.
nSamples = round(dur*fs);

% Frequency vector [Hz].
f = logspace(log10(fLow),log10(fHigh),nTones);    

% Time vector [s].
t = linspace(0,nSamples-1,nSamples) ./ fs;

% Phase [rad].
% In this case, we will consider that the phase of the ripple is
% the same for all tones. This can, of course, be changed (remember to
% keep it in rad!).
phi = ones(1,nTones) .* phi;

% Ripple velocity [cyc/s].
% For some weird reason, to match the Fig. 4A, this has to be in negative.
% Probably just a different convention while its definition.
omega = -omega;

% Tonotopic axis.
x = log2(f./fLow);


%% Create spectrogram.
rippleSpectrogram = zeros(nTones,nSamples);
for ii = 1:nTones
    for jj = 1:nSamples
        rippleSpectrogram(ii,jj) = L * (1 + deltaA .* sin(2.*pi.*(omega .* t(jj) + Omega2 .* x(ii)) + phi(ii)));
    end
end


%% Spectrogram.
hFig = figure();
hAxes = axes('Parent',hFig);
axis square;
imagesc(t, log10(f), rippleSpectrogram, 'Parent',hAxes); 
xlabel(hAxes, 'Time [s]');
ylabel(hAxes, 'Frequency [Hz]');
yTicks = [250, 500, 1000, 2000, 4000, 8000];
set(hAxes,'YTick',log10(yTicks));
set(hAxes,'YTickLabel',fliplr(yTicks));