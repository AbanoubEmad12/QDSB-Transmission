clc; clear; close all;
%% Step 1: Define Parameters
fs = 100e3;         %  sampling frequency 
fc = 18e3;        % Carrier frequency in Hz 
fs_audio = 44100; % Standard audio sampling rate 
filter_order = 51; % Filter order 
filter_cutoff = 10e3; %filter cutoff frequency 
%% Step 2: Load and Preprocess Audio File
[audio_data, fs_original] = audioread('/Users/abanoubkamel/Downloads/Audio Files/Amr_Diab.mp3');

% Convert to Mono if Stereo
if size(audio_data, 2) > 1
    audio_data = mean(audio_data, 2);
end
% Normalize the Audio Signal
audio_data = audio_data / max(abs(audio_data)); % Normalize between -1 and 1
% Resample to Match Modulation Sampling Rate
audio_resampled = resample(audio_data, fs, fs_original);
t = (0:length(audio_resampled)-1) / fs; 
%% Step 3: Modulation (AM DSB-SC)
carrier = cos(2 * pi * fc * t(:)); 
s_modulated = audio_resampled .* carrier; 
%% Step 4: Demodulation
demod_signal = s_modulated .* carrier; 

% Low-Pass Filter Design (Updated with filter order and cutoff frequency)
% lp_filter = fir1(filter_order, filter_cutoff/(fs/2)); % FIR filter design
% demod_filtered = 2 * filter(lp_filter, 1, demod_signal); % Apply low-pass filter
[b, a] = butter(1, filter_cutoff/(fs/2), 'low');
demod_filtered = 2 * filtfilt(b, a, demod_signal);
% [b, a] = butter(2, filter_cutoff/(fs/2), 's'); % Analog Butterworth filter design
% [num, den] = bilinear(b, a, fs); % Convert to digital equivalent using bilinear transform
% demod_filtered = 2 * filter(num, den, demod_signal); % Apply filter
% [sos, g] = butter(4, filter_cutoff/(fs/2), 'low', 's');
% demod_filtered = 2 * filtfilt(sos, g, demod_signal);
% Resample Demodulated Signal to Original Audio Rate
demod_resampled = resample(demod_filtered, fs_audio, fs);
modulated_playback = resample(s_modulated, fs_audio, fs);
modulated_playback = modulated_playback / max(abs(modulated_playback)); % Normalize
%% Step 5: Saving Audio Signals
[filename, pathname] = uiputfile({'*.wav', 'WAV files (*.wav)'}, 'Save Original Audio');
if ischar(filename)
    audiowrite(fullfile(pathname, filename), audio_data, fs_original);
    disp('Original Audio saved.');
else
    disp('Original Audio not saved.');
end
[filename, pathname] = uiputfile({'*.wav', 'WAV files (*.wav)'}, 'Save Modulated Audio');
if ischar(filename)
    audiowrite(fullfile(pathname, filename), modulated_playback, fs_audio);
    disp('Modulated Audio saved.');
else
    disp('Modulated Audio not saved.');
end
[filename, pathname] = uiputfile({'*.wav', 'WAV files (*.wav)'}, 'Save Demodulated Audio');
if ischar(filename)
    audiowrite(fullfile(pathname, filename), demod_resampled, fs_audio);
    disp('Demodulated Audio saved.');
else
    disp('Demodulated Audio not saved.');
end

%% Step 6: Time-Domain Plots

% Define Time Vectors
t_original = (0:length(audio_data)-1) / fs_original; % Time for original audio
t_resampled = (0:length(audio_resampled)-1) / fs;   % Time for resampled audio
t_modulated = (0:length(s_modulated)-1) / fs;       % Time for modulated signal
t_demodulated = (0:length(demod_filtered)-1) / fs;  % Time for demodulated signal

figure;
% Plot Original Message Signal (Time Domain)
subplot(4,1,1); plot(t_original, audio_data, 'b'); title('Original Message Signal (Time Domain)'); xlabel('Time (s)'); ylabel('Amplitude');% xlim([0 max(t_original)]);
subplot(4,1,2); time_zoom = t_modulated(1:500); carrier_zoom = carrier(1:500); 
plot(time_zoom, carrier_zoom, 'k'); title('Carrier Signal (Time Domain)'); xlabel('Time (s)'); ylabel('Amplitude');% xlim([min(time_zoom) max(time_zoom)]);
subplot(4,1,3);plot(t_modulated, s_modulated, 'r');title('Modulated Signal (AM DSB-SC)');xlabel('Time (s)'); ylabel('Amplitude');% xlim([0 max(t_modulated)]);
subplot(4,1,4);plot(t_demodulated, demod_filtered, 'g');title('Demodulated Signal');xlabel('Time (s)'); ylabel('Amplitude');% xlim([0 max(t_demodulated)]);

%% Step 7: Frequency-Domain Analysis (Magnitude Spectra)

% Compute FFT for Message Signal
N_message = 2^nextpow2(length(audio_resampled)); 
T_message = 1/fs; 
f_message = (-N_message/2:N_message/2-1) * (fs/N_message); 
message_spectrum = fftshift(abs(fft(audio_resampled, N_message))); 

% Compute FFT for Modulated Signal
N_modulated = 2^nextpow2(length(s_modulated));
T_modulated = 1/fs;
f_modulated = (-N_modulated/2:N_modulated/2-1) * (fs/N_modulated);
modulated_spectrum = fftshift(abs(fft(s_modulated, N_modulated)));

% Compute FFT for Demodulated Signal
N_demod = 2^nextpow2(length(demod_filtered));
T_demod = 1/fs;
f_demod = (-N_demod/2:N_demod/2-1) * (fs/N_demod);
demod_spectrum = fftshift(abs(fft(demod_filtered, N_demod)));

%% step 8: Frequency-Domain Plots (Magnitude Spectra)
figure;
subplot(3,1,1);plot(f_message/1e3, message_spectrum, 'b'); title('Frequency Spectrum of Message Signal'); xlabel('Frequency (kHz)'); ylabel('Magnitude');% xlim([-5 5]); % Focus on baseband frequencies
subplot(3,1,2); plot(f_modulated/1e3, modulated_spectrum, 'r'); title('Frequency Spectrum of Modulated Signal'); xlabel('Frequency (kHz)'); ylabel('Magnitude');% xlim([-fc/1e3 - 10, fc/1e3 + 10]); % Centered around carrier frequency
subplot(3,1,3); plot(f_demod/1e3, demod_spectrum, 'm'); title('Frequency Spectrum of Demodulated Signal'); xlabel('Frequency (kHz)');ylabel('Magnitude');% xlim([-5 5]); % Focus on baseband frequencies
grid on;

%% Step 9: Revised PSD Comparison Plot

% Compute PSD for Original Message Signal
[pxx_message, f_message_psd] = pwelch(audio_resampled, [], [], [], fs, 'centered', 'power');
pxx_message_dB = 10 * log10(pxx_message + eps);  % Convert to dB scale

% Compute PSD for Modulated Signal
[pxx_modulated, f_modulated_psd] = pwelch(s_modulated, [], [], [], fs, 'centered', 'power');
pxx_modulated_dB = 10 * log10(pxx_modulated + eps);  % Convert to dB scale

% Compute PSD for Demodulated Signal
[pxx_demod, f_demod_psd] = pwelch(demod_filtered, [], [], [], fs, 'centered', 'power');
pxx_demod_dB = 10 * log10(pxx_demod + eps);  % Convert to dB scale

% Convert frequency axes to kHz for better readability
f_message_psd_kHz = f_message_psd / 1e3;
f_modulated_psd_kHz = f_modulated_psd / 1e3;
f_demod_psd_kHz = f_demod_psd / 1e3;

%% step 10: Plotting PSDs
figure;

% PSD of Original Message Signal
subplot(3, 1, 1); plot(f_message_psd_kHz, pxx_message_dB, 'b', 'LineWidth', 1.5); title('Power Spectral Density of Original Message Signal (dB)'); xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)'); grid on;
subplot(3, 1, 2); plot(f_modulated_psd_kHz, pxx_modulated_dB, 'r', 'LineWidth', 1.5); title('Power Spectral Density of Modulated Signal (dB)'); xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)'); grid on;

% PSD of Demodulated Signal
subplot(3, 1, 3); plot(f_demod_psd_kHz, pxx_demod_dB, 'm', 'LineWidth', 1.5); title('Power Spectral Density of Demodulated Signal (dB)'); xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)'); grid on;
disp('PSD plots for Original Message, Modulated Signal, and Demodulated Signal generated successfully.');

%% Step 11: Testing the Nulling Effect of Quadrature Carriers

% Generate the sine carrier (quadrature carrier) with the same frequency
sine_carrier = sin(2 * pi * fc * t(:));  
demod_signal_sine = s_modulated .* sine_carrier;  % Apply quadrature demodulation
% Design a proper low-pass filter if not already defined
% FIR filter design
% lp_filter = fir1(filter_order, filter_cutoff / (fs / 2));  % Normalize cutoff to Nyquist frequency
% Improved filtering with zero-phase filtfilt for better performance
% demod_filtered_sine = 2 * filtfilt(lp_filter, 1, demod_signal_sine);
% [b, a] = butter(2, filter_cutoff/(fs/2), 's'); % Analog filter design
% [num, den] = bilinear(b, a, fs); % Convert to digital equivalent
% demod_filtered_sine = 2 * filtfilt(num, den, demod_signal_sine); % Zero-phase filtering

[b, a] = butter(1, filter_cutoff/(fs/2), 'low');
demod_filtered_sine = 2 * filtfilt(b, a, demod_signal_sine);


% [sos, g] = butter(4, filter_cutoff/(fs/2), 'low', 's');
% demod_filtered_sine = 2 * filtfilt(sos, g, demod_signal_sine);

% Apply windowing to reduce spectral leakage
window_func = hamming(length(demod_filtered_sine));
windowed_signal = demod_filtered_sine .* window_func;

% Create time vector for plotting
t_demodulated = (0:length(demod_filtered_sine)-1) / fs;

%% step 12: Time-Frequency Analysis
N_fft = 2^nextpow2(length(demod_filtered_sine));
f_axis = (-fs/2:fs/N_fft:fs/2-fs/N_fft)/1e3;  
% Frequency Domain Analysis without dB conversion
psd_nulled = abs(fftshift(fft(windowed_signal, N_fft)))/N_fft;

figure; 
subplot(2,1,1); plot(t(1:min(1000,length(demod_filtered_sine))), demod_filtered_sine(1:min(1000,length(demod_filtered_sine))), 'r'); title('Time Domain - Demodulated Signal (Quadrature)'); xlabel('Time (s)'); ylabel('Amplitude'); grid on;
subplot(2,1,2); plot(f_axis, psd_nulled, 'r');  title('Frequency Spectrum - Demodulated Signal (Quadrature)'); xlabel('Frequency (kHz)'); ylabel('Magnitude'); grid on;

demod_resampled_sine = resample(demod_filtered_sine, fs_audio, fs);

[file, path] = uiputfile('*.wav', 'Save Nulled Audio As');
if file ~= 0
    % Save without normalization to preserve true nulling effect
    audiowrite(fullfile(path, file), demod_resampled_sine, fs_audio);
    disp(['Nulled audio saved to: ' fullfile(path, file)]);

    % Play the proper nulled audio (should be silent)
    sound(demod_resampled_sine, fs_audio);
else
    disp('Audio save cancelled.');
end
% Prompt user to choose where to save the filtered audio
[filename, pathname] = uiputfile('audio_no_drums.wav', 'Save Filtered Audio As');
if isequal(filename, 0) || isequal(pathname, 0)
    disp('User canceled the save operation.');
else
    filtered_audio_path = fullfile(pathname, filename);
    audiowrite(filtered_audio_path, audio_filtered, fs_original);
    disp(['Filtered Audio (Drums Minimized) saved to: ' filtered_audio_path]);
end
