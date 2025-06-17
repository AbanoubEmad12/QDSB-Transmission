%% Step 1: PARAMETERS 
fs_audio = 44100;        % sample rate
fs = 100e3;              % Modulation sample rate
fc = 20e3;               % Single carrier frequency for QDSB
fb = 20e3;               % Quadrature frequency
filter_order = 51;       % FIR filter order
filter_cutoff = 1e3;    % LPF cutoff frequency

%% Step 2: LOAD AUDIO FILES
[audio1, Fs1] = audioread('/Users/abanoubkamel/Downloads/Audio Files/Amr_Diab.mp3');
[audio2, Fs2] = audioread('/Users/abanoubkamel/Downloads/Audio Files/Angham.mp3');

%% Step 3: Convert stereo to mono
if size(audio1, 2) > 1
    audio1 = mean(audio1, 2); % Convert to mono by averaging channels
end
if size(audio2, 2) > 1
    audio2 = mean(audio2, 2); % Convert to mono by averaging channels
end

%% Step 4: Resample to match modulation rate
audio1 = resample(audio1, fs, Fs1);
audio2 = resample(audio2, fs, Fs2);

% Triming 
L = min(length(audio1), length(audio2));
audio1 = audio1(1:L);
audio2 = audio2(1:L);
t = (0:L-1)/fs;

%% Step 4: QDSB MODULATION 
carrier_I = cos(2*pi*fc*t);  
carrier_Q = sin(2*pi*fb*t);  
s_I = audio1' .* carrier_I;  
s_Q = audio2' .* carrier_Q;  
composite = s_I + s_Q;       

%% Step 5: ADD NOISE 
No = 10;
Noise = sqrt(No) * randn(1, length(composite)); 
received_QDSB = composite + Noise;  
SNR_dB = 10 * log10(var(composite) / No);
fprintf('Current SNR (dB): %.2f\n', SNR_dB);

%% Step 6: FIR LPF 
lpFilt = fir1(filter_order, 2*filter_cutoff/fs);  

%% Step 7: DEMODULATION (Noiseless) 
demod_I = composite .* (2 * cos(2*pi*fc*t));  % Demodulate I signal
demod_Q = composite .* (2 * sin(2*pi*fb*t));  % Demodulate Q signal
rec1 = filter(lpFilt, 1, demod_I);  % Filter I signal
rec2 = filter(lpFilt, 1, demod_Q);  % Filter Q signal

%% Step 8: DEMODULATION (Noisy) 
demod_I_noisy = received_QDSB .* (2 * cos(2*pi*fc*t));  % Demodulate I signal with noise
demod_Q_noisy = received_QDSB .* (2 * sin(2*pi*fb*t));  % Demodulate Q signal with noise
rec1_noisy = filter(lpFilt, 1, demod_I_noisy);  % Filter I signal with noise
rec2_noisy = filter(lpFilt, 1, demod_Q_noisy);  % Filter Q signal with noise

%% Step 9: TIME DOMAIN PLOTS 
figure;
subplot(4,2,1); plot(t, audio1); title('Original Message 1 (Time)'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(4,2,2); plot(t, audio2); title('Original Message 2 (Time)'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(4,2,3); plot(t, s_I); title('Modulated Signal 1 (I, Time)'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(4,2,4); plot(t, s_Q); title('Modulated Signal 2 (Q, Time)'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(4,2,5); plot(t, composite); title('Composite QDSB Signal (Time)'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(4,2,6); plot(t, rec1); title('Recovered Message 1 (Time)'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(4,2,7); plot(t, rec2); title('Recovered Message 2 (Time)'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(4,2,8); plot(t, rec1 + rec2); title('Sum of Recovered Messages (Time)'); xlabel('Time (s)'); ylabel('Amplitude');

%% Step 10: FREQUENCY DOMAIN 
n = length(t);  % Length of time vector
f = (-n/2:n/2-1)*(fs/n);  % Frequency vector
M1 = abs(fftshift(fft(audio1, n)));
M2 = abs(fftshift(fft(audio2, n)));
S1 = abs(fftshift(fft(s_I, n)));
S2 = abs(fftshift(fft(s_Q, n)));
C  = abs(fftshift(fft(composite, n)));
R1 = abs(fftshift(fft(rec1, n)));
R2 = abs(fftshift(fft(rec2, n)));

figure;
subplot(4,2,1); plot(f/1e3, M1); title('Original Message 1 (Freq)'); xlabel('Frequency (kHz)'); ylabel('|M1(f)|');
subplot(4,2,2); plot(f/1e3, M2); title('Original Message 2 (Freq)'); xlabel('Frequency (kHz)'); ylabel('|M2(f)|');
subplot(4,2,3); plot(f/1e3, S1); title('Modulated Signal 1 (I, Freq)'); xlabel('Frequency (kHz)'); ylabel('|S1(f)|');
subplot(4,2,4); plot(f/1e3, S2); title('Modulated Signal 2 (Q, Freq)'); xlabel('Frequency (kHz)'); ylabel('|S2(f)|');
subplot(4,2,5); plot(f/1e3, C); title('Composite QDSB Signal (Freq)'); xlabel('Frequency (kHz)'); ylabel('|C(f)|');
subplot(4,2,6); plot(f/1e3, R1); title('Recovered Message 1 (Freq)'); xlabel('Frequency (kHz)'); ylabel('|R1(f)|');
subplot(4,2,7); plot(f/1e3, R2); title('Recovered Message 2 (Freq)'); xlabel('Frequency (kHz)'); ylabel('|R2(f)|');
subplot(4,2,8); plot(f/1e3, R1 + R2); title('Sum of Recovered Messages (Freq)'); xlabel('Frequency (kHz)'); ylabel('|R1 + R2|');

%% Step 11: PSD (Power Spectral Density) 
[psd_audio1, f_psd] = pwelch(audio1, [], [], [], fs, 'centered', 'power');
[psd_audio2, f_audio2] = pwelch(audio2, [], [], [], fs, 'centered', 'power');
[psd_s_I, f_s_I] = pwelch(s_I, [], [], [], fs, 'centered', 'power');
[psd_s_Q, f_s_Q] = pwelch(s_Q, [], [], [], fs, 'centered', 'power');
[psd_composite, f_composite] = pwelch(composite, [], [], [], fs, 'centered', 'power');
[psd_rec1, f_rec1] = pwelch(rec1, [], [], [], fs, 'centered', 'power');
[psd_rec2, f_rec2] = pwelch(rec2, [], [], [], fs, 'centered', 'power');

% Convert to dB
psd_audio1_dB = 10 * log10(psd_audio1 + eps);
psd_audio2_dB = 10 * log10(psd_audio2 + eps);
psd_s_I_dB = 10 * log10(psd_s_I + eps);
psd_s_Q_dB = 10 * log10(psd_s_Q + eps);
psd_composite_dB = 10 * log10(psd_composite + eps);
psd_rec1_dB = 10 * log10(psd_rec1 + eps);
psd_rec2_dB = 10 * log10(psd_rec2 + eps);
psd_combined_dB = 10 * log10(psd_rec1 + psd_rec2 + eps);

figure;
subplot(4,2,1); plot(f_psd/1e3, psd_audio1_dB); title('Original Message 1 PSD'); xlabel('Freq (kHz)'); ylabel('dB/Hz'); grid on;
subplot(4,2,2); plot(f_audio2/1e3, psd_audio2_dB); title('Original Message 2 PSD'); xlabel('Freq (kHz)'); ylabel('dB/Hz'); grid on;
subplot(4,2,3); plot(f_s_I/1e3, psd_s_I_dB); title('Modulated Signal 1 PSD (I)'); xlabel('Freq (kHz)'); ylabel('dB/Hz'); grid on;
subplot(4,2,4); plot(f_s_Q/1e3, psd_s_Q_dB); title('Modulated Signal 2 PSD (Q)'); xlabel('Freq (kHz)'); ylabel('dB/Hz'); grid on;
subplot(4,2,5); plot(f_composite/1e3, psd_composite_dB); title('Composite QDSB Signal PSD'); xlabel('Freq (kHz)'); ylabel('dB/Hz'); grid on;
subplot(4,2,6); plot(f_rec1/1e3, psd_rec1_dB); title('Recovered Message 1 PSD'); xlabel('Freq (kHz)'); ylabel('dB/Hz'); grid on;
subplot(4,2,7); plot(f_rec2/1e3, psd_rec2_dB); title('Recovered Message 2 PSD'); xlabel('Freq (kHz)'); ylabel('dB/Hz'); grid on;
subplot(4,2,8); plot(f_rec2/1e3, psd_combined_dB); title('Combined Recovered PSD'); xlabel('Freq (kHz)'); ylabel('dB/Hz'); grid on;

[psd_received_QDSB, f_received] = pwelch(received_QDSB, [], [], [], fs, 'centered', 'power');
psd_received_QDSB_dB = 10 * log10(psd_received_QDSB + eps);

figure;
subplot(2,1,1); plot(f_composite/1e3, psd_composite_dB); title('Transmitted QDSB Signal PSD'); xlabel('Freq (kHz)'); ylabel('dB/Hz'); grid on;
subplot(2,1,2); plot(f_received/1e3, psd_received_QDSB_dB); title('Received QDSB Signal PSD (with Noise)'); xlabel('Freq (kHz)'); ylabel('dB/Hz'); grid on;

%% Step 12: SAVE AUDIO FILES 
% Create output directory if it doesn't exist
output_dir = 'output_audio';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

audio1_original = resample(audio1, fs_audio, fs);
audio2_original = resample(audio2, fs_audio, fs);
% Recovered audio (without noise)
rec1_audio = resample(rec1', fs_audio, fs);
rec2_audio = resample(rec2', fs_audio, fs);
% Recovered audio (with noise)
rec1_noisy_audio = resample(rec1_noisy', fs_audio, fs);
rec2_noisy_audio = resample(rec2_noisy', fs_audio, fs);
% Normalize audio to prevent clipping
normalize = @(x) x / max(abs(x));
audio1_original = normalize(audio1_original);
audio2_original = normalize(audio2_original);
rec1_audio = normalize(rec1_audio);
rec2_audio = normalize(rec2_audio);
rec1_noisy_audio = normalize(rec1_noisy_audio);
rec2_noisy_audio = normalize(rec2_noisy_audio);

% % Ask the user to choose the output directory
output_dir = uigetdir('', 'Select Folder to Save Audio Files');

% Check if the user pressed Cancel (i.e., output_dir is 0)
if output_dir == 0
    disp('No directory selected. Exiting...');
    return;
end

% Save original audio files
audiowrite(fullfile(output_dir, 'original_message1.wav'), audio1_original, fs_audio);
audiowrite(fullfile(output_dir, 'original_message2.wav'), audio2_original, fs_audio);

% Save recovered audio files (without noise)
audiowrite(fullfile(output_dir, 'recovered_message1_clean.wav'), rec1_audio, fs_audio);
audiowrite(fullfile(output_dir, 'recovered_message2_clean.wav'), rec2_audio, fs_audio);

% Save recovered audio files (with noise)
audiowrite(fullfile(output_dir, 'recovered_message1_noisy.wav'), rec1_noisy_audio, fs_audio);
audiowrite(fullfile(output_dir, 'recovered_message2_noisy.wav'), rec2_noisy_audio, fs_audio);

fprintf('All audio files have been saved to the "%s" directory.\n', output_dir);


%% --- PHASE ERROR ANALYSIS ---
% More efficient implementation using a loop
phase_errors = [0, 10, 45, 90];
% for i = 1:length(phase_errors)
%     phi_rad = phase_errors(i) * (pi/180);
% 
%     % Demodulate with phase error
%     demod_I_phase = composite .* (2 * cos(2*pi*fc*t + phi_rad));
%     demod_Q_phase = composite .* (2 * sin(2*pi*fb*t + phi_rad));
% 
%     % Filter
%     rec1_phase = filter(lpFilt, 1, demod_I_phase);
%     rec2_phase = filter(lpFilt, 1, demod_Q_phase);
% 
%     % Create a new figure for each phase error
%     figure;
%     subplot(4,2,1); plot(t, audio1); title('Original Message 1 (Time)'); xlabel('Time (s)'); ylabel('Amplitude');
%     subplot(4,2,2); plot(t, audio2); title('Original Message 2 (Time)'); xlabel('Time (s)'); ylabel('Amplitude');
%     subplot(4,2,3); plot(t, s_I); title('Modulated Signal 1 (I, Time)'); xlabel('Time (s)'); ylabel('Amplitude');
%     subplot(4,2,4); plot(t, s_Q); title('Modulated Signal 2 (Q, Time)'); xlabel('Time (s)'); ylabel('Amplitude');
%     subplot(4,2,5); plot(t, composite); title('Composite QDSB Signal (Time)'); xlabel('Time (s)'); ylabel('Amplitude');
%     subplot(4,2,6); plot(t, rec1_phase); title(['Recovered Message 1 (φ = ', num2str(phase_errors(i)), '°)']); xlabel('Time (s)'); ylabel('Amplitude');
%     subplot(4,2,7); plot(t, rec2_phase); title(['Recovered Message 2 (φ = ', num2str(phase_errors(i)), '°)']); xlabel('Time (s)'); ylabel('Amplitude');
%     subplot(4,2,8); plot(t, rec1_phase + rec2_phase); title(['Sum of Recovered Messages (φ = ', num2str(phase_errors(i)), '°)']); xlabel('Time (s)'); ylabel('Amplitude');
% end

% More efficient implementation using a loop for frequency domain analysis


for i = 1:length(phase_errors)
    phi_rad = phase_errors(i) * (pi/180);

    % Demodulate with phase error
    demod_I_phase = composite .* (2 * cos(2*pi*fc*t + phi_rad));
    demod_Q_phase = composite .* (2 * sin(2*pi*fb*t + phi_rad));

    % Filter
    rec1_phase = filter(lpFilt, 1, demod_I_phase);
    rec2_phase = filter(lpFilt, 1, demod_Q_phase);

    % Compute frequency domain representations
    n = length(t);
    f = (-n/2:n/2-1)*(fs/n);  % Frequency vector

    M1 = abs(fftshift(fft(audio1, n)));
    M2 = abs(fftshift(fft(audio2, n)));
    S1 = abs(fftshift(fft(s_I, n)));
    S2 = abs(fftshift(fft(s_Q, n)));
    C  = abs(fftshift(fft(composite, n)));
    R1_phase = abs(fftshift(fft(rec1_phase, n)));
    R2_phase = abs(fftshift(fft(rec2_phase, n)));
    R_sum = abs(fftshift(fft(rec1_phase + rec2_phase, n)));

    % Create a new figure for each phase error
    figure;

    % Apply logarithmic scale to better visualize differences in the frequency domain
subplot(4,2,1); semilogy(f/1e3, M1); title('Original Message 1 (Freq)'); xlabel('Frequency (kHz)'); ylabel('|M1(f)|'); grid on;
subplot(4,2,2); semilogy(f/1e3, M2); title('Original Message 2 (Freq)'); xlabel('Frequency (kHz)'); ylabel('|M2(f)|'); grid on;
subplot(4,2,3); semilogy(f/1e3, S1); title('Modulated Signal 1 (I, Freq)'); xlabel('Frequency (kHz)'); ylabel('|S1(f)|'); grid on;
subplot(4,2,4); semilogy(f/1e3, S2); title('Modulated Signal 2 (Q, Freq)'); xlabel('Frequency (kHz)'); ylabel('|S2(f)|'); grid on;
subplot(4,2,5); semilogy(f/1e3, C); title('Composite QDSB Signal (Freq)'); xlabel('Frequency (kHz)'); ylabel('|C(f)|'); grid on;
subplot(4,2,6); semilogy(f/1e3, R1_phase); title(['Recovered Message 1 (φ = ', num2str(phase_errors(i)), '°)']); xlabel('Frequency (kHz)'); ylabel('|R1(f)|'); grid on;
subplot(4,2,7);semilogy(f/1e3, R2_phase);title(['Recovered Message 2 (φ = ', num2str(phase_errors(i)), '°)']);xlabel('Frequency (kHz)');ylabel('|R2(f)|');grid on;
subplot(4,2,8); semilogy(f/1e3, R_sum);title(['Sum of Recovered Messages (φ = ', num2str(phase_errors(i)), '°)']);xlabel('Frequency (kHz)');ylabel('|R1+R2(f)|');grid on;
end

