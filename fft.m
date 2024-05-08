% Calculate FFT
name_file = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/DATA/gait_without_exo/PROCESSED/SUBJECT_0X/SESSION_01/IRREGULAR_TERRAIN/EMG/CLEANED (median filter)/subject_02_cond_01_run_01_EMG.csv'

data = readtable(name_file);

Raw_EMG = data.TiAn_l';
time = data.time;
[a, b] = size(Raw_EMG);

Fs = 2048; % Sampling frequency
T = 1/Fs; % Sample time
L = b; % Length of signal
t = (0:L-1)*T; % Time vector

for i=1:a
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    f = Fs/2*linspace(0,1,NFFT/2+1);
    FFT_Raw_EMG(i,:) = fft(Raw_EMG(i,:),NFFT)/L;
end

psd = 2*abs(FFT_Raw_EMG(1,1:length(f)));

plot(f,psd)
title('Single-Sided Amplitude Spectrum of Raw EMG')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
hold on

[~,peaks] = findpeaks(psd, 'MinPeakDistance',20000);
freq_peaks = f(peaks);

plot(f(peaks), psd(peaks), 'x', 'MarkerSize', 10, 'DisplayName', 'Picos');

filtered_emg = Raw_EMG;
for freq = freq_peaks
    fNorm = freq/(Fs/2); %'HP' is the cutoff frequency 
    bw = fNorm/2;
    if fNorm < 1
        [b,a] = iirnotch(fNorm,bw);
        filtered_emg = filter(b, a, filtered_emg);
    end
end


[a, b] = size(filtered_emg);

Fs = 2048; % Sampling frequency
T = 1/Fs; % Sample time
L = b; % Length of signal
t = (0:L-1)*T; % Time vector

for i=1:a
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    f = Fs/2*linspace(0,1,NFFT/2+1);
    FFT_filtered_EMG(i,:) = fft(filtered_emg(i,:),NFFT)/L;
end

psd = 2*abs(FFT_filtered_EMG(1,1:length(f)));

figure
plot(f,psd)
title('Single-Sided Amplitude Spectrum of Raw EMG')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
hold on

figure
plot(time, Raw_EMG)
hold on
plot(time, filtered_emg)
