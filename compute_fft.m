% Calculate FFT
in_path = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/CUT_SIGNALS';
out_path = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/CLEANED_PSD_SIGNALS/';
files = dir(fullfile(in_path,'*.csv'));
muscles = {'time', 'ErSp_r', 'ReAb_r', 'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'ErSp_l', 'ReAb_l', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};


for i = 1:numel(files)
    name_file = fullfile(in_path, files(i).name);
    %Leemos los archivos csv
    data = readtable(name_file);
    [~, nombre_archivo, ~] = fileparts(name_file)
    time = data.time;
    filtered_signal = table();
    filtered_signal = addvars(filtered_signal,time,'NewVariableNames','time');

    for c=2:size(muscles,2)
        m = muscles{1,c};
        % m = 'TiAn_l';
        emg = data.(m)';
        [a, b] = size(emg);
      

        Fs = 2048; % Sampling frequency
        T = 1/Fs; % Sample time
        L = b; % Length of signal
        t = (0:L-1)*T; % Time vector

        FFT_Raw_EMG = [];
        for s=1:a
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            f = Fs/2*linspace(0,1,NFFT/2+1);
            FFT_Raw_EMG(s,:) = fft(emg(s,:),NFFT)/L;
        end

        psd = 2*abs(FFT_Raw_EMG(1,1:length(f)));

        figure;
        plot(f,psd)
        title('Single-Sided Amplitude Spectrum of Filtered EMG')
        xlabel('Frequency (Hz)')
        ylabel('|Y(f)|')
        hold on

        [~,peaks] = findpeaks(psd, 'MinPeakDistance',5000, 'MinPeakHeight', max(psd) * 0.05);
        freq_peaks = f(peaks);
        peaks_2 = peaks(freq_peaks <= 500);
        freq_peaks_2 = f(peaks_2);
        
        plot(f(peaks_2), psd(peaks_2), 'x', 'MarkerSize', 10, 'DisplayName', 'Picos');
        
        filtered_emg = emg;
        for freq = freq_peaks_2
            fNorm = freq/(Fs/2);
            bw = fNorm/10;
            if fNorm < 1
                [b,a] = iirnotch(fNorm,bw);
                filtered_emg = filter(b, a, filtered_emg);
            end
        end
        
        [a, b] = size(filtered_emg);
        
        L = b; % Length of signal
        t = (0:L-1)*T; % Time vector
        
        FFT_filtered_EMG = [];
        for x=1:a
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            f = Fs/2*linspace(0,1,NFFT/2+1);
            FFT_filtered_EMG(x,:) = fft(filtered_emg(x,:),NFFT)/L;
        end
        
        psd = 2*abs(FFT_filtered_EMG(1,1:length(f)));
        

        plot(f,psd)
        title('Single-Sided Amplitude Spectrum of Raw EMG')
        xlabel('Frequency (Hz)')
        ylabel('|Y(f)|')
        legend ('Raw','peaks','filtered')
        title(m)

        figure
        plot(time, emg)
        hold on
        plot(time, filtered_emg)
        title(m)

        filtered_signal = addvars(filtered_signal, filtered_emg', 'NewVariableNames',m);
        
    end
    
    writetable(filtered_signal,strcat(out_path,'/',nombre_archivo,'.csv'))
    clearvars filtered_signal
    %input('Presiona Enter para continuar...');
    close all
end