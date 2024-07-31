folder_in = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/';
folder_out_raw_data = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITH_EXO/EUROBENCH/';
folder_out_preprocessed_data = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITH_EXO/PREPROCESSED/';
folder_out_figures = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITH_EXO/MULTIMEDIA/MATLAB/';
load(strcat(folder_in,'RESULTS_ExoStim.mat'));
sf = 2048;
muscles = {'time', 'ErSp_r', 'ReAb_r', 'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'ErSp_l', 'ReAb_l', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};

% fieldSuffix = 'Gait';
fieldSuffix = 'Exo';
% cond = "01";
cond = "02";

%% Raw signals
for subject = 1:5
    rawField = ['rawEMG', fieldSuffix];
    emg_raw = RESULTS(subject).(rawField);
    samples = 0:size(emg_raw,1)-1;
    time = samples*(1/sf);

    emg = [time', emg_raw]; 
    raw_emg = array2table(emg, 'VariableNames', muscles);

    eventsField = ['indexes', fieldSuffix];
    events = RESULTS(subject).(eventsField)/sf;

    for i=2:size(muscles,2)
        m = muscles{1,i};
        raw_signal = raw_emg.(m);

        figure;
        plot(raw_emg.time, raw_signal, 'r', 'LineWidth', 1.5);
        xlabel('Tiempo (s)');
        ylabel('Amplitud');
        title(strcat('Raw ',m));
        grid on;
        xline(events, 'g--', 'LineWidth', 2)
        legend('RAW EMG', 'Trial separator');
        savefig(gcf,strcat(folder_out_figures,'RAW/',m,'_',num2str(subject),'_raw.fig'))
        close all
    end

    writetable(raw_emg,strcat(folder_out_raw_data, 'subject_0', num2str(subject),'_cond_', cond, '_run_01_EMG.csv'))
end

%% Filtered signals
for subject = 1:5
    rawField = ['emgFiltered', fieldSuffix];
    emg_raw = RESULTS(subject).(rawField);
    samples = 0:size(emg_raw,1)-1;
    time = samples*(1/sf);

    emg = [time', emg_raw]; 
    raw_emg = array2table(emg, 'VariableNames', muscles);

    eventsField = ['indexes', fieldSuffix];
    events = RESULTS(subject).indexesGait/sf;

    cleaned_signal = table();
    cleaned_signal = addvars(cleaned_signal,raw_emg.time,'NewVariableNames','time');

    for i=2:size(muscles,2)
        m = muscles{1,i};
        raw_signal = raw_emg.(m);

        figure;
        plot(raw_emg.time, raw_signal, 'r', 'LineWidth', 1.5);
        xlabel('Tiempo (s)');
        ylabel('Amplitud');
        title(strcat('Filtered ',m));
        grid on;
        xline(events, 'g--', 'LineWidth', 2)
        legend('FILTERED EMG', 'Trial separator');
        savefig(gcf,strcat(folder_out_figures,'FILTERED/',m,'_',num2str(subject),'_filtered.fig'))
        close all

        % Apply median filter to the whole filtered signals
        filtered_emg = medfilt1(raw_signal, 6);

        threshold_factor = 8; % Factor de umbral
        threshold = threshold_factor * std(filtered_emg);

        %Detectar y eliminar los picos que superan el umbral
        cleaned_emg = filtered_emg;
        cleaned_emg(abs(filtered_emg) > threshold) = 0;

        %Plotear cleaned signals
        figure;
        plot(raw_emg.time, cleaned_emg, 'r', 'LineWidth', 1.5);
        xlabel('Tiempo (s)');
        ylabel('Amplitud');
        title(strcat('Cleaned ',m));
        grid on;
        xline(events, 'g--', 'LineWidth', 2)
        legend('CLEANED EMG', 'Trial separator');
        savefig(gcf,strcat(folder_out_figures,'CLEANED/',m,'_',num2str(subject),'_cleaned.fig'))
        close all

        cleaned_signal = addvars(cleaned_signal,cleaned_emg,'NewVariableNames',m);

        % Plotear filtered vs cleaned
        
        figure;
        plot(raw_emg.time, raw_signal, 'b', 'LineWidth', 1.5);

        hold on;
        plot(raw_emg.time, cleaned_emg, 'r', 'LineWidth', 1.5);
        xlabel('Tiempo (s)');
        ylabel('Amplitud');
        title(strcat('Filtered vs Cleaned ',m));
        grid on;
        xline(events, 'g--', 'LineWidth', 2)
        legend('FILTERED EMG','CLEANED EMG', 'Trial separator');
        savefig(gcf,strcat(folder_out_figures,'FILTERED vs CLEANED/',m,'_',num2str(subject),'.fig'))
        close all

    end

    writetable(raw_emg,strcat(folder_out_preprocessed_data, 'FILTERED/', 'subject_0', num2str(subject),'_cond_', cond, '_run_01_EMG.csv'))
    writetable(cleaned_signal,strcat(folder_out_preprocessed_data, 'CLEANED/', 'subject_0', num2str(subject),'_cond_', cond, '_run_01_EMG.csv'))
end


