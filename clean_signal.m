% Code to clean artifacts from the signal
% Aplicar un filtro mediano para suavizar la señal y preservar los picos
folder = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/DATA/gait_with_exo/PROCESSED/SUBJECT_0X/SESSION_01/IRREGULAR_TERRAIN/EMG/FILTERED';
files = dir(fullfile(folder,'*.csv'))
muscles = {'time', 'ErSp_r', 'ReAb_r', 'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'ErSp_l', 'ReAb_l', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};


for i = 1:numel(files)
    name_file = fullfile(folder, files(i).name);
    original_emg = readtable(name_file);
    [~, nombre_archivo, ~] = fileparts(name_file)

    cleaned_signal = table();
    cleaned_signal = addvars(cleaned_signal,original_emg.time,'NewVariableNames','time');
    
    % Graficar la señal original y la señal limpiada
    for c=2:size(muscles,2)
        m = muscles{1,c};
        
        filtered_emg = medfilt1(original_emg.(m), 11);
    
        threshold_factor = 3; % Factor de umbral
        threshold = threshold_factor * std(filtered_emg);
    
        %Detectar y eliminar los picos que superan el umbral
        cleaned_emg = filtered_emg;
        cleaned_emg(abs(filtered_emg) > threshold) = 0;
    
    
        figure;
        plot(original_emg.time, original_emg.(m), 'b', 'LineWidth', 1.5);
        hold on;
        plot(original_emg.time, cleaned_emg, 'r', 'LineWidth', 1.5);
        xlabel('Tiempo (s)');
        ylabel('Amplitud');
        title(strcat('Cleaned ',m));
        grid on;
        %xline(events, 'g--', 'LineWidth', 2)
        legend('Original EMG', 'Cleaned EMG');
        savefig(gcf,strcat(m,'_',nombre_archivo,'.fig'))
    
        cleaned_signal = addvars(cleaned_signal,cleaned_emg,'NewVariableNames',m);
    end
    writetable(cleaned_signal, strcat(nombre_archivo, '.csv'))
end
