load('RESULTS_ExoStim.mat');
sf = 2048;
muscles = {'time', 'ErSp_r', 'ReAb_r', 'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'ErSp_l', 'ReAb_l', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};


for subject = 1
    emg_gait = RESULTS(subject).emgFilteredExo;
    samples = 0:size(emg_gait,1)-1;
    time = samples*(1/sf);
    
    emg = [time', emg_gait]; 

    raw_emg = array2table(emg, 'VariableNames', muscles);

    cleaned_signal = table();
    cleaned_signal = addvars(cleaned_signal,raw_emg.time,'NewVariableNames','time');
    

    %writetable(raw_emg, strcat('subject_0',num2str(subject),'_cond_02_run_01.csv'))

    events = RESULTS(subject).indexesExo/sf;

    % for i=2:size(muscles,2)
    %     m = muscles{1,i};
    %     figure;
    %     plot(raw_emg.time05, raw_emg.(m));
    %     title(m)
    %     %disp(strcat(m,'_',num2str(subject)))
    %     savefig(gcf,strcat(m,'_',num2str(subject),'.fig'))
    % end

    % Aplicar un filtro de media móvil para suavizar la señal
    % window_size = 0.1; % Tamaño de la ventana móvil en segundos
    % window_length = round(window_size * sf);
    % filtered_emg = movmean(raw_emg.ErSp_r, window_length); % Aplicar el filtro de media móvil
    % 
    % % Graficar la señal original y la señal filtrada
    % figure;
    % plot(raw_emg.time, raw_emg.ErSp_r, 'b', 'LineWidth', 1.5);
    % hold on;
    % plot(raw_emg.time, filtered_emg, 'r', 'LineWidth', 1.5);
    % xlabel('Tiempo (s)');
    % ylabel('Amplitud');
    % title('Eliminación de artefactos de una señal EMG');
    % legend('Señal EMG original con artefactos', 'Señal EMG filtrada');
    % grid on;

    % % Aplicar un filtro mediano para eliminar los artefactos bruscos
    % filtered_emg = medfilt1(raw_emg.ErSp_r, 11); % Aplicar filtro mediano con ventana de 11 muestras
    % 
    % % Graficar la señal original y la señal filtrada
    % figure;
    % plot(raw_emg.time, raw_emg.ErSp_r, 'b', 'LineWidth', 1.5);
    % hold on;
    % plot(raw_emg.time, filtered_emg, 'r', 'LineWidth', 1.5);
    % xlabel('Tiempo (s)');
    % ylabel('Amplitud');
    % title('Eliminación de artefactos bruscos de una señal EMG');
    % legend('Señal EMG original con artefactos bruscos', 'Señal EMG filtrada');
    % grid on;

        % Aplicar un filtro mediano para suavizar la señal y preservar los picos
    original_emg = raw_emg.GlMe_r;
    filtered_emg = medfilt1(original_emg, 11); % Aplicar filtro mediano con ventana de 11 muestras
    
    % Calcular el umbral adaptativo como un múltiplo de la desviación estándar
    threshold_factor = 3; % Factor de umbral
    threshold = threshold_factor * std(filtered_emg);
    
    % Detectar y eliminar los picos que superan el umbral
    cleaned_emg = filtered_emg;
    cleaned_emg(abs(filtered_emg) > threshold) = 0;
    
    % Graficar la señal original y la señal limpiada
    for i=2:size(muscles,2)
        m = muscles{1,i};
        original_emg = raw_emg.(m);
        filtered_emg = medfilt1(original_emg, 11);

        threshold_factor = 3; % Factor de umbral
        threshold = threshold_factor * std(filtered_emg);

        %Detectar y eliminar los picos que superan el umbral
        cleaned_emg = filtered_emg;
        cleaned_emg(abs(filtered_emg) > threshold) = 0;

        figure;
        plot(raw_emg.time, original_emg, 'b', 'LineWidth', 1.5);
        hold on;
        plot(raw_emg.time, cleaned_emg, 'r', 'LineWidth', 1.5);
        xlabel('Tiempo (s)');
        ylabel('Amplitud');
        title(strcat('Cleaned ',m));
        grid on;
        xline(events, 'g--', 'LineWidth', 2)
        % legend('Original EMG', 'Cleaned EMG');
        % savefig(gcf,strcat(m,'_',num2str(subject),'_cleaned.fig'))

        cleaned_signal = addvars(cleaned_signal,cleaned_emg,'NewVariableNames',m);
    end
    %writetable(raw_emg, strcat('subject_0',num2str(subject),'_cond_02_run_01_cleaned.csv'))
end


