load('RESULTS_ExoStim.mat');
sf = 2048;
muscles = {'time', 'ErSp_r', 'ReAb_r', 'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'ErSp_l', 'ReAb_l', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};


for subject = 1:5
    emg_gait = RESULTS(subject).emgFilteredExo;
    samples = 0:size(emg_gait,1)-1;
    time = samples*(1/sf);
    
    emg = [time', emg_gait]; 

    raw_emg = array2table(emg, 'VariableNames', muscles);

    cleaned_signal = table();
    cleaned_signal = addvars(cleaned_signal,raw_emg.time,'NewVariableNames','time');
    

    %writetable(raw_emg, strcat('subject_0',num2str(subject),'_cond_02_run_01.csv'))

    events = RESULTS(subject).indexesExo/sf;

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
        % plot(raw_emg.time, original_emg, 'b', 'LineWidth', 1.5);
        % hold on;
        plot(raw_emg.time, cleaned_emg, 'r', 'LineWidth', 1.5);
        xlabel('Tiempo (s)');
        ylabel('Amplitud');
        title(strcat('Cleaned ',m));
        grid on;
        xline(events, 'g--', 'LineWidth', 2)
        legend('Original EMG', 'Cleaned EMG');
        savefig(gcf,strcat(m,'_',num2str(subject),'_cleaned.fig'))

        cleaned_signal = addvars(cleaned_signal,cleaned_emg,'NewVariableNames',m);
    end
    % writetable(cleaned_signal, strcat('subject_0',num2str(subject),'_cond_02_run_01_cleaned.csv'))
end


