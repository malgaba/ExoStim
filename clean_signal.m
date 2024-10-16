% Code to clean artifacts from the signal
% Aplicar un filtro mediano para suavizar la señal y preservar los picos
folder = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/CLEANED_PSD_SIGNALS/';
folder_out = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/MEDIAN_FILTERED/';
files = dir(fullfile(folder,'*.csv'));
muscles = {'time', 'ErSp_r', 'ReAb_r', 'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'ErSp_l', 'ReAb_l', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};
b_plot = true;

for i = 1:numel(files)
    if ~contains(files(i).name,'gaitEvents')
        name_file = fullfile(folder, files(i).name);
        original_emg = readtable(name_file);
        [~, nombre_archivo, ~] = fileparts(name_file)
    
        cleaned_signal = table();
        cleaned_signal = addvars(cleaned_signal,original_emg.time,'NewVariableNames','time');
        
        % Graficar la señal original y la señal limpiada
        for c=2:size(muscles,2)
            m = muscles{1,c};
            data = original_emg.(m);
            % MEDIAN FILTER
            filtered_emg = medfilt1(original_emg.(m),6);

            threshold_factor = 8; % Factor de umbral
            threshold = threshold_factor * std(filtered_emg);

            %Detectar y eliminar los picos que superan el umbral
            cleaned_emg = filtered_emg;
            cleaned_emg(abs(filtered_emg) > threshold) = 0;
            
            % ICA
            % Aplicar ICA a la matriz de señales
            % [ica_results, A, W] = runica(original_emg.(m), 'lastEig', 3);
            % 
            % % Seleccionar la señal que contiene el artefacto
            % artifact = ica_results(1, :);
            % 
            % % Remover el artefacto de la señal EMG
            % cleaned_emg = original_emg.(m) - artifact;
    
            %UMBRALIZATION
            % umbral = std(original_emg.(m)) * 3;  % Multiplicado por un factor para ajustar el umbral
            % umbral = 0.0005;
            % cleaned_emg = original_emg.(m);  % Copia de seguridad de los datos originales
            % cleaned_emg(abs(original_emg.(m)) > umbral) = 0;  
    
            % % Calcular la diferencia de amplitud entre muestras consecutivas
            % senal_filtrada =  sgolayfilt(original_emg.(m), 3, 21);
            % dif_amplitud = diff(senal_filtrada);
            % 
            % % Definir un umbral basado en esta diferencia (por ejemplo, un múltiplo de la desviación estándar de la diferencia)
            % factor_umbral = 2;
            % umbral = std(dif_amplitud) * factor_umbral;
            % 
            % % Aplicar el umbral para identificar y eliminar el ruido
            % cleaned_emg = original_emg.(m); 
            % for i = 1:length(dif_amplitud)
            %     if abs(dif_amplitud(i)) > umbral
            %         cleaned_emg(i+1) = cleaned_emg(i);  % Reemplazar el siguiente valor con el valor actual si la diferencia supera el umbral
            %     end
            % end
    
            if b_plot
                figure;
                % signal = original_emg.(m)(1:length(original_emg.(m))-1,1);
                time = original_emg.time;
                plot(time, original_emg.(m), 'b');
                hold on;
                % % plot(original_emg.time, dif_amplitud, 'r', 'LineWidth', 1.5);
                plot(time, cleaned_emg, 'r', 'LineWidth', 1.5);
                xlabel('Tiempo (s)');
                ylabel('Amplitud');
                title(strcat('Filtered ',m));
                grid on;
                % %xline(events, 'g--', 'LineWidth', 2)
                legend('Filtered EMG', 'Cleaned EMG');
            end
            savefig(gcf,strcat(folder_out,nombre_archivo,'_',m,'.fig'))
           
    
            % figure;
            % plot(time, cleaned_emg);
            % title('Selecciona manualmente la plantilla');
            % xlabel('Tiempo');
            % ylabel('Amplitud');
            % 
            % fprintf('Haz clic en la gráfica para seleccionar los puntos de la plantilla.\n');
            % fprintf('Presiona Enter cuando hayas terminado.\n');
            % [x, ~] = ginput;
            % 
            % % Convertir las coordenadas x a índices de muestras
            % indices_muestras = round(x);
            % 
            % % Extraer la plantilla de la señal EMG
            % plantilla = cleaned_emg(indices_muestras(1):indices_muestras(end));
            % 
            % % Mostrar la plantilla seleccionada
            % figure;
            % plot(plantilla);
            % title('Plantilla seleccionada');
            % xlabel('Muestras');
            % ylabel('Amplitud');
        
            cleaned_signal = addvars(cleaned_signal,cleaned_emg,'NewVariableNames',m);
        end
    % disp('Pulse Enter para continuar...');
    % pause_input = input('','s');
    close all
    writetable(cleaned_signal, strcat(folder_out, nombre_archivo, '.csv'))
    end
end
