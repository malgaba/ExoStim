folder_in = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITH_EXO/PROCESSED/ENVELOPE/NORMALIZED/';
folder_out = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITH_EXO/PROCESSED/ENVELOPE/NORMALIZED/';
files = dir(fullfile(folder_in,'*.mat'));
low = 1;
sf = 2048;


for i = 1:numel(files)
    total_data = [];
    if ~contains(files(i).name,'normalized')
        name_file = fullfile(folder_in, files(i).name);
        [~, nombre_archivo, ~] = fileparts(name_file)
        muscle =  nombre_archivo(length(nombre_archivo)-5:end);
        data = load(name_file);
        % if contains(muscle,'_l')
        %     data = data.segments_emg_l;
        % elseif contains(muscle,'_r')
        %      data = data.segments_emg_r;
        % end
        
        data = data.data;

        concatenation_points = [];

        for c = 1:size(data,2)
            total_data = [total_data; data(:,c)];
            concatenation_points = [concatenation_points; length(total_data)];
        end

        plot(total_data);
        title(nombre_archivo(length(nombre_archivo)-5:end));

        hold on;
        for cp = 1:length(concatenation_points)-1 % No incluimos el último punto porque es el final del total_data
            xline(concatenation_points(cp), '--r', 'LineWidth', 2);
        end
        hold off;

        % env = load(strcat(folder_in,nombre_archivo,'_normalized.mat'));
        % if contains(name_file,'_l')
        %     data = data.segments_emg_l;
        %     env = env.individual_cycles_l;
        % else
        %     data = data.segments_emg_r;
        %     env = env.individual_cycles_r;
        % end

        % Añadir la leyenda con el total de líneas
        legend(['Total líneas: ' num2str(length(concatenation_points)-1)], 'Location', 'best');


        saveas(gcf,strcat(folder_out,string(nombre_archivo), '.png'))
        savefig(gcf,strcat(folder_out,string(nombre_archivo), '.fig'))
        close all
        % save(strcat(folder_out, string(nombre_archivo)), 'data');
        % save(strcat(folder_out, string(nombre_archivo), '_normalized'),'env');
    end
end
%%
        
        
% 
%         signal = data - mean(data);
%         absSignal = abs(signal);
% 
%         % env = [];
%         band = (2/sf)*low;
%         [B, A] = butter(3, band, 'low');
%         % 
%         for c = 1:size(data,2)
%             cut = absSignal(:,c);
%             env = filtfilt(B, A, cut);
%             total_env = [total_env, env];
%         end 
% 
% 
%         mean_env = mean(env,2);
% 
%         figure;
%         x_scaled = linspace(0,100,length(mean_env));
%         plot(x_scaled, mean_env, 'b', 'LineWidth', 2, 'DisplayName','Mean envelope');
%         xlabel('% of gait cycle');
%         ylabel('Amplitud (V)');
% 
%         hold on;
%         for h = 1:size(data, 2)
%             plot(x_scaled, data(:, h), 'Color', [0.7 0.7 0.7]); % Ploteo en gris claro
%         end
% 
%         title(replace(nombre_archivo,"_","-"));
%         legend('Mean envelope')
%         hold off;
% 
%         % Save envelopes figures without normalize
%         saveas(gcf,strcat(folder_out,string(nombre_archivo), '_', m, '.png'))
% 
%         max_individual_cycles = [];
% 
%         % EMG data normalized by the average of its peaks from all cycles
%         for j = 1:size(data,2)
%             max_individual_cycles(j) = max(data(:,j));
%         end
% 
%         media = mean(max_individual_cycles);
%         individual_cycles = data/media;
% 
%         mean_env = mean(individual_cycles,2);
% 
%         figure;
%         x_scaled = linspace(0,100,length(mean_env));
%         plot(x_scaled, mean_env, 'b', 'LineWidth', 2, 'DisplayName','Mean envelope');
%         xlabel('% of gait cycle');
%         ylabel('Amplitud (V)');
% 
%         hold on;
%         for h = 1:size(individual_cycles, 2)
%             plot(x_scaled, individual_cycles(:, h), 'Color', [0.7 0.7 0.7]); % Ploteo en gris claro
%         end
% 
%         title(replace(nombre_archivo,"_","-"));
%         legend('Mean envelope')
%         hold off;
% 
%         % EMG sin normalizar
%         save(strcat(folder_out, string(nombre_archivo), '_', m),'data');
%         % ENVELOPES sin normalizar
%         % save(strcat(folder_out, string(nombre_archivo).replace('.csv',''), '_', m),'segments_env_r');
% 
%         % ENVELOPES normalizados
%         saveas(gcf,strcat(folder_out,string(nombre_archivo), '_', m, '_normalized', '.png'))
%         save(strcat(folder_out, string(nombre_archivo), '_', m, '_normalized'),'individual_cycles');
%         % writematrix(segments_r,strcat(folder_out, string(nombre_archivo).replace('.csv',''), '_', m, '.csv'))
%     end
% end
