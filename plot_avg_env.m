%% Plot average envelopes
folder_in = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/CONCAT_SIGNALS/';
folder_out = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/ENVELOPE/';
files = dir(fullfile(folder_in,'*.csv'));
muscles = {'time', 'GaLa_l', 'TiAn_l'};
low = 1;
sf = 2048;

% retired_cycles = table();

for i = 1:numel(files)
    if ~contains(files(i).name,'estimatedGaitEvents')
        name_file = fullfile(folder_in, files(i).name);
        data = readtable(name_file);
        [~, nombre_archivo, ~] = fileparts(name_file)
    
        events = readtable(strcat(folder_in,nombre_archivo,'_estimatedGaitEvents.csv'));
        
        % temp_table = table();
        % temp_table = addvars(temp_table, string(nombre_archivo), 'NewVariableNames','Retired cycles');
    
        time = data.time;
    
        % t = table();
        % t = addvars(t,time,'NewVariableNames','time');

        t_m = table();

        segments_emg_l = [];
        segments_env_l = [];
        segments_emg_r = [];
        segments_env_r = [];

        for s=2:size(muscles,2)
            m = muscles{1,s};
            signal = data.(m);
    
            signal = signal - mean(signal);
            absSignal = abs(signal);
    
            band = (2/sf)*low;
            [B, A] = butter(3, band, 'low');
            env = filtfilt(B, A, absSignal);
            % env = absSignal;
    
            % t = addvars(t,env,'NewVariableNames',m);
            
    
            if contains(m, '_r')
                % figure;
                % plot(time, env, 'LineWidth', 2);
                % for c = 1:length(events.heel_strike_r)
                %     line([events.heel_strike_r(c), events.heel_strike_r(c)],ylim,'Color', 'r');
                % end
                % title(strcat(m, '(', string(nombre_archivo).replace('_','-'), ')'))
                % closest_curve = [];
                % closest_curve = input('Introduzca el número de segmento que quiere eliminar (separados por comas o espacios)...', 's');
                % 
                % closest_curve = str2num(closest_curve);
                % 
                % temp_table = addvars(temp_table, {closest_curve}, 'NewVariableNames',m);

                % close all

                for r=1:length(events.heel_strike_r)-1
                    mean_cycle = mean(events.heel_strike_r(r+1) - events.heel_strike_r(r));
                    std_cycle = std(events.heel_strike_r(r+1) - events.heel_strike_r(r));
                   if events.heel_strike_r(r+1) - events.heel_strike_r(r) <= (mean_cycle-2*std_cycle)
                        st_wind = find(time >= events.heel_strike_r(r), 1);
                        end_wind = find(time >= events.heel_strike_r(r+1), 1);           
                        new_env = env(st_wind:end_wind);
                        new_emg = signal(st_wind:end_wind);
                        % new_time = linspace(0,length(new_emg));

                        count = 1;
                        segments = [];

                        new_yi_length = length(new_emg);
                        xi = 1:new_yi_length;
                        ai = linspace(1,new_yi_length,200);
                        for ji = 1:size(ai,2)
                            t1i = ai(ji);
                            segments_env(ji,count) = interp1(xi',new_env(1:new_yi_length),t1i);
                            segments_emg(ji,count) = interp1(xi',new_emg(1:new_yi_length),t1i);
                        end
                        count = count+1;

                        % if length(new_emg)<size(segments_r,1)
                        %     new_emg = [new_emg; NaN(size(segments_r,1) - length(new_emg), 1)];
                        % elseif length(new_emg)>size(segments_r,1)
                        %     new_table = [];
                        %     for column = 1:size(segments_r,2)
                        %         segments_r_padded = [segments_r(:,column); NaN(length(new_emg) - size(segments_r,1), 1)];
                        %         new_table(:,column) = segments_r_padded;
                        %     end
                        %     segments_r = [];
                        %     segments_r = new_table;
                        % end
                        % 
                        segments_env_r(:,r) = segments_env;
                        segments_emg_r(:,r) = segments_emg;
                        % 
                        % plot(new_emg)
                        % title(strcat(m, '(', string(nombre_archivo).replace('_','-'), ')'))
                        % hold on;
                   end
                end
                % disp('Seleccione puntos en las curvas que desea eliminar. Pulse Enter cuando haya terminado.');
                % [x, y] = ginput(1); % Seleccionar puntos en la gráfica
                % 
                % columns_to_delete = [];
                % for p = 1:length(x)
                %     % Encontrar la curva más cercana al punto seleccionado
                %     min_dist = inf;
                %     closest_curve = -1;
                %     for j = 1:size(segments_r, 2)
                %         %dist = sum(abs(segments_r(:,j) - y(p)));
                %         dist = min(abs((1:size(segments_r, 1))' - x(p)) + abs(segments_r(:, j) - y(p)));
                %         if dist < min_dist
                %             min_dist = dist;
                %             closest_curve = j;
                %         end
                %     end
                % if closest_curve ~= -1
                %     segments_r(:, closest_curve) = [];
                % end

                
                
                % % close all
                mean_env = mean(segments_env_r,2);
                t_m = addvars(t_m, mean_env,'NewVariableNames',m);

                % Ploteo de la nueva media con una línea más gruesa
                figure;
                x_scaled = linspace(0,100,length(mean_env));
                plot(x_scaled, mean_env, 'b', 'LineWidth', 2, 'DisplayName','Mean envelope');
                xlabel('% of gait cycle');
                ylabel('Amplitud (V)');


                hold on;
                for h = 1:size(segments_env_r, 2)
                    plot(x_scaled, segments_env_r(:, h), 'Color', [0.7 0.7 0.7]); % Ploteo en gris claro
                end

                title(m);
                legend('Mean envelope')
                hold off;

                saveas(gcf,strcat(folder_out,string(nombre_archivo).replace('.csv',''), '_', m, '.png'))
                close all;

                % EMG data normalized by the average of its peaks from all cycles
                max_individual_cycles = [];
                for j = 1:size(segments_env_r,2)
                    max_individual_cycles(j) = max(segments_env_r(:,j));
                end

                % hold on; plot(mean(data,2), 'b', 'LineWidth', 2, 'DisplayName','Mean envelope');

                media = mean(max_individual_cycles);
                individual_cycles_r = segments_env_r/media;

                mean_env = mean(individual_cycles_r,2);

                figure;
                x_scaled = linspace(0,100,length(mean_env));
                plot(x_scaled, mean_env, 'b', 'LineWidth', 2, 'DisplayName','Mean envelope');
                xlabel('% of gait cycle');
                ylabel('Amplitud (V)');
            
                hold on;
                for h = 1:size(individual_cycles_r, 2)
                    plot(x_scaled, individual_cycles_r(:, h), 'Color', [0.7 0.7 0.7]); % Ploteo en gris claro
                end
                               
                title(strcat(replace(nombre_archivo,"_","-"),'-',m));
                legend('Mean envelope normalized')
                hold off;

                % savefig(gcf,strcat(folder_out,string(nombre_archivo).replace('.csv',''), '_', m, '.fig'))
                
                % EMG sin normalizar
                save(strcat(folder_out, string(nombre_archivo), '_', m),'segments_emg_r');
                % ENVELOPES sin normalizar
                % save(strcat(folder_out, string(nombre_archivo).replace('.csv',''), '_', m),'segments_env_r');

                % ENVELOPES normalizados
                saveas(gcf,strcat(folder_out,string(nombre_archivo), '_', m, '_normalized', '.png'))
                save(strcat(folder_out, string(nombre_archivo), '_', m, '_normalized'),'individual_cycles_r');
                % writematrix(segments_r,strcat(folder_out, string(nombre_archivo).replace('.csv',''), '_', m, '.csv'))

            end

            if contains(m, '_l')
                % figure;
                % plot(time, env, 'LineWidth', 2);
                % for c = 1:length(events.heel_strike_l)
                %     line([events.heel_strike_l(c), events.heel_strike_l(c)],ylim,'Color', 'r');
                % end
                % title(strcat(m, '(', string(nombre_archivo).replace('_','-'), ')'))
                % closest_curve = [];
                % closest_curve = input('Introduzca el número de segmento que quiere eliminar (separados por comas o espacios)...', 's');
                % 
                % closest_curve = str2num(closest_curve);
                % 
                % temp_table = addvars(temp_table, {closest_curve}, 'NewVariableNames',m);

                % close all

                for r=1:length(events.heel_strike_l)-1
                    mean_cycle = mean(events.heel_strike_l(r+1) - events.heel_strike_l(r));
                    std_cycle = std(events.heel_strike_l(r+1) - events.heel_strike_l(r));
                   if events.heel_strike_l(r+1) - events.heel_strike_l(r) <= (mean_cycle-2*std_cycle)
                        st_wind = find(time >= events.heel_strike_l(r), 1);
                        end_wind = find(time >= events.heel_strike_l(r+1), 1);           
                        new_env = env(st_wind:end_wind);
                        new_emg = signal(st_wind:end_wind);
                        % new_time = linspace(0,length(new_emg));

                        count = 1;
                        segments = [];

                        new_yi_length = length(new_emg);
                        xi = 1:new_yi_length;
                        ai = linspace(1,new_yi_length,200);
                        for ji = 1:size(ai,2)
                            t1i = ai(ji);
                            segments_env(ji,count) = interp1(xi',new_env(1:new_yi_length),t1i);
                            segments_emg(ji,count) = interp1(xi',new_emg(1:new_yi_length),t1i);
                        end
                        count = count+1;

                        % if length(new_emg)<size(segments_l,1)
                        %     new_emg = [new_emg; NaN(size(segments_l,1) - length(new_emg), 1)];
                        % elseif length(new_emg)>size(segments_l,1)
                        %     new_table = [];
                        %     for column = 1:size(segments_l,2)
                        %         segments_l_padded = [segments_l(:,column); NaN(length(new_emg) - size(segments_l,1), 1)];
                        %         new_table(:,column) = segments_l_padded;
                        %     end
                        %     segments_l = [];
                        %     segments_l = new_table;
                        % end
                        % 
                        segments_env_l(:,r) = segments_env;
                        segments_emg_l(:,r) = segments_emg;
                        % 
                        % plot(new_emg)
                        % title(strcat(m, '(', string(nombre_archivo).replace('_','-'), ')'))
                        % hold on;
                   end
                end
                % disp('Seleccione puntos en las curvas que desea eliminar. Pulse Enter cuando haya terminado.');
                % [x, y] = ginput(1); % Seleccionar puntos en la gráfica
                % 
                % columns_to_delete = [];
                % for p = 1:length(x)
                %     % Encontrar la curva más cercana al punto seleccionado
                %     min_dist = inf;
                %     closest_curve = -1;
                %     for j = 1:size(segments_l, 2)
                %         %dist = sum(abs(segments_l(:,j) - y(p)));
                %         dist = min(abs((1:size(segments_l, 1))' - x(p)) + abs(segments_l(:, j) - y(p)));
                %         if dist < min_dist
                %             min_dist = dist;
                %             closest_curve = j;
                %         end
                %     end
                % if closest_curve ~= -1
                %     segments_l(:, closest_curve) = [];
                % end

                
                
                % close all
                mean_env = mean(segments_env_l,2);
                t_m = addvars(t_m, mean_env,'NewVariableNames',m);

                % Ploteo de la nueva media con una línea más gruesa
                figure;
                x_scaled = linspace(0,100,length(mean_env));
                plot(x_scaled, mean_env, 'b', 'LineWidth', 2, 'DisplayName','Mean envelope');
                xlabel('% of gait cycle');
                ylabel('Amplitud (V)');


                hold on;
                for h = 1:size(segments_env_l, 2)
                    plot(x_scaled, segments_env_l(:, h), 'Color', [0.7 0.7 0.7]); % Ploteo en gris claro
                end

                title(m);
                legend('Mean envelope')
                hold off;

                saveas(gcf,strcat(folder_out,string(nombre_archivo), '_', m, '.png'))
                close all;

                % EMG data normalized by the average of its peaks from all cycles
                max_individual_cycles = [];
                for j = 1:size(segments_env_l,2)
                    max_individual_cycles(j) = max(segments_env_l(:,j));
                end

                % hold on; plot(mean(data,2), 'b', 'LineWidth', 2, 'DisplayName','Mean envelope');

                media = mean(max_individual_cycles);
                individual_cycles_l = segments_env_l/media;

                mean_env = mean(individual_cycles_l,2);

                figure;
                x_scaled = linspace(0,100,length(mean_env));
                plot(x_scaled, mean_env, 'b', 'LineWidth', 2, 'DisplayName','Mean envelope');
                xlabel('% of gait cycle');
                ylabel('Amplitud (V)');
            
                hold on;
                for h = 1:size(individual_cycles_l, 2)
                    plot(x_scaled, individual_cycles_l(:, h), 'Color', [0.7 0.7 0.7]); % Ploteo en gris claro
                end
                               
                title(strcat(replace(nombre_archivo,"_","-"),'-',m));
                legend('Mean envelope normalized')
                hold off;

                % savefig(gcf,strcat(folder_out,string(nombre_archivo).replace('.csv',''), '_', m, '.fig'))
                
                % EMG sin normalizar
                save(strcat(folder_out, string(nombre_archivo), '_', m),'segments_emg_l');
                % ENVELOPES sin normalizar
                % save(strcat(folder_out, string(nombre_archivo).replace('.csv',''), '_', m),'segments_env_r');

                % ENVELOPES normalizados
                saveas(gcf,strcat(folder_out,string(nombre_archivo), '_', m, '_normalized', '.png'))
                save(strcat(folder_out, string(nombre_archivo), '_', m, '_normalized'),'individual_cycles_l');
                % writematrix(segments_r,strcat(folder_out, string(nombre_archivo).replace('.csv',''), '_', m, '.csv'))
            end

        end
        % retired_cycles = vertcat(retired_cycles,temp_table);
        disp('Pulse Enter para continuar...');
        pause_input = input('','s');
        close all
    
    end

    % writetable(t,strcat(folder_out, string(nombre_archivo).replace('.csv',''), '_env.csv'))
    % writetable(t_m,strcat(folder_out, string(nombre_archivo).replace('.csv',''), '_meanEnv.csv'))
end

% writetable(retired_cycles,strcat(folder_out, 'RetiredCycles.csv'))