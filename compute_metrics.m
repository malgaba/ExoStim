folder_in = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/ENVELOPE/';
folder_out = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/RESULTS/';

files = dir(fullfile(folder_in,'*mat'));
muscles = {'time', 'ErSp_r', 'ReAb_r', 'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'ErSp_l', 'ReAb_l', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};
results  = table();
fs = 2048;
order = 8;

for i = 1:numel(files)
    if contains(files(i).name,'normalized')
        name_file = fullfile(folder_in, files(i).name);
        [~, nombre_archivo, ~] = fileparts(name_file)
        muscle = nombre_archivo(length(nombre_archivo)-16:length(nombre_archivo)-11);
        %Leemos los archivos csv
        data = load(name_file);
        if contains(muscle,'_l')
            data = data.individual_cycles_l;
        elseif contains(muscle,'_r')
            data = data.individual_cycles_r;
        end
    
        temp_table = table();
        temp_table = addvars(temp_table, string(nombre_archivo), 'NewVariableNames','Trial');
    
    
        rms_values = [];
        mav_values = [];
        mnf_values = [];
        mdf_values = [];
        for s=1:size(data,2)    
            % Root Mean Square
            signal = data(:,s);
            rms_value = sqrt(mean(signal.^2));
            rms_values = [rms_values, rms_value];

            % Mean Absolute Value
            mav_value = mean(signal);
            mav_values = [mav_values, mav_value];

            % Mean Frequency y median frequency (los ciclos son demasiado cortos para calcular
            % MNF y MDF con ventana)
            [Pxx, F] = pburg(signal, order, length(signal), fs);
            mnf_value = meanfreq(Pxx, F);
            mdf_value = medfreq(Pxx, F);

            mnf_values = [mnf_values mnf_value];
            mdf_values = [mdf_values mdf_value];

            % % Mean and median frequency
            % l = length(signal);
            % time_wind = fs/4;
            % st_point = 1;
            % end_point = st_point + time_wind;
            % 
            % n_wind = floor((1-st_point)/time_wind);
            % 
            % ntft_window = 2^(nextpow2(time_wind));
            % 
            % for c=1:n_wind
            %     [Pxx, F] = pburg(signal(st_point:end_point,:), order, nfft_window, fs);
            %     MNF(c,:) = meanfreq(Pxx, F);
            %     MDF(c,:) = medfreq(Pxx, F);
            % 
            %     st_point = end_point;
            %     end_point = end_point+time_wind;
            %     clearvars Pxx F
            %     if (end_point>1)
            %         end_point = 1;
            %     end
            % end
        end

        temp_table = addvars(temp_table, mean(rms_values, 'omitnan'), 'NewVariableNames','RMS');
        temp_table = addvars(temp_table, mean(mav_values, 'omitnan'), 'NewVariableNames','MAV');
        temp_table = addvars(temp_table, mean(mnf_values, 'omitnan'), 'NewVariableNames','MNF');
        temp_table = addvars(temp_table, mean(mdf_values, 'omitnan'), 'NewVariableNames','MDF');

        results = vertcat(results,temp_table);
    
    end
end

writetable(results,strcat(folder_out, 'Metrics_envelopes.csv'))