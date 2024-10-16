cond = 'WITH_EXO';
folder_in = strcat('/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_',cond,'/PROCESSED/ENVELOPE');
folder_out = strcat('/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_',cond,'/PROCESSED/RESULTS/');


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
        % mnf_values = [];
        % mdf_values = [];
        activation_ranges = [];

        for s=1:size(data,2)    
            % Root Mean Square
            signal = data(:,s);
            % signal = signal*1000; %Pasar a mV
            % rms_value = sqrt(mean(signal.^2));

            rms_value = rms(signal);
            rms_values = [rms_values, rms_value];

            


            % Mean Absolute Value
            mav_value = mean(abs(signal));
            mav_values = [mav_values, mav_value];

          
            % Mean Frequency y median frequency (los ciclos son demasiado cortos para calcular
            % MNF y MDF con ventana)
            % [Pxx, F] = pburg(signal, order, length(signal), fs);
            % mnf_value = meanfreq(Pxx, F);
            % mdf_value = medfreq(Pxx, F);
            % 
            % mnf_values = [mnf_values, mnf_value];
            % mdf_values = [mdf_values, mdf_value];

            activation_range = max(signal) - min(signal);
            activation_ranges = [activation_ranges, activation_range];

           
        end
        
        writematrix(rms_values', strcat(folder_out,strrep(files(i).name,'normalized.mat','RMS.csv')));
        writematrix(mav_values', strcat(folder_out,strrep(files(i).name,'normalized.mat','MAV.csv')));
        writematrix(activation_ranges', strcat(folder_out,strrep(files(i).name,'normalized.mat','RANGE.csv')));


       
        temp_table = addvars(temp_table, mean(rms_values, 'omitnan'), 'NewVariableNames','RMS');
        temp_table = addvars(temp_table, std(rms_values, 'omitnan'), 'NewVariableNames','STD_RMS');
        temp_table = addvars(temp_table, mean(mav_values, 'omitnan'), 'NewVariableNames','MAV');
        temp_table = addvars(temp_table, std(mav_values, 'omitnan'), 'NewVariableNames','STD_MAV');
        % temp_table = addvars(temp_table, mean(mnf_values, 'omitnan'), 'NewVariableNames','MNF');
        % temp_table = addvars(temp_table, mean(mdf_values, 'omitnan'), 'NewVariableNames','MDF');
        temp_table = addvars(temp_table, mean(activation_ranges, 'omitnan'), 'NewVariableNames','Range of activation');
        temp_table = addvars(temp_table, std(activation_ranges, 'omitnan'), 'NewVariableNames','STD_Range of activation');
    
        results = vertcat(results,temp_table);
    end
    
end


if strcmp(cond, 'WITHOUT_EXO')
    folder_in_iEMG = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/ENVELOPE';
elseif strcmp(cond, 'WITH_EXO')
    folder_in_iEMG = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITH_EXO/PROCESSED/ENVELOPE/NORMALIZED';
end

files_iEMG = dir(fullfile(folder_in_iEMG,'*mat'));

n_samples = 200; %Número de muestras de cada segmento
t = (0:n_samples-1) / fs;

results_2 = table();

for i = 1:numel(files_iEMG)
    % Porque el iEMG no se saca del envelope, sino de la señal normal
    if ~contains(files_iEMG(i).name,'normalized')
        name_file = fullfile(folder_in_iEMG, files_iEMG(i).name);
        [~, nombre_archivo, ~] = fileparts(name_file)
        muscle = nombre_archivo(length(nombre_archivo)-5:length(nombre_archivo));
        %Leemos los archivos csv
        data = load(name_file);

        if strcmp(cond, 'WITHOUT_EXO')
            if contains(muscle,'_l')
                data = data.segments_emg_l;
            elseif contains(muscle,'_r')
                data = data.segments_emg_r;
            end

        elseif strcmp(cond, 'WITH_EXO')
            data = data.data;
        end

        data = data * 1000; % Pasamos a mV

        % EMG data normalized by the average of its peaks from all cycles
        max_individual_cycles = [];
        for j = 1:size(data,2)
            max_individual_cycles(j) = max(data(:,j));
        end

        media = mean(max_individual_cycles);
        data = data/media;

        temp_table_2 = table();
        % temp_table = addvars(temp_table, string(nombre_archivo), 'NewVariableNames','Trial');

        iEMGs = [];
        for s=1:size(data,2)    
            signal = data(:,s);
            rectified_signal = abs(signal);
            % iEMG = sum(abs(signal))/length(signal); %mV/s
            iEMG = trapz(t, rectified_signal);
            
            iEMGs = [iEMGs iEMG];

        end

        writematrix(iEMGs', strcat(folder_out,strrep(files(i).name,'.mat','_iEMG.csv')));

        temp_table_2 = addvars(temp_table_2, mean(iEMGs, 'omitnan'), 'NewVariableNames','iEMG');
        temp_table_2 = addvars(temp_table_2, std(iEMGs, 'omitnan'), 'NewVariableNames','STD_iEMG');
        results_2 = vertcat(results_2,temp_table_2);

    end
end

final_results = horzcat(results,results_2);

writetable(final_results,strcat(folder_out, strcat('Metrics_envelopes_new_',cond,'.csv')))