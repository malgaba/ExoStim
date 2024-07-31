%% Concat signals
folder_in = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/MEDIAN_FILTERED/';
folder_events = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/TEMPLATES/';
folder_out = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/CONCAT_SIGNALS/';
muscles = {'time', 'ErSp_r', 'ReAb_r', 'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'ErSp_l', 'ReAb_l', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};
files_1 = dir(fullfile(folder_in,'*_1.csv'));
files_2 = dir(fullfile(folder_in,'*_2.csv'));
files = [files_1;files_2];

for i = 1:numel(files)
    name_file = fullfile(folder_in, files(i).name);
    [~, nombre_archivo, ~] = fileparts(name_file);
    nombre_archivo = extractBefore(string(nombre_archivo), strlength(nombre_archivo) - 1)

    data_1 = readtable(strcat(folder_in, nombre_archivo, '_1.csv'));
    data_2 = readtable(strcat(folder_in, nombre_archivo, '_2.csv'));
    data = vertcat(data_1(:,1:15),data_2(:,1:15));
    
    

    events_1 = readtable(strcat(folder_events, nombre_archivo, '_1_estimatedGaitEvents.csv'));
    events_2 = readtable(strcat(folder_events, nombre_archivo, '_2_estimatedGaitEvents.csv'));
    events = vertcat(events_1, events_2);
    for s=2:size(muscles,2)
        m = muscles{1,s};
        figure;
        plot(data.time, data.(m))
        if contains(m,'l')
            for c = 1:length(events.heel_strike_l)
                line([events.heel_strike_l(c), events.heel_strike_l(c)],ylim,'Color', 'r', 'LineWidth', 2);
            end
        else
            for c = 1:length(events.heel_strike_r)
                line([events.heel_strike_r(c), events.heel_strike_r(c)],ylim,'Color', 'r', 'LineWidth', 2);
            end
        end

        title(m)
    end
    

    disp('Pulse Enter para continuar...');
    pause_input = input('','s');
    close all
    % 
    writetable(data,strcat(folder_out,nombre_archivo,'.csv'))
    writetable(events,strcat(folder_out,nombre_archivo,'_gaitEvents.csv'))

end




