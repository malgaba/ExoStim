%% CUT SIGNALS
% Cargamos el csv con los datos filtrados
folder = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PREPROCESSED/FILTERED';
folder_out = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/FILTERED_cut/';
files = dir(fullfile(folder,'*.csv'));
muscles = {'time', 'ErSp_r', 'ReAb_r', 'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'ErSp_l', 'ReAb_l', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};
load ("/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/RESULTS_ExoStim.mat")

fieldSuffix = 'Gait';
% fieldSuffix = 'Exo';
cond = "01";
% cond = "02";

for i = 1:numel(files)
    name_file = fullfile(folder, files(i).name);
    % Leemos el archivo csv
    original_emg = readtable(name_file);
    [~, nombre_archivo, ~] = fileparts(name_file)
    
    % Leemos los index de los eventos
    
    eventsField = ['indexes', fieldSuffix];
    events = RESULTS(i).(eventsField);

    % Cut the signal and save with name
    cut_signal_1 = original_emg(events(1):events(2),:);
    writetable(cut_signal_1, strcat(folder_out,'subject_0', num2str(i),'_cond_', cond, '_run_01.csv'));

    cut_signal_2 = original_emg(events(3):events(4),:);
    writetable(cut_signal_2, strcat(folder_out,'subject_0', num2str(i),'_cond_', cond, '_run_02.csv'));

    if length(events) > 6
        cut_signal_3 = original_emg(events(5):events(6),:);
        writetable(cut_signal_3, strcat(folder_out,'subject_0', num2str(i),'_cond_', cond, '_run_03.csv'));
        cut_signal_4 = original_emg(events(7):end,:);
        writetable(cut_signal_4, strcat(folder_out,'subject_0', num2str(i),'_cond_', cond, '_run_02_2.csv'));
    elseif length(events) == 6
        cut_signal_3 = original_emg(events(5):end,:);
        writetable(cut_signal_3, strcat(folder_out,'subject_0', num2str(i),'_cond_', cond, '_run_03.csv'));
    end
end