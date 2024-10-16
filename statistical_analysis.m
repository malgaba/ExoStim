% STATISTICAL ANALYSIS WITHOUT EXO
cond = 'WITH_EXO';
folder_in = strcat('/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_',cond,'/PROCESSED/RESULTS/');
folder_out = strcat('/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_',cond,'/PROCESSED/STATISTICAL_ANALYSIS/');

%RMS SUBJECT_01 WITHOUT EXO
muscles_01 = {'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};

results = table([], [], [], [], [], [], [], [], [], 'VariableNames', {'Trial', 'RMS', 'h_RMS', 'MAV', 'h_MAV', 'iEMG', 'h_iEMG', 'Range', 'h_Range'});


for i = 1:length(muscles_01)
    name_muscle = muscles_01(i);
    file_pre = strcat(folder_in,'subject_01_cond_01_run_01_',name_muscle{1,1}); 
    file_post = strcat(folder_in,'subject_01_cond_01_run_03_',name_muscle{1,1}); 
    % RMS
    % x = RMS OFF PRE
    % y = RMS OFF POST
    x_RMS = load(strcat(file_pre,'_RMS.csv'));
    y_RMS = load(strcat(file_post,'_RMS.csv'));
    [p_RMS, h_RMS] = ranksum(x_RMS,y_RMS);

    x_MAV = load(strcat(file_pre,'_MAV.csv'));
    y_MAV = load(strcat(file_post,'_MAV.csv'));
    [p_MAV, h_MAV] = ranksum(x_MAV,y_MAV);

    x_iEMG = load(strcat(file_pre,'_iEMG.csv'));
    y_iEMG = load(strcat(file_post,'_iEMG.csv'));
    [p_iEMG, h_iEMG] = ranksum(x_iEMG,y_iEMG);

    x_range = load(strcat(file_pre,'_RANGE.csv'));
    y_range = load(strcat(file_post,'_RANGE.csv'));
    [p_range, h_range] = ranksum(x_range,y_range);


    [~, name_file_pre, ~] = fileparts(file_pre);
    [~, name_file_post, ~] = fileparts(file_post);
    new_row = {strcat(name_file_pre,'-VS-', name_file_post), p_RMS, h_RMS, p_MAV, h_MAV, p_iEMG, h_iEMG, p_range, h_range};
    results = [results; new_row];  % Añadir nueva fila a la tabla
end

muscles_03 = {'GlMe_r', 'ReFe_r', 'BiFe_r', 'TiAn_r', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'TiAn_l'};

for i = 1:length(muscles_03)
    name_muscle = muscles_03(i);
    file_pre = strcat(folder_in,'subject_03_cond_01_run_01_',name_muscle{1,1}); 
    file_post = strcat(folder_in,'subject_03_cond_01_run_03_',name_muscle{1,1}); 
    % RMS
    % x = RMS OFF PRE
    % y = RMS OFF POST
    x_RMS = load(strcat(file_pre,'_RMS.csv'));
    y_RMS = load(strcat(file_post,'_RMS.csv'));
    [p_RMS, h_RMS] = ranksum(x_RMS,y_RMS);

    x_MAV = load(strcat(file_pre,'_MAV.csv'));
    y_MAV = load(strcat(file_post,'_MAV.csv'));
    [p_MAV, h_MAV] = ranksum(x_MAV,y_MAV);

    x_iEMG = load(strcat(file_pre,'_iEMG.csv'));
    y_iEMG = load(strcat(file_post,'_iEMG.csv'));
    [p_iEMG, h_iEMG] = ranksum(x_iEMG,y_iEMG);

    x_range = load(strcat(file_pre,'_RANGE.csv'));
    y_range = load(strcat(file_post,'_RANGE.csv'));
    [p_range, h_range] = ranksum(x_range,y_range);


    [~, name_file_pre, ~] = fileparts(file_pre);
    [~, name_file_post, ~] = fileparts(file_post);
    new_row = {strcat(name_file_pre,'-VS-', name_file_post), p_RMS, h_RMS, p_MAV, h_MAV, p_iEMG, h_iEMG, p_range, h_range};
    results = [results; new_row];  % Añadir nueva fila a la tabla
end

muscles_04 = {'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};

for i = 1:length(muscles_04)
    name_muscle = muscles_04(i);
    file_pre = strcat(folder_in,'subject_04_cond_01_run_01_',name_muscle{1,1}); 
    file_post = strcat(folder_in,'subject_04_cond_01_run_03_',name_muscle{1,1}); 
    % RMS
    % x = RMS OFF PRE
    % y = RMS OFF POST
    x_RMS = load(strcat(file_pre,'_RMS.csv'));
    y_RMS = load(strcat(file_post,'_RMS.csv'));
    [p_RMS, h_RMS] = ranksum(x_RMS,y_RMS);

    x_MAV = load(strcat(file_pre,'_MAV.csv'));
    y_MAV = load(strcat(file_post,'_MAV.csv'));
    [p_MAV, h_MAV] = ranksum(x_MAV,y_MAV);

    x_iEMG = load(strcat(file_pre,'_iEMG.csv'));
    y_iEMG = load(strcat(file_post,'_iEMG.csv'));
    [p_iEMG, h_iEMG] = ranksum(x_iEMG,y_iEMG);

    x_range = load(strcat(file_pre,'_RANGE.csv'));
    y_range = load(strcat(file_post,'_RANGE.csv'));
    [p_range, h_range] = ranksum(x_range,y_range);


    [~, name_file_pre, ~] = fileparts(file_pre);
    [~, name_file_post, ~] = fileparts(file_post);
    new_row = {strcat(name_file_pre,'-VS-', name_file_post), p_RMS, h_RMS, p_MAV, h_MAV, p_iEMG, h_iEMG, p_range, h_range};
    results = [results; new_row];  % Añadir nueva fila a la tabla  
end

writetable(results,strcat(folder_out, 'comparison_', cond, '.csv'))

%% STATISTICAL ANALYSIS WITH EXO
cond = 'WITH_EXO';
folder_in = strcat('/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_',cond,'/PROCESSED/RESULTS/');
folder_out = strcat('/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_',cond,'/PROCESSED/STATISTICAL_ANALYSIS/');

%RMS SUBJECT_01 WITHOUT EXO
muscles_01 = {'ReFe_r', 'BiFe_r', 'GaLa_r', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};

results = table([], [], [], [], [], [], [], [], [], 'VariableNames', {'Trial', 'RMS', 'h_RMS', 'MAV', 'h_MAV', 'iEMG', 'h_iEMG', 'Range', 'h_Range'});


for i = 1:length(muscles_01)
    name_muscle = muscles_01(i);
    file_pre = strcat(folder_in,'subject_01_cond_02_run_01_',name_muscle{1,1}); 
    file_post = strcat(folder_in,'subject_01_cond_02_run_03_',name_muscle{1,1}); 
    % RMS
    % x = RMS OFF PRE
    % y = RMS OFF POST
    x_RMS = load(strcat(file_pre,'_RMS.csv'));
    y_RMS = load(strcat(file_post,'_RMS.csv'));
    [p_RMS, h_RMS] = ranksum(x_RMS,y_RMS);

    x_MAV = load(strcat(file_pre,'_MAV.csv'));
    y_MAV = load(strcat(file_post,'_MAV.csv'));
    [p_MAV, h_MAV] = ranksum(x_MAV,y_MAV);

    x_iEMG = load(strcat(file_pre,'_iEMG.csv'));
    y_iEMG = load(strcat(file_post,'_iEMG.csv'));
    [p_iEMG, h_iEMG] = ranksum(x_iEMG,y_iEMG);

    x_range = load(strcat(file_pre,'_RANGE.csv'));
    y_range = load(strcat(file_post,'_RANGE.csv'));
    [p_range, h_range] = ranksum(x_range,y_range);


    [~, name_file_pre, ~] = fileparts(file_pre);
    [~, name_file_post, ~] = fileparts(file_post);
    new_row = {strcat(name_file_pre,'-VS-', name_file_post), p_RMS, h_RMS, p_MAV, h_MAV, p_iEMG, h_iEMG, p_range, h_range};
    results = [results; new_row];  % Añadir nueva fila a la tabla
end

muscles_03 = {'ReFe_r', 'BiFe_r', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'TiAn_l'};

for i = 1:length(muscles_03)
    name_muscle = muscles_03(i);
    file_pre = strcat(folder_in,'subject_03_cond_02_run_01_',name_muscle{1,1}); 
    file_post = strcat(folder_in,'subject_03_cond_02_run_03_',name_muscle{1,1}); 
    % RMS
    % x = RMS OFF PRE
    % y = RMS OFF POST
    x_RMS = load(strcat(file_pre,'_RMS.csv'));
    y_RMS = load(strcat(file_post,'_RMS.csv'));
    [p_RMS, h_RMS] = ranksum(x_RMS,y_RMS);

    x_MAV = load(strcat(file_pre,'_MAV.csv'));
    y_MAV = load(strcat(file_post,'_MAV.csv'));
    [p_MAV, h_MAV] = ranksum(x_MAV,y_MAV);

    x_iEMG = load(strcat(file_pre,'_iEMG.csv'));
    y_iEMG = load(strcat(file_post,'_iEMG.csv'));
    [p_iEMG, h_iEMG] = ranksum(x_iEMG,y_iEMG);

    x_range = load(strcat(file_pre,'_RANGE.csv'));
    y_range = load(strcat(file_post,'_RANGE.csv'));
    [p_range, h_range] = ranksum(x_range,y_range);


    [~, name_file_pre, ~] = fileparts(file_pre);
    [~, name_file_post, ~] = fileparts(file_post);
    new_row = {strcat(name_file_pre,'-VS-', name_file_post), p_RMS, h_RMS, p_MAV, h_MAV, p_iEMG, h_iEMG, p_range, h_range};
    results = [results; new_row];  % Añadir nueva fila a la tabla
end

muscles_04 = {'GlMe_r', 'BiFe_r', 'GaLa_r', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};

for i = 1:length(muscles_04)
    name_muscle = muscles_04(i);
    file_pre = strcat(folder_in,'subject_04_cond_02_run_01_',name_muscle{1,1}); 
    file_post = strcat(folder_in,'subject_04_cond_02_run_03_',name_muscle{1,1}); 
    % RMS
    % x = RMS OFF PRE
    % y = RMS OFF POST
    x_RMS = load(strcat(file_pre,'_RMS.csv'));
    y_RMS = load(strcat(file_post,'_RMS.csv'));
    [p_RMS, h_RMS] = ranksum(x_RMS,y_RMS);

    x_MAV = load(strcat(file_pre,'_MAV.csv'));
    y_MAV = load(strcat(file_post,'_MAV.csv'));
    [p_MAV, h_MAV] = ranksum(x_MAV,y_MAV);

    x_iEMG = load(strcat(file_pre,'_iEMG.csv'));
    y_iEMG = load(strcat(file_post,'_iEMG.csv'));
    [p_iEMG, h_iEMG] = ranksum(x_iEMG,y_iEMG);

    x_range = load(strcat(file_pre,'_RANGE.csv'));
    y_range = load(strcat(file_post,'_RANGE.csv'));
    [p_range, h_range] = ranksum(x_range,y_range);


    [~, name_file_pre, ~] = fileparts(file_pre);
    [~, name_file_post, ~] = fileparts(file_post);
    new_row = {strcat(name_file_pre,'-VS-', name_file_post), p_RMS, h_RMS, p_MAV, h_MAV, p_iEMG, h_iEMG, p_range, h_range};
    results = [results; new_row];  % Añadir nueva fila a la tabla  
end

writetable(results,strcat(folder_out, 'comparison_', cond, '.csv'))
