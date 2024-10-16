% Cargar datos del archivo CSV
cond = 'WITH_EXO';
folder = strcat('/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_',cond,'/PROCESSED/RESULTS/');
metrics_envelopes = readtable(strcat(folder,strcat('Metrics_envelopes_new_',cond,'.csv')), 'Delimiter',',');

if strcmp(cond,'WITH_EXO')
    column_names = {'OFF-pre', 'ON', 'OFF-post'};
else
    column_names = {'OFF-pre', 'OFF-post'};
end

% Inicializar tablas para cada métrica
if strcmp(cond,'WITH_EXO')
    results_RMS_col = array2table(zeros(0, 3), 'VariableNames', column_names);
    std_RMS_col = array2table(zeros(0, 3), 'VariableNames', column_names);
    results_MAV_col = array2table(zeros(0, 3), 'VariableNames', column_names);
    std_MAV_col = array2table(zeros(0, 3), 'VariableNames', column_names);
    % results_MNF_col = array2table(zeros(0, 3), 'VariableNames', column_names);
    % results_MDF_col = array2table(zeros(0, 3), 'VariableNames', column_names);
    results_range_col = array2table(zeros(0, 3), 'VariableNames', column_names);
    std_range_col = array2table(zeros(0, 3), 'VariableNames', column_names);
    results_iEMG_col = array2table(zeros(0, 3), 'VariableNames', column_names);
    std_iEMG_col = array2table(zeros(0, 3), 'VariableNames', column_names);
    
    results_RMS_col_cum = array2table(zeros(0, 3), 'VariableNames', column_names);
    std_RMS_col_cum = array2table(zeros(0, 3), 'VariableNames', column_names);
    results_MAV_col_cum = array2table(zeros(0, 3), 'VariableNames', column_names);
    std_MAV_col_cum = array2table(zeros(0, 3), 'VariableNames', column_names);
    % results_MNF_col_cum = array2table(zeros(0, 3), 'VariableNames', column_names);
    % results_MDF_col_cum = array2table(zeros(0, 3), 'VariableNames', column_names);
    results_range_col_cum = array2table(zeros(0, 3), 'VariableNames', column_names);
    std_range_col_cum = array2table(zeros(0, 3), 'VariableNames', column_names);
    results_iEMG_col_cum = array2table(zeros(0, 3), 'VariableNames', column_names);
    std_iEMG_col_cum = array2table(zeros(0, 3), 'VariableNames', column_names);

else
    results_RMS_col = array2table(zeros(0, 2), 'VariableNames', column_names);
    std_RMS_col = array2table(zeros(0, 2), 'VariableNames', column_names);
    results_MAV_col = array2table(zeros(0, 2), 'VariableNames', column_names);
    std_MAV_col = array2table(zeros(0, 2), 'VariableNames', column_names);
    % results_MNF_col = array2table(zeros(0, 2), 'VariableNames', column_names);
    % results_MDF_col = array2table(zeros(0, 2), 'VariableNames', column_names);
    results_range_col = array2table(zeros(0, 2), 'VariableNames', column_names);
    std_range_col = array2table(zeros(0, 2), 'VariableNames', column_names);
    results_iEMG_col = array2table(zeros(0, 2), 'VariableNames', column_names);
    std_iEMG_col = array2table(zeros(0, 2), 'VariableNames', column_names);
    
    results_RMS_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
    std_RMS_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);    
    results_MAV_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
    std_MAV_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
    % results_MNF_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
    % results_MDF_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
    results_range_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
    std_range_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
    results_iEMG_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
    std_iEMG_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
end

total_RMS = table();
total_std_RMS = table();
total_MAV = table();
total_std_MAV = table();
% total_MNF = table();
% total_MDF = table();
total_range = table();
total_iEMG = table();
total_std_range = table();
total_std_iEMG = table();

% Crear tabla de títulos
results_RMS = table();
std_RMS = table();
results_MAV = table();
std_MAV = table();
% results_MNF = table();
% results_MDF = table();
results_range = table();
results_iEMG = table();
std_range = table();
std_iEMG = table();

% Recorrer cada fila del CSV y procesar los datos
for i = 1:height(metrics_envelopes)
    % Extraer información del Trial
    trial = metrics_envelopes.Trial{i};
    tokens = split(trial, '_');
    subject = str2double(tokens{2});
    condition = strcat(tokens{3}, '_', tokens{4});
    run = strcat(tokens{5}, '_', tokens{6});
    muscle = strcat(tokens{7}, '_', tokens{8});
    
    % Crear filas para cada métrica
    newRow_RMS = {subject, condition, muscle};
    newRow_std_RMS = {subject, condition, muscle};
    newRow_MAV = {subject, condition, muscle};
    newRow_std_MAV = {subject, condition, muscle};
    % newRow_MNF = {subject, condition, muscle};
    % newRow_MDF = {subject, condition, muscle};
    newRow_range = {subject, condition, muscle};
    newRow_std_range = {subject, condition, muscle};
    newRow_iEMG = {subject, condition, muscle};
    newRow_std_iEMG = {subject, condition, muscle};

    
    % Añadir filas a las tablas correspondientes
    results_RMS = [results_RMS; newRow_RMS];
    std_RMS = [std_RMS; newRow_std_RMS];
    results_MAV = [results_MAV; newRow_MAV];
    std_MAV = [std_MAV; newRow_std_MAV];
    % results_MNF = [results_MNF; newRow_MNF];
    % results_MDF = [results_MDF; newRow_MDF];
    results_range = [results_range; newRow_range];
    std_range = [std_range; newRow_std_range];
    results_iEMG = [results_iEMG; newRow_iEMG];
    std_iEMG = [std_iEMG; newRow_std_iEMG];

end


results_RMS = unique(results_RMS, 'rows');
std_RMS = unique(std_RMS, 'rows');
results_MAV = unique(results_MAV, 'rows');
std_MAV = unique(std_MAV, 'rows');
% results_MNF = unique(results_MNF, 'rows');
% results_MDF = unique(results_MDF, 'rows');
results_range = unique(results_range, 'rows');
std_range = unique(std_range, 'rows');
results_iEMG = unique(results_iEMG, 'rows');
std_iEMG = unique(std_iEMG, 'rows');


% Crear tabla de resultados
c = 1;
sub = 1;

for f = 1:height(metrics_envelopes)
    % Extraer información del Trial
    trial = metrics_envelopes.Trial{f};
    tokens = split(trial, '_');
    subject = str2double(tokens{2});
    condition = strcat(tokens{3}, '_', tokens{4});
    run = strcat(tokens{5}, '_', tokens{6});
    muscle = strcat(tokens{7}, '_', tokens{8});

    if sub ~= subject
        results_RMS_col_cum = [results_RMS_col_cum; results_RMS_col];
        std_RMS_col_cum = [std_RMS_col_cum; std_RMS_col];
        results_MAV_col_cum = [results_MAV_col_cum; results_MAV_col];
        std_MAV_col_cum = [std_MAV_col_cum; std_MAV_col];
        % results_MNF_col_cum = [results_MNF_col_cum; results_MNF_col];
        % results_MDF_col_cum = [results_MDF_col_cum; results_MDF_col];
        results_range_col_cum = [results_range_col_cum; results_range_col];
        std_range_col_cum = [std_range_col_cum; std_range_col];
        results_iEMG_col_cum = [results_iEMG_col_cum; results_iEMG_col];
        std_iEMG_col_cum = [std_iEMG_col_cum; std_iEMG_col];
    end

    % Extraer valores de las métricas
    RMS = metrics_envelopes.RMS(f);
    RMS_STD = metrics_envelopes.STD_RMS(f);
    MAV = metrics_envelopes.MAV(f);
    MAV_STD = metrics_envelopes.STD_MAV(f);
    % MNF = metrics_envelopes.MNF(f);
    % MDF = metrics_envelopes.MDF(f);
    range = metrics_envelopes.RangeOfActivation(f);
    range_STD = metrics_envelopes.STD_RangeOfActivation(f);
    iEMG = metrics_envelopes.iEMG(f);
    iEMG_STD = metrics_envelopes.STD_iEMG(f);


    switch run
        case 'run_01'
            col = 'OFF-pre';
        case 'run_02'
            col = 'ON';
        case 'run_03'
            col = 'OFF-post';
    end

    results_RMS_col.(col)(c) = RMS;
    std_RMS_col.(col)(c) = RMS_STD;
    results_MAV_col.(col)(c) = MAV;
    std_MAV_col.(col)(c) = MAV_STD;
    % results_MNF_col.(col)(c) = MNF;
    % results_MDF_col.(col)(c) = MDF;
    results_range_col.(col)(c) = range;
    std_range_col.(col)(c) = range_STD;
    results_iEMG_col.(col)(c) = iEMG;
    std_iEMG_col.(col)(c) = iEMG_STD;

    c = c + 1;

    if c > 14
        c = 1;
    end
    

    sub = subject;

end

results_RMS_col_cum = [results_RMS_col_cum; results_RMS_col];
std_RMS_col_cum = [std_RMS_col_cum; std_RMS_col];
results_MAV_col_cum = [results_MAV_col_cum; results_MAV_col];
std_MAV_col_cum = [std_MAV_col_cum; std_MAV_col];
% results_MNF_col_cum = [results_MNF_col_cum; results_MNF_col];
% results_MDF_col_cum = [results_MDF_col_cum; results_MDF_col];
results_range_col_cum = [results_range_col_cum; results_range_col];
std_range_col_cum = [std_range_col_cum; std_range_col];
results_iEMG_col_cum = [results_iEMG_col_cum; results_iEMG_col];
std_iEMG_col_cum = [std_iEMG_col_cum; std_iEMG_col];

total_RMS = [results_RMS, results_RMS_col_cum];
total_std_RMS = [std_RMS, std_RMS_col_cum];
total_MAV = [results_MAV, results_MAV_col_cum];
total_std_MAV = [std_MAV, std_MAV_col_cum];
% total_MNF = [results_MNF, results_MNF_col_cum];
% total_MDF = [results_MDF, results_MDF_col_cum];
total_range = [results_range, results_range_col_cum];
total_std_range = [std_range, std_range_col_cum];
total_iEMG = [results_iEMG, results_iEMG_col_cum];
total_std_iEMG = [std_iEMG, std_iEMG_col_cum];


% Asignar nombres de columnas a las tablas
if strcmp(cond,'WITH_EXO')
    total_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    total_std_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    total_MAV.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    total_std_MAV.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % total_MNF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % total_MDF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    total_range.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    total_std_range.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    total_iEMG.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    total_std_iEMG.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    
    % % Asignar nombres de columnas a las tablas
    % total_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % total_std_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % total_MAV.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % total_std_MAV.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % % total_MNF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % % total_MDF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % total_range.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % total_std_range.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % total_iEMG.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
    % total_std_iEMG.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};


else
    total_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    total_std_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    total_MAV.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    total_std_MAV.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    % total_MNF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    % total_MDF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    total_range.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    total_std_range.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    total_iEMG.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    total_std_iEMG.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    
    % % Asignar nombres de columnas a las tablas
    % total_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    % total_std_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    % total_MAV.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    % % total_MNF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    % % total_MDF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    % total_range.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
    % total_iEMG.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
end

% Guardar tablas en un archivo Excel
writetable(total_RMS, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'RMS');
writetable(total_std_RMS, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'STD_RMS');
writetable(total_MAV, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'MAV');
writetable(total_std_MAV, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'STD_MAV');
% writetable(total_MNF, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'MNF');
% writetable(total_MDF, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'MDF');
writetable(total_range, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'Range');
writetable(total_std_range, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'STD_Range');
writetable(total_iEMG, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'iEMG');
writetable(total_std_iEMG, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'STD_iEMG');

% writetable(total_RMS, strcat(folder,'RMS_envelope.csv'));
% writetable(total_MAV, strcat(folder,'MAV_envelope.csv'));
% % writetable(total_MNF, strcat(folder,'MNF_envelope.csv'));
% % writetable(total_MDF, strcat(folder,'MDF_envelope.csv'));
% writetable(total_range, strcat(folder,'Range_envelope.csv'));
% writetable(total_iEMG, strcat(folder,'iEMG_envelope.csv'));

