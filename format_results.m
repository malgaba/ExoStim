% Cargar datos del archivo CSV
folder = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/RESULTS/';
metrics_envelopes = readtable(strcat(folder,'Metrics_envelopes.csv'), 'Delimiter',',');

% column_names = {'OFF-pre', 'ON', 'OFF-post'};

column_names = {'OFF-pre','OFF-post'};

% Inicializar tablas para cada métrica
% results_RMS_col = array2table(zeros(0, 3), 'VariableNames', column_names);
results_RMS_col = array2table(zeros(0, 2), 'VariableNames', column_names);
results_MAV_col = array2table(zeros(0, 2), 'VariableNames', column_names);
results_MNF_col = array2table(zeros(0, 2), 'VariableNames', column_names);
results_MDF_col = array2table(zeros(0, 2), 'VariableNames', column_names);
results_RMS_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
results_MAV_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
results_MNF_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
results_MDF_col_cum = array2table(zeros(0, 2), 'VariableNames', column_names);
total_RMS = table();
total_MAV = table();
total_MNF = table();
total_MDF = table();

% Crear tabla de títulos
results_RMS = table();
results_MAV = table();
results_MNF = table();
results_MDF = table();

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
    newRow_MAV = {subject, condition, muscle};
    newRow_MNF = {subject, condition, muscle};
    newRow_MDF = {subject, condition, muscle};
    
    % Añadir filas a las tablas correspondientes
    results_RMS = [results_RMS; newRow_RMS];
    results_MAV = [results_MAV; newRow_MAV];
    results_MNF = [results_MNF; newRow_MNF];
    results_MDF = [results_MDF; newRow_MDF];

end


results_RMS = unique(results_RMS, 'rows');
results_MAV = unique(results_MAV, 'rows');
results_MNF = unique(results_MNF, 'rows');
results_MDF = unique(results_MDF, 'rows');


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
        results_MAV_col_cum = [results_MAV_col_cum; results_MAV_col];
        results_MNF_col_cum = [results_MNF_col_cum; results_MNF_col];
        results_MDF_col_cum = [results_MDF_col_cum; results_MDF_col];
    end

    % Extraer valores de las métricas
    RMS = metrics_envelopes.RMS(f);
    MAV = metrics_envelopes.MAV(f);
    MNF = metrics_envelopes.MNF(f);
    MDF = metrics_envelopes.MDF(f);

    switch run
        case 'run_01'
            col = 'OFF-pre';
        case 'run_02'
            col = 'ON';
        case 'run_03'
            col = 'OFF-post';
    end

    results_RMS_col.(col)(c) = RMS;
    results_MAV_col.(col)(c) = MAV;
    results_MNF_col.(col)(c) = MNF;
    results_MDF_col.(col)(c) = MDF;
    c = c + 1;

    if c > 14
        c = 1;
    end
    

    sub = subject;

end

results_RMS_col_cum = [results_RMS_col_cum; results_RMS_col];
results_MAV_col_cum = [results_MAV_col_cum; results_MAV_col];
results_MNF_col_cum = [results_MNF_col_cum; results_MNF_col];
results_MDF_col_cum = [results_MDF_col_cum; results_MDF_col];

total_RMS = [results_RMS, results_RMS_col_cum];
total_MAV = [results_MAV, results_MAV_col_cum];
total_MNF = [results_MNF, results_MNF_col_cum];
total_MDF = [results_MDF, results_MDF_col_cum];

% Asignar nombres de columnas a las tablas
% total_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'ON', 'OFF-post'};
total_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
total_MAV.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
total_MNF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
total_MDF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};

% Asignar nombres de columnas a las tablas
total_RMS.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
total_MAV.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
total_MNF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};
total_MDF.Properties.VariableNames = {'Subject', 'Condition', 'Muscle', 'OFF-pre', 'OFF-post'};

% Guardar tablas en un archivo Excel
writetable(total_RMS, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'RMS');
writetable(total_MAV, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'MAV');
writetable(total_MNF, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'MNF');
writetable(total_MDF, strcat(folder,'Results_envelope.xlsx'), 'Sheet', 'MDF');

writetable(total_RMS, strcat(folder,'RMS_envelope.csv'));
writetable(total_MAV, strcat(folder,'MAV_envelope.csv'));
writetable(total_MNF, strcat(folder,'MNF_envelope.csv'));
writetable(total_MDF, strcat(folder,'MDF_envelope.csv'));
