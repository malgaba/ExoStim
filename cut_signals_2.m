%% Cut partes pochas de la señal
folder = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/FILTERED_cut/';
folder_out = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/CUT_SIGNALS/';
files = dir(fullfile(folder,'*.csv'));
muscles = {'time', 'ErSp_r', 'ReAb_r', 'GlMe_r', 'ReFe_r', 'BiFe_r', 'GaLa_r', 'TiAn_r', 'ErSp_l', 'ReAb_l', 'GlMe_l', 'ReFe_l', 'BiFe_l', 'GaLa_l', 'TiAn_l'};
cond = "01"; % 1 = without_exo, 2 = with_exo
subject = 3;
nom_file = strcat('subject_0',num2str(subject),'_cond_',cond,'_run_03.csv');

% for i = 1:numel(files)
%     name_file = fullfile(folder, files(i).name);
% Leemos el archivo csv
data = readtable(strcat(folder,nom_file));
%[~, nombre_archivo, ~] = fileparts(name_file)

for i=2:size(muscles,2)
    m = muscles{1,i};
    original_emg = data.(m);
    time = data.time;
    figure;
    plot(time, original_emg)
    xlabel('Tiempo');
    ylabel('Amplitud');
    title(m)
end

%% Once you have selected which muscle do you want to cut
m = 'BiFe_r';
original_emg = data.(m);

plot(time, original_emg)
xlabel('Tiempo');
ylabel('Amplitud');

% Habilitar la funcionalidad de zoom
zoom on;

% Esperar hasta que el usuario realice el zoom
disp('Haz zoom en la región de interés y luego presiona Enter.');
pause;

% Deshabilitar la funcionalidad de zoom
zoom off;

fprintf('Haz clic en la gráfica para seleccionar los puntos de la plantilla.\n');
fprintf('Presiona Enter cuando hayas terminado.\n');
[x, ~] = ginput(2);

if ~isempty(x)
    indice_inicio = find(time >= x(1), 1); % Índice del primer punto seleccionado
    indice_fin = find(time <= x(2), 1, 'last'); % Índice del segundo punto seleccionado

    % Divide la señal en dos partes: antes y después de la sección seleccionada
    before = table();

    after = table();
    for i=1:size(muscles,2)
        m = muscles{1,i};
        before = addvars(before,data.(m)(1:indice_inicio),'NewVariableName',m);
        after = addvars(after,data.(m)(indice_fin:end),'NewVariableName',m);
    end
    % Guarda cada parte en un archivo separado
    writetable(before, strcat(folder_out,nom_file.replace('.csv', ''), '_1.csv'));
    writetable(after, strcat(folder_out,nom_file.replace('.csv', ''), '_2.csv'));
else
    disp('No se ha seleccionado nada. La señal no se cortará.');
end