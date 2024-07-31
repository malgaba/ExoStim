% COMPROBAR EVENTOS
subject = 2;
cond = "01"; %"" = with, "out" = without
folder = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/MEDIAN_FILTERED/';
path_events = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/TEMPLATES/';
file = strcat('subject_0',num2str(subject),'_cond_',cond,'_run_01.csv');

nom_file = strcat(folder,file)
data = readtable(nom_file);
FS = 2048;

events = readtable(strcat(path_events, file.replace('.csv',''), '_estimatedGaitEvents.csv'));


%files = dir(fullfile(folder,'*.csv'));
muscles = {'time', 'ErSp_', 'ReAb_', 'GlMe_', 'ReFe_', 'BiFe_', 'GaLa_', 'TiAn_'};

time = data.time;

for i=2:size(muscles,2)
    m = muscles{1,i};
    m_left = strcat(m,"l");
    m_right = strcat(m,"r");

    XIl = data.(m_left);

    % Xabs = abs(emg_left);
    % band = (2/FS)*1;
    % [B,A] = butter(3,band,'low');
    % XIl = filtfilt(B,A,Xabs);
    
    figure('WindowState', 'maximized');
    subplot(1,2,1);
    plot(time, XIl)
    xlabel('Tiempo');
    ylabel('Amplitud');
    title(m_left)
    for c = 1:length(events.heel_strike_l)
        line([events.heel_strike_l(c), events.heel_strike_l(c)],ylim,'Color', 'r', 'LineWidth', 2);
    end
    grid on;
    hold on;

    XIr = data.(m_right);

    % Xabs = abs(emg_right);
    % band = (2/FS)*1;
    % [B,A] = butter(3,band,'low');
    % XIr = filtfilt(B,A,Xabs);
    
    subplot(1,2,2)
    plot(time,XIr)
    xlabel('Tiempo');
    ylabel('Amplitud');
    title(m_right)
    for c = 1:length(events.heel_strike_r)
        line([events.heel_strike_r(c), events.heel_strike_r(c)],ylim,'Color', 'r', 'LineWidth', 2);
    end
    grid on;
    % 
    % if contains(m, '_l')
    %     for c = 1:length(events.heel_strike_l)
    %         line([events.heel_strike_l(c), events.heel_strike_l(c)],ylim,'Color', 'r', 'LineWidth', 2);
    %     end
    % elseif contains(m, '_r')
    %     for c = 1:length(events.heel_strike_r)
    %         line([events.heel_strike_r(c), events.heel_strike_r(c)],ylim,'Color', 'r', 'LineWidth', 2);
    %     end
    % end

    linkaxes([subplot(1, 2, 1), subplot(1, 2, 2)], 'x');
    linkaxes([subplot(1, 2, 1), subplot(1, 2, 2)], 'y');

    legend('EMG', 'Events');
    hold off;

end

%% Add events manually
close all
muscle = 'TiAn';
subplot(2,1,1)
h1 = plot(time,data.(strcat(muscle,'_l')));
for c = 1:length(events.heel_strike_l)
        line([events.heel_strike_l(c), events.heel_strike_l(c)],ylim,'Color', 'r', 'LineWidth', 2);
end

subplot(2,1,2)
h2 = plot(time,data.(strcat(muscle,'_r')));
for c = 1:length(events.heel_strike_r)
        line([events.heel_strike_r(c), events.heel_strike_r(c)],ylim,'Color', 'r', 'LineWidth', 2);
end

fprintf('Seleccione los eventos que falten y pulse enter cuando termine')
[x, ~] = ginput();

list_l = events.heel_strike_l;
list_r = events.heel_strike_r;

if ismember(gca, findall(gcf, 'Type', 'axes'))
    ax = gca;

    if ax == h1.Parent
        % Añadir punto al subplot 1
        for i=1:length(x)
            list_l = [list_l; x(i)];
            list_l = sort(list_l);
        end
    elseif ax == h2.Parent
        for i=1:length(x)
            list_r = [list_r; x(i)];
            list_r = sort(list_r);
        end
    end
end

len1 = height(list_l);
len2 = height(list_r);
maxLen = max(len1, len2);

% Rellenar con NaN para igualar las longitudes
if len1 < maxLen
    list_l(len1+1:maxLen,1) = NaN;
    
elseif len2 < maxLen
    list_r(len2+1:maxLen,1) = NaN;
end

new_events = table(list_l, list_r,'VariableNames',{'heel_strike_l','heel_strike_r'});

writetable(new_events, strcat(path_events, file.replace('.csv',''), '_gaitEvents.csv'));

close all
