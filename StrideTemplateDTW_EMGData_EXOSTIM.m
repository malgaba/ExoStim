%% This code has been written in order to take data collected from free-
%  living axivity sensor accelerometer data collected as part of the
%  greater Ortho robotics study and apply multi-dimensional dynamic time
%  warping to the data in order to detect individual steps taken as part of
%  steady-state straight walking bouts
%
%  The MD-DTW method requires templates to be created of each input signal,
%  and as such walking bouts were initially identified from the data for 
%  22 participants and manually segmented in order to identify a single 
%  stride using the vertical acceleration signal from the shank using a 
%  similar method to that proposed by (Jasiewicz et al. Gait & Posture, 
%  2006.) based on the foot acceleration curves and from there the points 
%  identified were corrected based on local maxima detection to ensure 
%  proper segmentation.
%  This process was updated in March of 2024 following comparison with
%  optical motion capture data in which the proper timing of IC (initial
%  contact) was determined. Based on this update, new templates will be
%  created to be used for stride identification moving forward. 
%
%  These templates must be time and scale normalized such that they
%  ultimately can be used to identify similar patterns in continuous
%  free-living data. 
%
%  This code will take these templates and determine the ultimate
%  requirements for identifying strides from free-living data (hopefully
%  a method which does not require steps to be manually detected for each
%  participant). Then strides will be pulled from the entirety of the
%  collection period to ultimately examine stride variability across a 7
%  day period
%
%Created by: Annemarie Laudanski            on October 26, 2023
%Updated by: Annemarie Laudanski            on May 27, 2024
%Updated by: Marina Algaba Vidoy            on May 31, 2024


clear;clc;%close all

%% Initialization
Datapath = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/MEDIAN_FILTERED/';
total_events = table();

% Clinical Sensor Data
% Participants_ALL = {'02','05','07','08','10','11','12','13','14','15','16','18', ...
%     '19','20','21','22','23','24','25','26','32','34'}; 
% Data Collection Session
% Session_ALL = {'1','1','1','1','1','1','1','1','1','1','1','1','1','1','1',...
%     '1','1','1','1','1','1','1'};
% % Day of collection from which steady-state walking bout was identified
% Days_ALL = {'3','4','4','4','4','4','4','4','3','3','4','3','3','3','3','3',...
%     '3','3','4','4','4','4'};
% % Surgical Leg
% Leg_ALL = {'L','R','R','R','L','R','L','L','R','R','L','L', ...
%     'R','L','R','L','L','R','R','L','R','R'};
% 
% % Walking Bout Subgroup from which templates were created
% 
% Participants = {'02','05','07','08','10','11','12','13','14','15','16','18', ...
%     '19','20','21','22','23','24','25','26','32','34'}; 
% Session = {'1','1','1','1','1','1','1','1','1','1','1','1','1','1','1',...
%     '1','1','1','1','1','1','1'};
% 
% % Days = {'4','6','7','6','5','7','6','4','6','2','2','7','4','7','7','2',...
% %     '7','2','6','7','6'};
% 
% Leg = {'L','R','R','R','L','R','L','L','R','R','L','L', ...
%     'R','L','R','L','L','R','R','L','R','R'};

% WalkingBoutStart = [5992350,4465740,4724240,4584570,7645420,5622300,...
%     5999970,4732020,5171600,3390090,3689270,3867940,5888690,3249380,...
%     4390600,5499900,4547280,3230700,5363060,2770380,3942420];
% 
% WalkingBoutEnd = [6020450,4494670,4762070,4611440,7675000,5654390,6040160,...
%     4761130,5202900,3425520,3713570,3951060,5913470,3280940,4423980,...
%     5531740,4577060,3264240,5405500,2862000,3968750];

% TrlsForTemps = [1,3,4;2,6,7;1,NaN,NaN;3,4,5;4,5,6;2,7,8;...
%     1,5,6;5,8,10;2,NaN,NaN;5,6,10;1,4,5;2,NaN,NaN;1,5,7;4,8,NaN;...
%     1,NaN,NaN;2,5,9;2,3,4;2,4,7;2,5,7;2,3,11;11,12,14;3,4,5];

% Filter Parameters
% CoF = 10;  % cutoff frequency
FS = 2048;
% [b,a] = butter(2,CoF/(FS/2),'low');
% bpl = 1; %lower bandpass frequency (Hz)
% bpu = 15; %upper bandpass frequency (Hz)

% CoFt = 15;  % cutoff frequency
% [bt,at] = butter(2,CoFt/(FS/2),'low');

% Threshold value for minima detection
ThreshVals = 60;
%NOTE - Participants 12 and 25 seem to have minima sitting much higher than
%       the others, 12 would need a min of at least 16, and 25 needs 10 at
%       least... so for now we will set the min at 16 to see, otherwise it 
%       seems like a minimum below 5 would work for everyone else.

%% Processing - first save all stride templates created using the: 
%  "ManuallyDetectHeelStrike_Updated2024.m" code
path_templates = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/TEMPLATES/';
name_file = 'subject_02_cond_01_run_01';
templates = dir(fullfile(path_templates,strcat(name_file, '*.csv')));
clear x_start x_end

SubDataPath = [Datapath '/'];
GaitDataFile = [SubDataPath strcat(name_file, '.csv')];
muscle_template = 'GaLa_l';

data = readtable(GaitDataFile);

% for i = 1:length(templates)
    % strideNums = double.empty(length(Participants_ALL),0);
for sub = 1%:length(Participants_ALL)
    tic
    % First read in the data and the segmentation times
    
%     StanceDataFile = [SubDataPath 'QtStd.csv'];
    % HeelStrikeEventsFile = [path_templates templates(i).name];
    HeelStrikeEventsFile = [path_templates strcat(name_file,'_template_',muscle_template,'.csv')];

%     HeelStrikeEventsFileb = [SubDataPath 'HSIndex_Vert_ISTA24e.csv'];%
    
%     raw_GAcc = readtable(GaitDataFile);

    % muscle_template = templates(i).name(length(templates(i).name)-9:length(templates(i).name)-4);
    
    XRaw = data.(muscle_template);
    Xabs = abs(XRaw);
    band = (2/FS)*1;
    [B,A] = butter(3,band,'low');
    X = filtfilt(B,A,Xabs);
    ShankAcc = X;
    clear raw_GAcc
    % raw_SAcc = readcell(StanceDataFile);
    % ShankAccStnd = cell2mat(raw_SAcc(:,7:9));

%     if isfile(HeelStrikeEventsFileb)
%         raw_GE = readcell(HeelStrikeEventsFileb);
%     else
    raw_GE = readtable(HeelStrikeEventsFile);
%     end
    GaitEvents = table2array(raw_GE);

    % Now that we have the data imported, we can use the GaitEvents to 
    % segment the ShankAcc (the walking bout) into single strides for each
    % of the 3 axes.
    % save the surgical leg along with the participant data for future uses
    % ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).Side = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Segment and save the Stride Templates
    
    count = 1;

    for strd = 1:size(GaitEvents,1)-1
            % We need to make sure that any longer intervals which would 
            % contain more than one step is taken out. We can do this by 
            % calculating the mean and std of GaitEvents and removing any 
            % interval where the spread is bigger than Mean+SD
%         if (GaitEvents(strd,2)-GaitEvents(strd,1)) > mean(GaitEvents(:,2)-GaitEvents(:,1))+...
%                 std(GaitEvents(:,2)-GaitEvents(:,1))
%             %skip this interval
%         elseif (GaitEvents(strd,2)-GaitEvents(strd,1)) < 25 
%             %skip this interval
%         else
            %             if sub == 3
            %                 %for sub ROB007 the segmentation was done from frame 500 so
            %                 %500 needs to be added to each of the gait event indices
            %                 sgmntx = ShankAcc(GaitEvents(strd)+500:GaitEvents(strd+1)+500,1);
            %                 sgmnty = ShankAcc(GaitEvents(strd)+500:GaitEvents(strd+1)+500,2);
            %                 sgmntz = ShankAcc(GaitEvents(strd)+500:GaitEvents(strd+1)+500,3);
            %             else
            sgmntx = ShankAcc(GaitEvents(strd):GaitEvents(strd+1),1);
            %             plot(sgmntx)
            %             hold on
%             sgmnty = ShankAcc(GaitEvents(strd,1):GaitEvents(strd,2),2);
%             sgmntz = ShankAcc(GaitEvents(strd,1):GaitEvents(strd,2),3);
            %             end
            % Interpolate the data to 201 points.
            new_yi_length = length(sgmntx);
            xi = 1:new_yi_length;
            ai = linspace(1,new_yi_length,200);
            for ji = 1:size(ai,2)
                t1i = ai(ji);
                ShankAccTemplateStrides.Strides.x(ji,count) = interp1(xi',sgmntx(1:new_yi_length),t1i);
%                 ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).Strides.y(ji,count) = interp1(xi,sgmnty(1:new_yi_length),t1i);
%                 ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).Strides.z(ji,count) = interp1(xi,sgmntz(1:new_yi_length),t1i);
            end
            count = count+1;
%         end
    end
    toc
    clear ShankAcc ShankAccStnd GaitEvents sgmntx sgmnty sgmntz new_yi_length xi ai 
end
% end
%% Plot Templates and select those most appropriate for use moving forward
for n = 3
    plotTo = size(ShankAccTemplateStrides.Strides.x,2);
    trials = [1:plotTo]; %[4,8]; %
    figure; hold on; 
    plot(ShankAccTemplateStrides.Strides.x(:,trials))
%     plot(mean(ShankAccTemplateStrides.(['ROB0' Participants_ALL{n}]).(['Session' Session_ALL{n}]).Strides.x(:,trials),2),'k','LineWidth',3)
end
% set(gcf, 'Position', get(0, 'Screensize'));
%NOTES 
% 13467, 26789, 1--,345,45679,278,156711,458910(5810),2--,
% 5610,145,2--,157,48-,1--,259,234,247,257,2311,111214,345
% (3,9,12,14,15(ends dont align), are ok but not for mean curves)
%% Take a look at the newly formed templates for each participant
n_template = 10;
close all
% for sub = 1:length(Participants_ALL)
%     % figure; hold on;
%     for temp = 1:4
%         ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).BPStrides.x(:,temp) = ...
%             bandpass(ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).Strides.x(:,TrlsForTemps(sub,temp)),[bpl bpu],FS);
%         ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).BPStrides.y(:,temp) = ...
%             bandpass(ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).Strides.y(:,TrlsForTemps(sub,temp)),[bpl bpu],FS);
%         ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).BPStrides.z(:,temp) = ...
%             bandpass(ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).Strides.z(:,TrlsForTemps(sub,temp)),[bpl bpu],FS);
% 
%         % plot(bandpass(ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).Strides.x(:,TrlsForTemps(sub,temp)),[bpl bpu],FS)/8)
%     end
%     ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).MeanBPStrides.x(:,1) = ...
%         mean(ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).BPStrides.x,2);
%     ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).MeanBPStrides.y(:,1) = ...
%         mean(ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).BPStrides.y,2);
%     ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).MeanBPStrides.z(:,1) = ...
%         mean(ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).BPStrides.z,2);
% end
% for sub = 1:22
%     figure; hold on;
%     plot(ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).BPStrides.x);
%     plot(ShankAccTemplateStrides.(['ROB0' Participants_ALL{sub}]).(['Session' Session_ALL{sub}]).MeanBPStrides.x,'k','LineWidth',3)
% end

%% Now process the walking bouts- eventually we want to extend this to data from the
%  Entire free-living period however start with just the walking bouts so I
%  can compare outcomes with those previously achieved for OARSI
tic
low = 1;
FullDataPath = Datapath;
for sub = 1 %:length(Participants) %3:length(Participants)
    for FLDay = 1
        % Load in the Walking Bout Data
    
        % XI is the data we used from which to identify the templates but it's
        % possible that a combination of all axes will be useful for the
        % identification from free-living 
    
%         FullSubDataPath = [FullDataPath 'ROB0' (Participants{sub}) '\Session' (Session{sub}) '\'];
        FullDayGaitDataFile = [FullDataPath strcat(name_file, '.csv')]; %' (Days{sub}) '
        GaitDataFile = [SubDataPath strcat(name_file, '.csv')];
        FDShankAll = readtable(FullDayGaitDataFile);
        FDShankAcc = FDShankAll.(muscle_template);
        time = FDShankAll.time;

        subplot(2,1,1)
        plot(time,FDShankAcc)
        title('Filtered signal')
        grid on;
       


        %ENVELOPE the data
        Xabs = abs(FDShankAcc);
        band = (2/FS)*low;
        [B,A] = butter(3,band,'low');
        XIc = filtfilt(B,A,Xabs);
        
        % Try band-pass filtering since we applied this same filter to the data
        % prior to selecting the strides for templating.
    
%         XIc = bandpass(FDShankAcc(:,1),[bpl bpu],FS);
%         YIc = bandpass(FDShankAcc(:,2),[bpl bpu],FS);
%         ZIc = bandpass(FDShankAcc(:,3),[bpl bpu],FS);
        clear FDShankAcc
        % toc
    
        % We will save the templates as XJ,YJ,ZJ
        % March 20, 2024 - try running the calculation for each walking bout
        % using the mean template generated from each participant (i.e. 22
        % permutations of calculations) - this way we can see if there is an
        % optimal subset of templates which would generate the most reliable
        % stride identification.
        % NOTE THIS WILL BE POST ISTA ABSTRACT!!!
    
    %     for tmpltCount = 1:size
%         TrialChoice = TrlsForTemps(sub,:);
%         TrialChoice = TrialChoice(~isnan(TrialChoice));
        XJ = ShankAccTemplateStrides.Strides.x(:,n_template);
%         YJ = filtfilt(bt,at,mean(ShankAccTemplateStrides.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).Strides.y(:,TrlsForTemps(sub,:)),2,"omitnan"));
%         ZJ = filtfilt(bt,at,mean(ShankAccTemplateStrides.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).Strides.z(:,TrlsForTemps(sub,:)),2,"omitnan"));
%         clear TrialChoice
        % figure; plot(ShankAccTemplateStrides.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).Strides.x(:,1:6))
        % hold on; plot(XJ,'k','LineWidth',2)
        XJ = normalize(XJ,'range',[-1 1]);
%         YJ = normalize(YJ,'range',[-1 1]);
%         ZJ = normalize(ZJ,'range',[-1 1]);
        % figure; plot(XJ)
    
        % Window the data and process each window to find the distance
        % functions for each scaled segment. Let's start by trying a 30 second
        % window with a 15 second overlap. Create "ind" a variable which
        % contains the index of each frame contained in each window.

        % March 23 - bump this up to try a 2 minute window with 1 min overlap?
%         blksze = 120*FS; 
%         ind = bsxfun(@plus, (1 : blksze), (0:blksze/2:numel(XIc) - (blksze/2)).');
%         if ind(end) > numel(XIc)
%             ind(end,:) = [];
%             ind(end+1,1:blksze) = [ind(end)+1:1:numel(XIc), NaN(1,(blksze-(numel(XIc)-(ind(end)))),'single')];
%         end
    
        % Now loop by the number of rows in ind to calculate the distance
        % matrix from the normalized data for each window. 
        for wnd = 1
            % tic
            n_samples = linspace(1,size(XIc,1),size(XIc,1))';
            XI = normalize(XIc,'range',[-1 1]);
            % toc
            subplot(2,1,2)
            plot(time,XI)
            grid on;

            linkaxes([subplot(2, 1, 1), subplot(2, 1, 2)], 'x','y');
            set(gcf, 'Position', get(0, 'Screensize'));

            % figure; plot(n_samples,XI)
            if ~exist('x_start', 'var')
                fprintf('Haz clic en la gráfica para seleccionar los puntos de la parte de marcha.\n');
                fprintf('Presiona Enter cuando hayas terminado.\n');
                [x_index, ~] = ginput(2);
    
                
               if ~isempty(x_index)
                    x_start = x_index(1);
                    x_end = x_index(2);
                    XI = XI(x_start*FS-time(1)*FS:x_end*FS-time(1)*FS);
               end
            else
                XI = XI(x_start*FS-time(1)*FS:x_end*FS-time(1)*FS);
                close all
            end

            % close all

            %SELECT WALKING TRIALS
        
            [N,~] = size(XI);
            [M,~] = size(XJ);
    
            % We want to do the DTW in both dimensions, so test each point of the
            % template against each point of the unknown signal
            % tic
            fun = @(a,b) sqrt((a-b).^2);
            Dx = flipud(bsxfun(fun,XI(1:N,1),XJ(1:M,1)')');
            % Dy = flipud(bsxfun(fun,YI(1:N,1),YJ(1:M,1)'))';
            % Dz = flipud(bsxfun(fun,ZI(1:N,1),ZJ(1:M,1)')');
            % toc
        
            % Test Adding the Distance matrices together 
            D =  Dx; %+ Dy + Dz ; %
            clear Dx %Dy Dz
            % Build the accumulated cost matrix C from the Distance Matrix D to search
            % for the optimal warping path p according to Barth et al. Sensors 2015
            Ctemp = zeros(M,N);
        
            % The bottom row of C is equal to the bottom row of D
            Ctemp(M,:) = D(M,:);
        
            % The first column of C is calculated by summing the values from the bottom
            % to the top of the first column of D
            for i = 1:M-1
                Ctemp(M-i,1) = sum(D(M-i:M,1),1);
            end
        
            % The remaining elements of C are calculated by adding the minimum value of
            % the three neighboring elements to the matrix element of D
        
            for j = 2:N %length of the movement sequence  
                for i = 1:M-1 %length of the template-1
                    Ctemp(M-i,j) = min([Ctemp((M+1)-i,(j-1)),Ctemp((M-i),(j-1)),Ctemp(((M+1)-i),j)]) + D(M-i,j);
                end
            end
            % HeatMap(Ctemp)
            % figure; plot(Ctemp(1,:))
            
            % Save Ctemp as the cost-matrix C for each window and clear Ctemp
            DTWResults_multTemp.C = Ctemp;
            clear Ctemp D
        
            % Save the Distance Function which is the accumulated distances along
            % the top row of C (at the end of any possible stride cycle) for
            % each window
            DTWResults_multTemp.DistFnct = DTWResults_multTemp.C(1,:);
            % figure; plot(DTWResults_multTemp.DistFnct)
    
            % Define the threshold value Tau for stride segmentation where each
            % selected minimum representing the end of a stride has to be smaller than
            % the threshold. 
                % ThreshVal = ThresVals; 
            
            % Find the local minima which are smaller than Tau - 
            % Because Matlab doesn't have a minima function, we have to flip Tau and
            % the signal based on the maxima and then find any peaks above this new
            % threshold. Additionally, because we are trying to use a generic minima
            % rather than a patient specific value, we can set a minimum distance 
            % between peaks to aid in the differentiation of strides. We will chose 50 
            % as the data was collected at 100Hz (ie 100 frames per second) so a
            % participant would have to be taking >2 steps in less than half a second
            % to not satisfy this threshold which is highly unlikely given the
            % end-stage OA population.
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flip_data = 1.01*max(DTWResults_multTemp.DistFnct)-DTWResults_multTemp.DistFnct;
            figure; plot(flip_data); title('Select the minimum threshold');
            Tau = input('Introduzca el umbral para los mínimos: ', 's');
            Tau = str2double(Tau);
            close all


            % Encuentra los mínimos locales
            [minPeaks, minLocs] = findpeaks(-flip_data, 'MinPeakHeight', -Tau, 'MinPeakDistance',5000);%, 'SortStr', 'descend');%, 'MinPeakDistance',1000);
            minPeaks = -minPeaks;

            figure;
            plot(flip_data)
            hold on;
            % plot(maxIndexBeforeMin, localMaxBeforeMin, 'x', 'MarkerSize', 10, 'DisplayName', 'Picos');
            plot(minLocs, minPeaks, 'x', 'MarkerSize', 10, 'DisplayName', 'Picos');
            hold on;

            % figure; plot(flip_data); title('Select the maximum threshold');
            Tau = input('Introduzca el umbral para los máximos: ', 's');
            Tau = str2double(Tau);

            [maxPeaks, maxLocs] = findpeaks(flip_data, 'MinPeakHeight', Tau, 'MinPeakDistance',600);
            % figure;
            % plot(flip_data)
            % hold on;
            % % plot(maxIndexBeforeMin, localMaxBeforeMin, 'x', 'MarkerSize', 10, 'DisplayName', 'Picos');
            plot(maxLocs, maxPeaks, 'x', 'MarkerSize', 10, 'DisplayName', 'Picos');

            % Inicializar las variables para almacenar los resultados
            localMaxBeforeMin = [];
            maxIndexBeforeMin = [];
            validMax = [];

            % Busca los máximos antes de los mínimos
            % start = 1
            % for s = 1:length(maxLocs):
            %     [MinimaTemp, locs] = findpeaks(flip_data(, Tau);

            for i = 1:length(minLocs)
                % Encuentra el mínimo local que viene después del máximo actual
                
                if i > 1
                    prevMaxIdx = find(maxLocs < minLocs(i) & maxLocs > minLocs(i-1));
                else
                    prevMaxIdx = find(maxLocs < minLocs(i));
                end
                % [valor_maximo, prevMaxIdx] = max(maxLocs < minLocs(i));

                if ~isempty(prevMaxIdx)
                    % Guarda el máximo y su índice si cumple la condición
                    indexDifference = minLocs(i) - maxLocs(prevMaxIdx);
                    validMax = prevMaxIdx(indexDifference <= 6000);
                    max_validMax = max(maxPeaks(validMax));
                    %prevMaxIdx_one = find(maxPeaks == max_validMax);
                    % Guarda el máximo y su índice si cumple la condición
                    % localMaxBeforeMin = [localMaxBeforeMin; maxPeaks(prevMaxIdx_one)];
                    maxIdx = find(maxPeaks == max_validMax);
                    % valor_maximo = max(maxLocs(prevMaxIdx));
                    localMaxBeforeMin = [localMaxBeforeMin; maxPeaks(maxIdx)];
                    % maxIndexBeforeMin = [maxIndexBeforeMin; maxLocs(prevMaxIdx_one)];
                    maxIndexBeforeMin = [maxIndexBeforeMin; maxLocs(maxIdx)];
                    % maxLocs = maxLocs(size(prevMaxIdx)+1:end);
                end
                prevMaxIdx = [];
                % prevMaxIdx_one = []; 

            end

            figure;
            plot(flip_data)
            hold on;
            plot(maxIndexBeforeMin, localMaxBeforeMin, 'x', 'MarkerSize', 10, 'DisplayName', 'Picos');
            plot(minLocs, minPeaks, 'x', 'MarkerSize', 10, 'DisplayName', 'Picos');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % We sum to time the part we have removed from the start of the
            % signal
            if exist("x_start")
                gait_events = table(sort(time(maxIndexBeforeMin)+(x_start-time(1))),'VariableNames',{strcat('heel_strike_',muscle_template(end))});
            else
                gait_events = table(sort(time(maxIndexBeforeMin)),'VariableNames',{strcat('heel_strike_',muscle_template(end))});
            end
            first_time = false;
            if isempty(total_events)
                first_time = true;
            end

            
            len1 = height(gait_events);
            len2 = height(total_events);
            maxLen = max(len1, len2);

            % Rellenar con NaN para igualar las longitudes
            if ~first_time
                if len1 < maxLen
                    gait_events{len1+1:maxLen,1} = NaN;
                    
                elseif len2 < maxLen
                    total_events{len2+1:maxLen,1} = NaN;
                end
            end

            total_events = horzcat(gait_events, total_events);

            x_index = [];
            clearvars x_start x_end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Here I suggest you pick the maximum peak from each of the
            % identified blocks (you'll see you have repetitions of values
            % that are very close together). Either you can find the max or
            % you can smooth the data and the template before the DTW and
            % then the peak should just be 1 value.
            % Then you can use the points identified as your segmentation
            % points for start and end.
            % Basically ignore everything after this. 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

%% 
close all
writetable(total_events,strcat(path_templates,name_file,'_estimatedGaitEvents.csv'))
    %% 
    %% 
    %% 
    
%         % Added Nov14, in order to hone in on the actual steps being take, we
%         % want to find a non-generic threshold once we've done the first pass
%         % with the generic value. To do this, we wil enforce the cutoff
%         % specifying that steps must be longer than 50frames. 
% 
%         % The beauty here is that the minima in the distance function 
%         % calculated should be the smallest when a step is completed (since the
%         % first row of the distance matrix equates to the cost of warping the
%         % last point of the template to the trial signal). So while we started
%         % with a generic template, we can use the highest peaks (which equate
%         % to the lowest mins in the OG data) to find those cutpoints. 
% 
%         % Since the MinimaTemp (which is actually peaks in the flipped data)
%         % have been sorted in decending order, we can find the height of the
%         % last peak which occurs before the cutoff of an unrealistic number of
%         % steps (calculated by dividing the length of the trial by 50)
%             maxNumSteps = round(length(...
%                 DTWResults_multTemp.DistFnct)/50)+1;
%             if length(MinimaTemp)> maxNumSteps
%                 newTau = MinimaTemp(round(length(...
%                     DTWResults_multTemp.DistFnct)/50)+1);
%             else
%                 newTau = Tau;        
%             end
%         % Then we find new peaks and peak indices without sorting using this
%         % new Tau value. 
%             [MinimaTempNew, MinIdxTemp] = findpeaks((1.01*max(...
%                 DTWResults_multTemp.DistFnct)-...
%                 DTWResults_multTemp.DistFnct),...
%                 'MinPeakHeight',newTau);
% 
%             clear MinimaTemp  
%             MinimaTemp = MinimaTempNew;
%             if ~isempty(MinIdxTemp)
% 
%             % Instead of using MinPeakDistance I think we need to go back through 
%             % the peaks once they're identified and then say, if they are within 
%             % 50 frames of eachother, pick the max peak. 
%                 threshold = 50;
%                 transitions = diff([0, [(MinIdxTemp(2:end) - MinIdxTemp(1:end-1)), 0] < threshold]);
% 
%                 % find the runs where values all fell below the threshold
%                 runstarts = find(transitions == 1);
%                 runends = find(transitions == -1);
% 
%             % Right now this works if the last run goes to the last step, but there
%             % was an issue where the only transition==1 happened at the start of
%             % the trial in which case we don't want to use the loop... if it
%             % becomes a problem it'll need to be revisited
%                 if length(runends)<length(runstarts)
%                     runends(1,length(runends)+1) = length(transitions);
%                 end
% 
%             % Use the runs to identify blocks of peaks occuring within 50 frames of
%             % eachother
%                 blocks = [];
%                 blocks = cell2mat(arrayfun(@(s,e) vertcat(blocks, s:1:e), runstarts, runends, 'UniformOutput', false));
% 
%             % Find the peak of all peaks for each block
%                 PeakInBlocks = (arrayfun(@(s,e) (find(MinimaTemp(s:e) == max(MinimaTemp(s:e)))+s-1), ...
%                 runstarts, runends, 'UniformOutput', false));
% 
%             % Assign the index of each true peak to MinIdx
%                 MinIdx = sort([MinIdxTemp(~ismember(1:1:length(MinIdxTemp),blocks)), MinIdxTemp(1,[PeakInBlocks{:}])]);
% 
%             %  Calculate the warping path by tracing back from the minima 
%             %  identified along the lowest cost path (the warping path) until we 
%             %  reach the starting point of each stride
% 
%             %  We note that the start of the warping path is in the bottom row of C 
%             %  and represents the end of the stride cycle while the end of the 
%             %  warping path is in the top row of C and represents the start of the 
%             %  stride cycle. Further, the warping path p has to be a monotonic
%             %   function where only neighboring elements are added, and at every 
%             %   new element, at least one index must decrease.
% 
%                 % tic
%                 for i = 1:length(MinIdx)
%                     temppath = [];
%                     j = M;
%                     k = MinIdx(1,i);
%                     cnt = 2;
%                     % in order to constrain the warping path, we want to make sure we
%                     % are not deviating too far from the diagonal, let's start by
%                     % allowing a repetition in the same direction no more than 3 times
%                     cntj = 0;
%                     cntk = 0;
%                     temppath(1,1:2) = [j,k];
%                     while j >= 1
%                         if j >1
%                             if k == 1
%                                 temppath(cnt,[1,2]) = [1,1];
%                                 j = 0;
%                             else
%                                 [~,testidx] = min([...
%                                     DTWResults_multTemp.C(201-j,(k-1)),...
%                                     DTWResults_multTemp.C((202-j),(k-1)),...
%                                     DTWResults_multTemp.C((202-j),k,1)]);
%                                 if testidx == 1
%                                     k = k-1;
%                                     cntk = cntk+1;
%                                     countj = 0;
%                                     if cntk == 10
%                                         j = j-1;
%                                         cntk = 0;
%                                     end
%                                 elseif testidx == 2
%                                     k = k-1; 
%                                     j = j-1;
%                                     cntk = 0; cntj = 0;
%                                 else
%                                     j = j-1;
%                                     cntj = cntj+1;
%                                     countk = 0;
%                                     if cntj == 3
%                                         k = k-1;
%                                         cntj = 0;
%                                     end
%                                 end
%                                 temppath(cnt,[1,2]) = [j,k];
%                                 cnt = cnt+1;
%                             end
%                         elseif j == 1 && k == 1
%                             % the loop should end here since we've reached the start
%                             % of the window and can't move any further in any
%                             % direction.
%                             j = 0; k = 0;
%                         else
%                             [~,testidx] = min([...
%                                 DTWResults_multTemp.C(201-j,(k-1)),...
%                                 DTWResults_multTemp.C(201-j,k,1)]);
%                             if testidx == 1
%                                 if k == 1
%                                     temppath(cnt,[1,2]) = [j,k];
%                                     j = 0;
%                                 else
%                                     k = k-1;
%                                     temppath(cnt,[1,2]) = [j,k];
%                                 end
%                             else
%                                 j = j-1;
%                             end
%                             cnt = cnt+1;
%                         end
%                     end
%                     % size(temppath)
%                     % Save the start and end of the warping path for each identified 
%                     % stride in temppath, making sure to add in the appropriate value 
%                     % for the start of the window such that these data will correspond 
%                     % to their appropriate position in XIc
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(1:2,i) =...
%                         [temppath(end,2)+ind(wnd,1);temppath(1,2)+ind(wnd,1)];
%             %         plot(temppath(:,2),temppath(:,1),'g')
%             %         hold on
% 
%                     % On Nov 19 tried saving the warping path and the cost of
%                     % warping summed along the path to compare between steps which
%                     % appear to have reasonable lengths but start at very different
%                     % times in the stride cycle but this didn't seem to help detect
%                     % the correct stride. 
%                     % Idea 2 is to save the minimum cost at the final frame of the
%                     % stride cycle and to compare those.
%                     % tempcost = 0;
%                     % tic
%                     % for npath = 1:length(temppath)
%                     %     tempcost = tempcost + ...
%                     %         DTWResults_multTemp.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['Window' num2str(wnd)]).C(202-temppath(npath,1),temppath(npath,2));
%                     % end
%                     % toc
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCosts(1,i) = ...
%                         DTWResults_multTemp.C(202-temppath(1,1),temppath(1,2));
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCostsStartofStride(1,i) = ...
%                         DTWResults_multTemp.C(temppath(end,1)+100,temppath(end,2));
% 
%                     clear temppath cnt testidx 
%                 end
%                 % toc
%                 clear blocks PeakInBlocks runstarts runends Tau MinimaTemp ...
%                          MinIdxTemp MinIdx transitions
% 
%                 % Check the warping paths that were detected to ensure they are all
%                 % capturing only 1 stride, and that the stride captured is
%                 % complete. 
% 
%                 % To do this we will first ensure that if the start of any path
%                 % isn't equal to the end of the previous path, but they are within
%                 % 50 frames of each other, that the values chosen are appropriate 
%                 % or if they should be adjusted
%                 WPStartEndDiffs = (WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(1,2:end)-...
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,1:end-1));
%                 ChekWPStartInds = ismember(WPStartEndDiffs,WPStartEndDiffs(WPStartEndDiffs~=0 & WPStartEndDiffs<50));
%                 WPchecks = find(ChekWPStartInds>0);
%                 for ichecks = 1:length(WPchecks)
%                     if WPStartEndDiffs(WPchecks(ichecks)) >0
%                         [~,pkloc] = findpeaks(XIc((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,WPchecks(ichecks))-5):...
%                             (WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(1,WPchecks(ichecks)+1)+5)));
%                         if isempty(pkloc)
%                             pkloc = 1;
%                         end
%                         maxInd = pkloc(end) + WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,WPchecks(ichecks))-6;
%                     else
%                         [~,pkloc] = findpeaks(XIc((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(1,WPchecks(ichecks)+1)-5):...
%                             (WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,WPchecks(ichecks))+5)));
%                         if isempty(pkloc)
%                             pkloc = 1;
%                         end
%                         maxInd = pkloc(end) + WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(1,WPchecks(ichecks)+1)-6;
%                     end
%                     if maxInd < ind(wnd,1)
%                         maxInd = ind(wnd,1);
%                     end
%                     % update both the end of the previous stride and the start of
%                     % the current stride with the max value between them
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,WPchecks(ichecks)) = maxInd;
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(1,WPchecks(ichecks)+1) = maxInd;
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCostsStartofStride(1,WPchecks(ichecks)+1) = ...
%                         DTWResults_multTemp.C(200,(maxInd - ind(wnd,1)+1));
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCostsEndofStride(1,WPchecks(ichecks)) = ...
%                         DTWResults_multTemp.C(1,(maxInd - ind(wnd,1)+1));
%                     clear pkloc maxInd 
%                 end
% 
%                 % Now ensure that the length of each warping path is longer than
%                 % the mean+2stdev of all calculated warping paths.
%                 WPlengths = WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,:)-...
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(1,:);
%                 meanWPs = mean(WPlengths);
%                 stdWPs = std(WPlengths);
% 
%                 if stdWPs > meanWPs/4
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(:,...
%                         (WPlengths<(meanWPs-stdWPs))) = [];
%                      WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCosts(:,...
%                         (WPlengths<(meanWPs-stdWPs))) = [];
%                      WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCostsStartofStride(:,...
%                         (WPlengths<(meanWPs-stdWPs))) = [];
%                 else
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(:,...
%                         (WPlengths<(meanWPs-2*stdWPs))) = [];
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCosts(:,...
%                         (WPlengths<(meanWPs-2*stdWPs))) = [];
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCostsStartofStride(:,...
%                         (WPlengths<(meanWPs-2*stdWPs))) = [];
%                 end
% 
%                 % Lastly, we need to combine the strides found in each window into
%                 % a single variable for the entire walking bout while avoiding
%                 % duplicates or overlap!
%                 if wnd == 1 || wndcond == 1
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides = ...
%                         WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(:,2:end);
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.WarpingCosts = ...
%                         WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCosts(:,2:end);
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.WarpingCostsStartofStride = ...
%                         WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCostsStartofStride(:,2:end);
%                 else
%                     %%%%%%%%%%%%%%%%%%%%% NOV 19 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     % sub=3 length mismatch wp b and t in window 6
% 
%                     % find all non-duplicated points in the stride endings
%                     [WPb,ia,ib] = ...
%                         union(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:),...
%                         WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,2:end));
%                     % if there were no strides found in this window (or no
%                     % new strides) ib will be empty and we should dump out
%                     % of this loop now. 
%                     if isempty(ib)
%                         % do nothing
%                     else
%                         % because we only want to start at the second stride in the new
%                         % window we must add one to ib
%                         ib = ib+1;
%                         % find all non-duplicated points in the stride beginnings to
%                         % ensure that there isn't a mismatch between rows
%                         % [WPttemp,iat,ibt] = ...
%                         %     union(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).AccX.Strides(1,:),...
%                         %     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).AccX.(['Window' num2str(wnd)]).Strides(1,:));
% 
%                         % Nov 21 - What it seems like we're missing here is that there
%                         % are some duplicated end points and right now the above code 
%                         % "union" favours the 1st input (which is the strides previously 
%                         % compiled from other windows) rather than selecting the
%                         % appropriate stride between the 2.
%                         % Instead, we should also find for all instances of the
%                         % bottom being the same, in which instances are the top values 
%                         % (the start of the stride) different, and in those instances, 
%                         % compare the strides and determine which is ultimately
%                         % correct. We can do this using the cost of warping the start
%                         % of the template to the detected start of the stride cycle.
%                         % Whichever start point has the lowest cost should be the value
%                         % moved forward. 
%                         bIndDup = ismember(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:),...
%                             WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,2:end));
%                         tIndDup = ismember(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(1,:),...
%                             WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(1,2:end));
%                         % Specify the elements of ia which might need to be replaced by 
%                         % new elements in ib
%                         ConflictStrides_inA = find(bIndDup==1 & bIndDup~=tIndDup);
%                         % Find the corresponding elements in ib to compare with
%                         ConflictStrides_inB = find(ismember(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,2:end),...
%                             WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,(bIndDup==1 & bIndDup~=tIndDup))))+1; %+1 b/c we are ignoring the first step
% 
%                         for fix = 1:length(ConflictStrides_inA)
%                             if WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.WarpingCostsStartofStride(1,ConflictStrides_inA(fix))>...
%                                    WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCostsStartofStride(1,ConflictStrides_inB(fix))
%                                 ia(ia==ConflictStrides_inA(fix)) = [];
%                                 ib = sort([ib;ConflictStrides_inB(fix)]);
%                                 % NOTE we don't need an else statement because if the
%                                 % stride in ia was the correct choice, then we just
%                                 % leave things the way they were!
%                             end
%                         end
% 
%                         % find any points which are within 50 frames of eachother
%                         % (following our previous thresholding)
%                         rpts = diff([0, [(WPb(2:end) - WPb(1:end-1)), 0] < threshold]);
%                         rptS = find(rpts == 1);
%                         rptE = find(rpts == -1);
%                         cntWP = sum(abs(rpts(1:end-1)))/2;
% 
%                         for cntr = 1:cntWP
%                             if abs(rptE(cntr)-rptS(cntr))>1
%                                 rptE(cntr) = rptE(cntr)-1;
%                             end
%                             % First check if the distance between the detected start
%                             % point is non-trivial (ie that we are starting from a
%                             % differnet peak in the stride cycle rather than 2 points
%                             % along the same peak
%                             if abs(WPb(rptS(cntr)-cntr+1) - WPb(rptE(cntr)-cntr+1)) > 10
%                                 if any(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:)==WPb(rptS(cntr)-cntr+1))
%                                     Cost1 = WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.WarpingCosts(...
%                                         ia(ia==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:)==...
%                                             WPb(rptS(cntr)-cntr+1))~=0))));
%                                     Cost2 = WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCosts(...
%                                         ib(ib==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,:)==...
%                                             WPb(rptE(cntr)-cntr+1))~=0))));
%                                     if ~isempty(Cost1)
%                                         if Cost1 < Cost2
%                                             ib(ib==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,:)==...
%                                                 WPb(rptE(cntr)-cntr+1))~=0))) = [];
%                                             WPb(rptE(cntr)-cntr+1) = [];
%                                         else
%                                             ia(ia==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:)==...
%                                                 WPb(rptS(cntr)-cntr+1))~=0))) = [];
%                                             WPb(rptS(cntr)-cntr+1) = [];
%                                         end
%                                     else
%                                         ib(ib==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,:)==...
%                                             WPb(rptE(cntr)-cntr+1))~=0))) = [];
%                                         WPb(rptE(cntr)-cntr+1) = [];
%                                     end
%                                 else
%                                     Cost1 = WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.WarpingCosts(...
%                                         ia(ia==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:)==...
%                                             WPb(rptE(cntr)-cntr+1))~=0))));
%                                     Cost2 = WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCosts(...
%                                         ib(ib==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,:)==...
%                                             WPb(rptS(cntr)-cntr+1))~=0))));
%                                     if ~isempty(Cost1)
%                                         if Cost1 < Cost2 
%                                             ib(ib==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,:)==...
%                                                 WPb(rptS(cntr)-cntr+1))~=0))) = [];
%                                             WPb(rptS(cntr)-cntr+1) = [];
%                                         else
%                                             ia(ia==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:)==...
%                                                 WPb(rptE(cntr)-cntr+1))~=0))) = [];
%                                             WPb(rptE(cntr)-cntr+1) = [];
%                                         end
%                                     else
%                                         ib(ib==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,:)==...
%                                             WPb(rptS(cntr)-cntr+1))~=0))) = [];
%                                         WPb(rptS(cntr)-cntr+1) = [];
%                                     end
%                                 end
%                             % Keep the point which represents a smaller value in XIc and
%                             % remove the other from WPb and from the appropriate ia/ib
%                             elseif abs(XIc(WPb(rptS(cntr)-cntr+1),1)) > abs(XIc(WPb(rptE(cntr)-cntr+1),1))
%                                 % here the value which started the close pair is larger
%                                 % so we get rid of that one
%                                 if any(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:)==WPb(rptS(cntr)-cntr+1))
%                                     ia(ia==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:)==...
%                                         WPb(rptS(cntr)-cntr+1))~=0))) = [];
%                                 else
%                                     ib(ib==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,:)==...
%                                         WPb(rptS(cntr)-cntr+1))~=0))) = [];
%                                 end
%                                 WPb(rptS(cntr)-cntr+1) = [];
%                             else
%                                 % here the value which ended the close pair is larger
%                                 % so we ger rid of that one.
%                                 if any(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:)==WPb(rptE(cntr)-cntr+1))
%                                     ia(ia==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:)==...
%                                         WPb(rptE(cntr)-cntr+1))~=0))) = [];
%                                 else
%                                     ib(ib==(find((WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(2,:)==...
%                                         WPb(rptE(cntr)-cntr+1))~=0))) = [];
%                                 end
%                                 WPb(rptE(cntr)-cntr+1) = [];
%                             end
%                         end
% 
%                         % Now use ia/ib to define the top row and save the new warping
%                         % path for these combined windows
%                         WPt = sort([WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(1,ia'),...
%                             WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides(1,ib')]);
%                         WPcost = [WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.WarpingCosts(1,ia'),...
%                             WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCosts(1,ib')];
%                         WPcostS = [WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.WarpingCostsStartofStride(1,ia'),...
%                             WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).WarpingCostsStartofStride(1,ib')];
%                         if length(WPb) < length(WPt)
%                             WPt = WPt(1:length(WPb));
%                             WPcost = WPcost(1:length(WPb));
%                             WPcostS = WPcostS(1:length(WPb));
%                         end
%                         % Finally, replace the previous Strides array with this new longer one
%                          WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides = [WPt;WPb];
%                          WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.WarpingCosts = WPcost;
%                          WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.WarpingCostsStartofStride = WPcostS;
%                          clear WPt WPb WPcostS WPcost WPttemp rptE rptS 
%                     end
%                 end
%                 % Save the number of strides taken in this window
%                 strideNums(sub,1) = size(...
%                     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.(['Window' num2str(wnd)]).Strides,2);
%                 clear XI
%                 wndcond = 0;
%                 WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX = ...
%                 rmfield(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX,(['Window' num2str(wnd)]));
%             else
%                 if exist('WarpingPaths','var')
%                     if isfield(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]),(['FreeLiving_day' num2str(FLDay)]))
%                         if isempty(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides)
%                             wndcond = 1;
%                         else
%                             %do nothing
%                         end
%                     else
%                         wndcond = 1;
%                     end
%                 else
%                     wndcond = 1;
%                 end
%             end
%         end
%         clear DTWResults_multTemp
%% 
% To keep only the "right" strides, we could run a basic dtw here and
% keep only those steps with a low cost maybe? 
% DTWResults_Strides.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Distance = ...
%     zeros(1,length(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides));
DTWResults_Strides.Distance = zeros(1,length(ShankAccTemplateStrides.Strides)); %%PAIRED VARIABLE
% for nStrides = 1:length(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides)
%     DTWResults_Strides.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Distance(1,nStrides) = ...
%         dtw(normalize(XIc(WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(1,nStrides):...
%         WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,nStrides),1),'range',[-1 1]),XJ');
% end
for nStrides = 1:length(ShankAccTemplateStrides.Strides)
    DWTResults_multTemp.DistFnct(1,nStrides) = dwt(normalize(XIc(ShankAccTemplateStrides.Strides.x(1,nStrides):ShankAccTemplateStrides.Strides(2,nStrides),1),'range',[-1 1]),XJ');
end
%Discard any steps whose dtw distance exceeds the mean+stdev
% meanWPsDTW = mean(DTWResults_Strides.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Distance,"omitnan");
% stdWPsDTW = std(DTWResults_Strides.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Distance,"omitnan");

meanWPsDTW = mean(DWTResults_multTemp.DistFnct,"omitnan");
stdWPsDTW = std(DWTResults_multTemp.DistFnct,"omitnan");

% WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(:,...
%     (DTWResults_Strides.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Distance>(meanWPsDTW+stdWPsDTW))) = [];

ShankAccTemplateStrides.Strides(:,...
    (DWTResults_multTemp.DistFnct>(meanWPsDTW+stdWPsDTW))) = [];

% figure; hold on; 
% for i = 1:size(WarpingPaths.(['ROB0' Participants{sub}]).Session1.(['FreeLiving_day' num2str(FLDay)]).AccX.Strides,2) 
%     plot(XIc(WarpingPaths.(['ROB0' Participants{sub}]).Session1.(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(1,i):...
%         WarpingPaths.(['ROB0' Participants{sub}]).Session1.(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,i)))
% end

figure; hold on; 
for i = 1:size(ShankAccTemplateStrides.Strides,2) 
    plot(XIc(ShankAccTemplateStrides.Strides(1,i):...
        ShankAccTemplateStrides.Strides(2,i)))
end

% 
% %         % As a final step, remove any strides which fall below mean-2*stdev of
% %         % the collective strides found in case any windows resulted in strides
% %         % which didn't match the quality of the others
% WPlengthsC = WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).AccX.Strides(2,:)-...
%     WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).AccX.Strides(1,:);
% meanWPsC = mean(WPlengthsC);
% stdWPsC = std(WPlengthsC);
% 
% WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).AccX.Strides(:,...
%     (WPlengthsC<(meanWPsC-2*stdWPsC))) = [];
% WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).AccX.WarpingCosts(:,...
%     (WPlengthsC<(meanWPsC-2*stdWPsC))) = [];

WPlengthsC = ShankAccTemplateStrides.Strides(2,:)-...
    ShankAccTemplateStrides.Strides(1,:);
meanWPsC = mean(WPlengthsC);
stdWPsC = std(WPlengthsC);

ShankAccTemplateStrides.Strides(:,...
    (WPlengthsC<(meanWPsC-2*stdWPsC))) = [];
WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).AccX.WarpingCosts(:,...
    (WPlengthsC<(meanWPsC-2*stdWPsC))) = [];
% 
%         % Save the numer of strides taken across the entire walking bout
% 
%         strideNums(sub,1) = size(...
%             WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides,2);
% 
%         clear XIc DTWResults_Strides MinimaTempNew ind
% 
%         HeelStrikeIndsISTA = [WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(1,:)',...
%             WarpingPaths.(['ROB0' Participants{sub}]).(['Session' Session{sub}]).(['FreeLiving_day' num2str(FLDay)]).AccX.Strides(2,:)'];
% 
%         DPclinical = ['C:\Users\annem\OneDrive - Dalhousie University\ORS\OrthoRoboticsIMUData\ROB0' ...
%                 (Participants{sub}) '\Session' (Session{sub}) '\'];
% 
%         writematrix(HeelStrikeIndsISTA, [DPclinical 'HeelStrikeIndex_FLday' num2str(FLDay) 'ISTA2025.csv'])
% 
%         clear WarpingPaths HeelStridIndsISTA
%         toc
%     %     end
%         mail = ("an309604@dal.ca");
%         t_body = (join(["Code Finished Running Sub" Participants{sub} "Day" num2str(FLDay)]));
%         subj = ("Matlab Update");
%         sendolmail(mail,subj,t_body)
%     end    
%     mail = ("an309604@dal.ca");
%     t_body = (join(["Code Finished Running Sub" Participants{sub}]));
%     subj = ("Matlab Update");
%     sendolmail(mail,subj,t_body)
% end
% toc
% mail = ("an309604@dal.ca");
% t_body = ("Code Finished Running");
% subj = ("Matlab Update");
% sendolmail(mail,subj,t_body)
% 
% %% Now that we have these indices, we can save files for each participant
% %  and align them in Python prior to doing the variability comparison. 
% tic
% DPclinical = 'C:\Users\annem\OneDrive - Dalhousie University\ORS\OrthoRoboticsIMUData\ROB0';
% FullDataPath = 'C:\Users\annem\OneDrive - Dalhousie University\RoboticsStudyData_UpToOarsiAnalysis\IMUData\';
% 
% for sub = 1:length(Participants)
%     for FLDay = 2:7
%         FullSubDataPath = [FullDataPath 'ROB0' (Participants{sub}) '\Session' (Session{sub}) '\'];
%         FullDayGaitDataFile = [FullSubDataPath 'FreeLiving_day' num2str(FLDay) '.csv']; %' (Days{sub}) '
%         FDShankAcc = readmatrix(FullDayGaitDataFile);
%         HSIndPath = [DPclinical (Participants{sub}) '\Session' (Session{sub}) '\HeelStrikeIndex_FLday' num2str(FLDay) 'ISTA2025.csv'];
%         HSInd = readmatrix(HSIndPath);
% 
%         XRawISTA = []; YRawISTA = []; ZRawISTA = [];
% 
%         for stps = 1:length(HSInd)
%             XRawISTA = [XRawISTA;FDShankAcc(HSInd(stps,1):HSInd(stps,2),1)];
%             YRawISTA = [YRawISTA;FDShankAcc(HSInd(stps,1):HSInd(stps,2),2)];
%             ZRawISTA = [ZRawISTA;FDShankAcc(HSInd(stps,1):HSInd(stps,2),3)];
%         end
% 
%         GaitISTA = zeros(length(XRawISTA),42);
%         GaitISTA(:,7:9) = [XRawISTA, YRawISTA, ZRawISTA];
%         writematrix(GaitISTA, [DPclinical (Participants{sub}) '\Session' (Session{sub}) '\Gait_FLday' num2str(FLDay) 'ISTA2024.csv'])
% 
% 
%         clear HeelStrikeIndsISTA
%         clear XRawISTAYRawISTA ZRawISTAGaitISTA2025
%     end
%     mail = ("an309604@dal.ca");
%     t_body = (join(["Code Finished Running Sub" Participants{sub}]));
%     subj = ("Matlab Update");
%     sendolmail(mail,subj,t_body)
% end
% %% last save the stride indices for the new gait sessions for python
% tic
% for sub = 1:length(Participants)
%     for FLDay = 2:7
%         HSIndPath = [DPclinical (Participants{sub}) '\Session' (Session{sub}) '\HeelStrikeIndex_FLday' num2str(FLDay) 'ISTA2025.csv'];
%         HSInd = readmatrix(HSIndPath);
% 
%         for stps = 1:length(HSInd)
%             if stps == 1
%                 NewHSInd = [1,HSInd(stps,2)-HSInd(stps,1)];
%             else
%                 NewHSInd = [NewHSInd;NewHSInd(stps-1,2),HSInd(stps,2)-HSInd(stps,1)+NewHSInd(stps-1,2)];
%             end
% 
%         end
% 
%         writematrix(NewHSInd, [DPclinical (Participants{sub}) '\Session' (Session{sub}) '\HSPythonIndex_FLday' num2str(FLDay) 'ISTA2024.csv'])
% 
%         clear HSInd NewHSInd HSIndPath
%     end
%     mail = ("an309604@dal.ca");
%     t_body = (join(["Code Finished Running Sub" Participants{sub}]));
%     subj = ("Matlab Update");
%     sendolmail(mail,subj,t_body)
% end
% 
% toc
% %% functions
% 
% function sendolmail(to,subject,body,attachments) 
% %function sendolmail(to,subject,body,attachments)
% %Sends email using MS Outlook. The format of the function is  
% %Similar to the SENDMAIL command. 
% 
% % Create object and set parameters. 
% h = actxserver('outlook.Application'); 
% mail = h.CreateItem('olMail'); 
% mail.Subject = subject; 
% 
% % multiple recipients 
% if length(to) > 1 
%       to = strjoin(to,';'); 
% end
% mail.To = to; 
% 
% mail.BodyFormat = 'olFormatHTML'; 
% mail.HTMLBody = body; 
% 
% % Add attachments, if specified. 
% if nargin == 4 
%     for i = 1:length(attachments) 
%         mail.attachments.Add(attachments{i}); 
%     end 
% end 
% 
% % Send message and release object. 
% mail.Send; 
% h.release; 
% 
% end
