%% This code has been written in order to take data collected from free-
%  living axivity sensor accelerometer data collected as part of the
%  greater Ortho robotics study and identify individual HS (heel strike)
%  events from identified walking bouts such that these indices can be used
%  to create templates from which more stride cycles can be identified from
%  the complete free-living data
%  
%
% This Program was written by:   Kaitlin Genge
% Created on:                    2023
% Last updated by:               Annemarie Laudanski
% Last updated on:               May 27, 2024
% Last updated by:               Marina Algaba Vidoy
% Last updated on:               May 31, 2024
clc; clear; %close all

%% Setup

Datapath = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/MEDIAN_FILTERED/';
folder_out = '/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/GAIT_WITHOUT_EXO/PROCESSED/TEMPLATES/';
% DPclinical = ['C:\Users\annem\Documents\ORS\OrthoRoboticsIMUData\'];
% DPclinical = 'C:\Users\annem\OneDrive - Dalhousie University\Matlab\CSIC\Marina - DTW\';


% % Clinical Data
% Participants = {'02','05','07','08','10','11','12','13','14','15','16','18','19',... %1-13
%     '20','21','22','23','24','25','26','32','34'}; %14-22 participants
% % Data Collection Session
% Session = {'1','1','1','1','1','1','1','1','1','1','1','1','1','1','1',...
%     '1','1','1','1','1','1','1'};
% % Day of collection from which steady-state walking bout was identified
% Days = {'3','4','4','4','4','4','4','4','3','3','4','3','3','3','3','3',...
%     '3','3','4','4','4','4'};


% % Gait walking bouts start and stop
% WindowStartG = [5661000, 3662340, 4332560, 3479630, 5441630, 5811340, 5192810, 4191710, 6330020, 4425130,...
%     3958201, 3635390, 3800280, 3272180, 5933400, 4345450,5138060, 4312840, 5162710, 5597080, 6054760, 4302650]; 
% WindowEndG = [5664360, 3666930, 4336060, 3483150, 5444720, 5814360, 5197060, 4194160, 6333760, 4428790,...
%     3961690, 3638310, 3803610, 3274440, 5935910, 4347800, 5139900, 4315120, 5165320, 5599540, 6057610, 4305920];
% 
% % Standing bout start and stop
% WindowStartS = [6332380, 4050660, 3899720, 5649210, 4796070, 3769580, 3963690, 5479180, 7132260, 4805040,...
%     3963140, 3638500, 3798110, 3270640, 5931140, 4348280, 5140150, 4316630, 5155810, 5579520, 6049290, 4306340]; 
% WindowEndS = [6338660, 4055640, 3904260, 5649540, 4805410, 3774290, 3965830, 5484690, 7140010, 4810200,...
%     3963390, 3638780, 3798320, 3271010, 5931350, 4348540, 5140410, 4316860, 5156040, 5579800, 6049520, 4306560];

%Sampling frequency
fs = 2048; %Sample freq (Hz)

% Functions
StrideTooLongfun = @(a,b) (b-a) > mean(b-a)+2*std(b-a);
StrideTooShortfun = @(a,b) (b-a) < mean(b-a)-2*std(b-a);
Stridefun2 = @(a,b) [a(b)] ;
%% Look at the data and select an appropriate walking bout and standing 
% for d = 1:length(Participants)
%     Datapath_clinical = [DPclinical 'ROB0' (Participants{d}) '\Session' (Session{sub}) '\FreeLiving_day' (Day{d}) '.csv'];
%     data = csvread(Datapath_clinical);
% 
%     XRaw = data(:,1);
%     YRaw = data(:,2);
%     ZRaw = data(:,3);
%     data_size = length(XRaw)-1;
% 
%     %Creating time vector
%     timeRaw = 0:1:data_size;
%     timeRaw = timeRaw';
%     bpl = 1; %lower bandpass frequency (Hz)
%     bpu = 15; %upper bandpass frequency (Hz)
% 
%     X = bandpass(XRaw,[bpl bpu],fs);
%     Y = bandpass(YRaw,[bpl bpu],fs);
%     Z = bandpass(ZRaw,[bpl bpu],fs);
% 
%     %Removing samples before and after walk trial (if necessary)
%     startPoint = 1;
%     endPoint = length(X);
%     X = X(startPoint:endPoint,:);
%     Y = Y(startPoint:endPoint,:);
%     Z = Z(startPoint:endPoint,:);
%     time = timeRaw(1:(endPoint+1-startPoint),:);
% 
% 
%     % Root sum square
%     RawData = (X.^2+Y.^2+Z.^2).^(1/2);
% 
%     %Smoothing data
%     fc = 10; %butterworth cutoff frequency
%     [b,a] = butter(4,fc/(fs/2));
%     FilteredData = filter(b,a,RawData);
%     FilteredData = smoothdata(FilteredData,'sgolay','SmoothingFactor',0.2,'Degree',4);
% 
% 
% 
%     figure
%     plot(X) %bandpass filtered X-axis (vertical) data
%     hold on
%     plot(FilteredData) %combined root sum squared signal
% 
% end

%% Now manually select the heel strike moments
%Change this value one at a time to run for each participant

% files = dir(fullfile(Datapath,'*.csv'));
% 
% for i = 1:numel(files) %1:length(Participants)
%     % tic
%     if ~contains(files(i).name,'gaitEvents')
%         Datapath_clinical = fullfile(Datapath, files(i).name);
%         [~, name_file, ~] = fileparts(Datapath_clinical)

Datapath_clinical = fullfile(Datapath, 'subject_02_cond_01_run_01.csv');
[~, name_file, ~] = fileparts(Datapath_clinical)
data = readtable(Datapath_clinical);

%     fs = 100; %Sample freq (Hz)
    % Save the Standing Data
    % XRaw = data(WindowStartS(d):WindowEndS(d),1);
    % YRaw = data(WindowStartS(d):WindowEndS(d),2);
    % ZRaw = data(WindowStartS(d):WindowEndS(d),3);

    %Export Standing data for python from these participants 
    % savefile_StandingPython = [DPclinical 'ROB-' (Participants{d}) '\QtStd.csv'];
    % Standing = zeros(length(XRaw),42);
    % 
    % Standing(:,7:9) = [XRaw,YRaw,ZRaw];
    % writematrix(Standing, savefile_StandingPython)
    % clear XRaw YRaw ZRaw 
    % toc

    % Now save data split for the gait segment 
    %     XRaw = data(:,7); YRaw = data(:,8); ZRaw = data(:,9);
for c=2:size(data.Properties.VariableNames,2)
    m = data.Properties.VariableNames{1,c};

    figure;
    plot(data.time, data.(m))
    title(m)
end
    
segment_muscle = input('Introduzca el músculo que quiere usar para la segmentación', 's');
segment_muscle = strsplit(string(segment_muscle), ' ');
close all

for s = 1:size(segment_muscle,2)
    selected_muscle = string(segment_muscle{s}); 

    XRaw = data.(selected_muscle);

%     YRaw = data(WindowStartG(d):WindowEndG(d),2);
%     ZRaw = data(WindowStartG(d):WindowEndG(d),3);

%Export the walking data for python from these participants
%     savefile_Python = [DPclinical 'ROB-' (Participants{d}) '\GaitEvent.csv'];
% Gait = zeros(length(XRaw),42);

% Gait(:,7:9) = [XRaw,YRaw,ZRaw];
% writematrix(Gait, savefile_Python)
% toc

%Finding size of dataset
    data_size = size(XRaw, 1)-1;

    %Creating time vector
    timeRaw = 0:1:data_size;
    timeRaw = timeRaw';

%     R = (XRaw.^2+YRaw.^2+ZRaw.^2).^(1/2);

% Filtering
%     bpl = 1; %lower bandpass frequency (Hz)
%     bpu = 15; %upper bandpass frequency (Hz)
%   
%     X = bandpass(XRaw,[bpl bpu],fs);
%     Y = bandpass(YRaw,[bpl bpu],fs);
%     Z = bandpass(ZRaw,[bpl bpu],fs);

    Xabs = abs(XRaw);
    band = (2/fs)*1;
    [B,A] = butter(3,band,'low');
    X = filtfilt(B,A,Xabs);

%Removing samples before and after walk trial (if necessary)
    startPoint = 1;
    endPoint = size(X,1);
    X = X(startPoint:endPoint,:);
%     Y = Y(startPoint:endPoint,:);
%     Z = Z(startPoint:endPoint,:);
    time = timeRaw(1:(endPoint+1-startPoint),:);

% % Root sum square
% RawData = (X.^2+Y.^2+Z.^2).^(1/2);
% 
% %Smoothing data
% fc = 10; %butterworth cutoff frequency
% [b,a] = butter(4,fc/(fs/2));
% FilteredData = filter(b,a,RawData);
% FilteredData = smoothdata(FilteredData,'sgolay','SmoothingFactor',0.2,'Degree',4);

% plot the data, and prompt a manual segmentation where the user will
    % identify the heel strike instances 
    figure(1)
    plot(X)
    set(gcf, 'Position', get(0, 'Screensize'));
    title(selected_muscle)
    [xi,yi] = getpts(figure(1));
    % round the x indices to the nearest int and keep only non-duplicates
    
    
    if contains(selected_muscle,'ReFe') || contains(selected_muscle,'GlMe')
        percent = 0.10;
    elseif contains(selected_muscle,'GaLa')
        percent = 0.40;
    elseif contains(selected_muscle,'ErSp')
        percent = 0.60;
    elseif contains(selected_muscle,'BiFe')
        percent = 0.95;
    elseif contains(selected_muscle,'TiAn')
        percent = 0.05;
    end

    for i = 1:length(xi)-1
        cycle_duration = xi(i+1) - xi(i);
        % for i = 1:size(xi,1)
        xi(i) = xi(i) - percent * cycle_duration;
        % end
    end

    % StrideTimeTemp = round(unique(xi));
    close figure 1

    % Now we want to ensure that the points which were selected are
    % correct, and to do this we will examine the 5 points before and after
    % the identified min and clean up the true HS events

    % now search within 35/10 frames before and after identified points of 
    % the peak and chose the nearest true peak
%     StrideTime = zeros(length(StrideTimeTemp),1);
%     for pkpts = 1:length(StrideTimeTemp)
%         % for debugging, plot the current segment
% %         figure; hold on; plot(X(StrideTimeTemp(pkpts)-35:StrideTimeTemp(pkpts+1),1))
% %         plot(X(StrideTimeTemp(pkpts)-35:StrideTimeTemp(pkpts)+10,1))
%         
%         % we need to consider the data around the identified potential
%         % HS instance. We must look in a window large enough to capture the
%         % acceleration peak occuring before HS and the ensuing rise in 
%         % acceleration post HS
%         [~,lcsV] = findpeaks(-X(StrideTimeTemp(pkpts)-35:StrideTimeTemp(pkpts)+10,1));
%         % also identify the valey closest to the manual marker
%         [~,inputInd] = min(abs(lcsV-35));
%         lcsInputV = lcsV(inputInd);
% 
%         if length(lcsV) == 1
%             % note that findpeaks assumes your first cell of data is
%             % datapoint 1, so to save the correct peak location you
%             % have to remove 1 additional value from the pkpts index
%             % to result in the correct peak being saved.
%             StrideTime(pkpts) = StrideTimeTemp(pkpts)-36 + lcsV;
%         elseif length(lcsV)>1
%             % if this is the case, we need to broaden the search
%             % area and determine if there are any peaks occuring
%             % between the max and min accelerations leading up to
%             % HS as they will represent HS
% 
%             % find the acceleration peak occuring just before HS - Note
%             % there might be a local max between the peak and the HS min
%             % therefore we have to compare the values of the nearest points
%             % to find the true peak closest to the min. 
%             [~,lcsPtemp] = findpeaks(X(StrideTimeTemp(pkpts)-35:StrideTimeTemp(pkpts),1));
%             
%             if ~isempty(lcsPtemp)
%                 [~,inputPInd] = min(abs(lcsPtemp - lcsInputV));
%                 lcsInputP = lcsPtemp(inputPInd);
%                 
%                 % Now, test to see if the temp peak prior to this peak is
%                 % greater than it, and that there is no min between the 2 that
%                 % is more than 2 points from this peak (which could occur in an
%                 % inflexion point but shouldn't occur between true peaks)
%                 if inputPInd-1 > 0
%                     lcsPrevInputP = lcsPtemp(inputPInd-1);
%                     if X(StrideTimeTemp(pkpts)-36+lcsPrevInputP) > X(StrideTimeTemp(pkpts)-36+lcsInputP)
%                         if inputPInd-2 > 0
%                             %check one more level up...if it exists
%                             lcsPrevPrevInputP = lcsPtemp(inputPInd-2);
%                             if X(StrideTimeTemp(pkpts)-36+lcsPrevPrevInputP) > X(StrideTimeTemp(pkpts)-36+lcsPrevInputP)
%                                 if abs(lcsV(lcsV > lcsPrevPrevInputP & lcsV < lcsPrevInputP) - lcsPrevInputP) <= 2
%                                     lcsP = lcsPrevPrevInputP;
%                                 elseif abs(lcsV(lcsV > lcsPrevPrevInputP & lcsV < lcsPrevInputP) - lcsPrevInputP) <= 3 && (lcsInputV-lcsPrevPrevInputP) > 10 
%                                     lcsP = lcsPrevPrevInputP;
%                                 elseif abs(lcsV(lcsV > lcsPrevInputP & lcsV < lcsInputP) - lcsInputP) <= 2
%                                     lcsP = lcsPrevInputP;
%                                 elseif abs(lcsV(lcsV > lcsPrevInputP & lcsV < lcsInputP) - lcsInputP) <= 3 && (lcsInputV-lcsPrevInputP) > 10 
%                                     lcsP = lcsPrevInputP;
%                                 else
%                                     lcsP = lcsInputP;
%                                 end
%                             else
%                                 if abs(lcsV(lcsV > lcsPrevInputP & lcsV < lcsInputP) - lcsInputP) <= 2
%                                     lcsP = lcsPrevInputP;
%                                 elseif abs(lcsV(lcsV > lcsPrevInputP & lcsV < lcsInputP) - lcsInputP) <= 3 && (lcsInputV-lcsPrevInputP) > 10 
%                                     lcsP = lcsPrevInputP;
%                                 else
%                                     lcsP = lcsInputP;
%                                 end
%                             end
%                         % if the previous peak wasn't bigger than this peak (and all 
%                         % the other criteria), then test if the difference between the 
%                         % peak closest to the min identified is separated by at least
%                         % 10 frames
%                         else
%                             if abs(lcsV(lcsV > lcsPrevInputP & lcsV < lcsInputP) - lcsInputP) <= 2
%                                 lcsP = lcsPrevInputP;
%                             elseif abs(lcsV(lcsV > lcsPrevInputP & lcsV < lcsInputP) - lcsInputP) <= 3 && (lcsInputV-lcsPrevInputP) > 10 
%                                 lcsP = lcsPrevInputP;
%                             else
%                                 lcsP = lcsInputP;
%                             end
%                         end
%                     elseif (lcsInputV - lcsInputP) > 10
%                         lcsP = lcsInputP;
%                     else
%                         lcsP = lcsPtemp(end);
%                     end
%                 else
%                     lcsP = lcsInputP;
%                 end
%             else
%                 lcsP = 1;
%             end
% 
%             % now we can find any peaks which occur between the acceleration
%             % peak and the acceleration min
% 
%             % first, find the minmum acceleration following the peak
%             lcsVpotmin = lcsV(lcsV > lcsP);
%             % determine if any of these which fall after the min nearest to 
%             % the identified point are separated by a peak greater than 1/2
%             % the distance between the peak and that min. 
%             if any(lcsVpotmin(lcsVpotmin > lcsInputV))
%                 % find the peak between these points
%                 [~,lcspeakbtwnminsTemp] = findpeaks(X(StrideTimeTemp(pkpts)-36+lcsInputV : ...
%                     StrideTimeTemp(pkpts)-36+max(lcsVpotmin),1));
%                 lcspeakbtwnmins = lcspeakbtwnminsTemp + lcsInputV;
%                 if X(StrideTimeTemp(pkpts)-36+lcspeakbtwnmins) > (X(StrideTimeTemp(pkpts)-36 + lcsInputV,1) + ...
%                     (X(StrideTimeTemp(pkpts)-36 + lcsP,1) - X(StrideTimeTemp(pkpts)-36 + lcsInputV,1))/2)
% 
%                     lcsVpotmin = lcsVpotmin(~(lcsVpotmin > lcspeakbtwnmins));
%                 end
%             end
% 
%             lcsVmin = lcsVpotmin(X(StrideTimeTemp(pkpts)-35+lcsVpotmin,1) ==...
%                 min(X(StrideTimeTemp(pkpts)-35+lcsVpotmin(lcsVpotmin > lcsP),1)));
%             if lcsVmin - lcsP < 4
%                 lcsVtemp = lcsV(lcsV~= lcsVmin);
%                 lcsVmin = lcsVtemp(X(StrideTimeTemp(pkpts)-35+lcsVtemp,1) ==...
%                     min(X(StrideTimeTemp(pkpts)-35+lcsVtemp(lcsVtemp > lcsP),1)));
%             end
% 
%             % next, find any potential peaks between the min and the
%             % max accelerations in this segment
%             if isempty(lcsVmin)
%                 StrideTime(pkpts) = StrideTimeTemp(pkpts)-36 + lcsInputV;
%             else
%                 [~,lcsPotPeaksTemp] = findpeaks(X(StrideTimeTemp(pkpts)-35 + lcsP:StrideTimeTemp(pkpts)-35 + lcsVmin,1));
%                 if ~isempty(lcsPotPeaksTemp)
%                     lcsPotPeaks = lcsPotPeaksTemp + lcsP;
%                     if any(lcsPotPeaks(lcsPotPeaks > lcsInputV))
%                         lcsPotPeaks = lcsPotPeaks(~(X(StrideTimeTemp(pkpts)-35 + lcsPotPeaks,1) > ...
%                             (X(StrideTimeTemp(pkpts)-36 + lcsInputV,1) + (X(StrideTimeTemp(pkpts)-36 + lcsP,1) - ...
%                             X(StrideTimeTemp(pkpts)-36 + lcsInputV,1))/2)));
%                     end
%     
%                     % We only want to consider peaks which happen in the second
%                     % half of acceleration decrease
%                     lcsConsideredPeaks = lcsPotPeaks(lcsPotPeaks >= lcsP + (lcsVmin-lcsP)/2);
%                     if ~isempty(lcsConsideredPeaks)
%                         StrideTime(pkpts) = StrideTimeTemp(pkpts)-36 + lcsConsideredPeaks(1);
%                     else
%                         StrideTime(pkpts) = StrideTimeTemp(pkpts)-36 + lcsVmin;
%                     end
%                 elseif length(lcsV) > find(lcsV == lcsVmin)
%                     [~,lcsPotPeakAfterMin] = findpeaks(X(StrideTimeTemp(pkpts)-36 +...
%                         lcsVmin:StrideTimeTemp(pkpts)-36 + lcsV(find(lcsV == lcsVmin)+1),1));
%     
%                     % we need to check the value of this potential peak to
%                     % dicern if it occurs in the value or is in fact the next
%                     % peak following HS
%                     if X(StrideTimeTemp(pkpts)- 36 + lcsVmin + lcsPotPeakAfterMin(1)-1) > ...
%                         (X(StrideTimeTemp(pkpts)-36 + lcsVmin) + (X(StrideTimeTemp(pkpts)-36 + lcsP) - ...
%                         X(StrideTimeTemp(pkpts)-36 + lcsVmin))/2)
%                         
%                         StrideTime(pkpts) = StrideTimeTemp(pkpts)-36 + lcsVmin -1;
%                     else
%                         StrideTime(pkpts) = StrideTimeTemp(pkpts)- 36 + lcsVmin + lcsPotPeakAfterMin(1)-1;
%                     end
%                 elseif length(lcsV) == find(lcsV == lcsVmin)
%                     StrideTime(pkpts) = StrideTimeTemp(pkpts)- 36 + lcsVmin;
%                 else 
%                     lcsV = find(find(X(StrideTimeTemp(pkpts)-10:StrideTimeTemp(pkpts)+10,1) == ...
%                         min(Data(StrideTimeTemp(pkpts)-10:StrideTimeTemp(pkpts)+10,1))),1,'last');
%                     StrideTime(pkpts) = StrideTimeTemp(pkpts)-11 + lcsV;
%                 end
%             end
%             clear lcsV lcsPotPeakAfterMin lcsConsideredPeaks lcsPotPeaks lcsPotPeaksTemp lcsVmin lcsP
% %             close figure 1
%         end
% 
%     notStridesa = bsxfun(StrideTooLongfun,StrideTime(1:length(StrideTime)-1),StrideTime(2:length(StrideTime)));
%     notStridesb = bsxfun(StrideTooShortfun,StrideTime(1:length(StrideTime)-1),StrideTime(2:length(StrideTime)));
%     notStrides = or(notStridesa,notStridesb);
%     HSInd = [bsxfun(Stridefun2,StrideTime(1:end-1),~notStrides),bsxfun(Stridefun2,StrideTime(2:end),~notStrides)];
    HSInd = round(unique(xi));

    writematrix(HSInd, strcat(folder_out, name_file, '_template_', selected_muscle,'.csv'))

%     clear xi yi XRaw YRaw ZRaw X Y Z RawData FilteredData data HSInd notStridesa notStridesb notStrides StrideTime StrideTimeTemp
end


