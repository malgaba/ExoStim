% Compute speed
data = load('/Users/neuralrehabilitationgroup/PycharmProjects/DATA_ExoStim/EXOSTIM_DATA/RESULTS_ExoStim.mat');

SUBJECTS = [1, 3, 4];
sf = 2048;
speeds = table();

% Distance/Time (m/s)

% Velocidad subject_01_cond_01
v111 = 76 / (data.RESULTS(1).indexesGait(2)/sf - data.RESULTS(1).indexesGait(1)/sf);
v112 = 82 / (data.RESULTS(1).indexesGait(4)/sf - data.RESULTS(1).indexesGait(3)/sf);
v113 = 74 / (data.RESULTS(1).indexesGait(6)/sf - data.RESULTS(1).indexesGait(5)/sf);

% Velocidad subject_03_cond_01
v311 = 22.5 / (data.RESULTS(3).indexesGait(2)/sf - data.RESULTS(3).indexesGait(1)/sf);
v312 = 30 / (data.RESULTS(3).indexesGait(4)/sf - data.RESULTS(3).indexesGait(3)/sf);
v313 = 33 / (data.RESULTS(3).indexesGait(6)/sf - data.RESULTS(3).indexesGait(5)/sf);

% Velocidad subject_04_cond_01
v411 = 55 / (data.RESULTS(4).indexesGait(2)/sf - data.RESULTS(4).indexesGait(1)/sf);
v412 = 65.2 / (data.RESULTS(4).indexesGait(4)/sf - data.RESULTS(4).indexesGait(3)/sf);
v413 = 61.5 / (data.RESULTS(4).indexesGait(6)/sf - data.RESULTS(4).indexesGait(5)/sf);

% Velocidad subject_01_cond_02
v121 = 33 / (data.RESULTS(1).indexesExo(2)/sf - data.RESULTS(1).indexesExo(1)/sf);
v122 = 33 / ((data.RESULTS(1).indexesExo(4)/sf - data.RESULTS(1).indexesExo(3)/sf) + (data.RESULTS(1).indexesExo(8)/sf - data.RESULTS(1).indexesExo(7)/sf));
v123 = 32.25 / (data.RESULTS(1).indexesExo(6)/sf - data.RESULTS(1).indexesExo(5)/sf);

% Velocidad subject_03_cond_02
v321 = 31 / (data.RESULTS(3).indexesExo(2)/sf - data.RESULTS(3).indexesExo(1)/sf);
v322 = 22.5 / (data.RESULTS(3).indexesExo(4)/sf - data.RESULTS(3).indexesExo(3)/sf);
v323 = 29.25 / (data.RESULTS(3).indexesExo(6)/sf - data.RESULTS(3).indexesExo(5)/sf);

% Velocidad subject_04_cond_02
v421 = 35 / (data.RESULTS(4).indexesExo(2)/sf - data.RESULTS(4).indexesExo(1)/sf);
v422 = 37 / (data.RESULTS(4).indexesExo(4)/sf - data.RESULTS(4).indexesExo(3)/sf);
v423 = 33 / (data.RESULTS(4).indexesExo(6)/sf - data.RESULTS(4).indexesExo(5)/sf);

speed = {'COND_01', v111, v112, v113, v311, v312, v313, v411, v412, v413, ...
          'COND_02', v121, v122, v123, v321, v322, v323, v421, v422, v423}

% speeds = ['COND_01', v111, v112, v113, v311, v312, v313, v411, v412, v413, 'COND_02', v121, v122, v123, v321, v322, v323, v421, v422, v423]
