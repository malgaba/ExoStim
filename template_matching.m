% TEMPLATE MATCHING
time = data.time;
cleaned_emg = data.TiAn_l;


figure;
plot(time, cleaned_emg);
title('Selecciona manualmente la plantilla');
xlabel('Tiempo');
ylabel('Amplitud');

fprintf('Haz clic en la gráfica para seleccionar los puntos de la plantilla.\n');
fprintf('Presiona Enter cuando hayas terminado.\n');
[x, ~] = ginput;

% Convertir las coordenadas x a índices de muestras
indices_muestras = round(x);

% Extraer la plantilla de la señal EMG
plantilla = cleaned_emg(indices_muestras(1):indices_muestras(end));

% Mostrar la plantilla seleccionada
figure;
plot(plantilla);
title('Plantilla seleccionada');
xlabel('Muestras');
ylabel('Amplitud');