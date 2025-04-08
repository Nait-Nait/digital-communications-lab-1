clc;
close all;
clear all;

% parámetros de configuración
fm = 100000; 
tm = 1/fm;
ls = 200; 
f_c = 1000; 
f_s = 5000; 
t_s = 1/f_s; 
tau = 0.5*t_s; 
d = tau/t_s; 

% vectores
t = (0:ls-1)*tm;
m_t = sin(2*pi*f_c*t);

% auxiliaries
r = floor(t_s/tm);
s = floor(tau/tm);
disp(r)

% muestreo natural
s_nat = zeros(1,length(t));
for i=1:length(m_t)
    if mod(i,r)==0
        s_nat(i:i+s) = 1;
    end
end
s_nat = s_nat(1:length(t));
m_t_nat = m_t.*s_nat;


% muestreo instantaneo
m_t_inst = zeros(1,length(t));
for i=1:length(m_t)
    if mod(i,r)==0
        m_t_inst(i:i+s) = m_t(i);
    end
end

m_t_inst = m_t_inst(1:length(t));

figure;
plot(t, m_t); 
hold on;
plot(t, m_t_nat, 'r'); 
plot(t, m_t_inst, 'g'); 
grid on;

xlabel('Tiempo (s)');
ylabel('Amplitud');
legend;
title('Comparación de Muestreos: Natural vs. Instantáneo');


% Fourier Transform
M_t = fft(m_t);
M_t_nat = fft(m_t_nat);
M_t_inst = fft(m_t_inst);

% Frequency axis for plotting
f_axis = (0:(length(t) - 1)) * (1 / (ls * tm));

% Create the first figure with original, natural, and instantaneous signals
figure;
subplot(3, 1, 1);
plot(t, m_t);
title('Señal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(3, 1, 2);
plot(t, m_t_nat);
title('Muestreo Natural');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(3, 1, 3);
plot(t, m_t_inst);
title('Muestreo Instantáneo');
xlabel('Tiempo (s)');
ylabel('Amplitud');

% Create the second figure with Fourier Transforms
figure;
subplot(3, 1, 1);
plot(f_axis, abs(M_t));
title('Transformada de Fourier Señal Original');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');

subplot(3, 1, 2);
plot(f_axis, abs(M_t_nat));
title('Transformada de Fourier Muestreo Natural');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');

subplot(3, 1, 3);
plot(f_axis, abs(M_t_inst));
title('Transformada de Fourier Muestreo Instantáneo');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');

% PCM Parameters
bit_depth = 8; 
pcm_levels = 2^bit_depth;

% Quantize the instantaneous signal using PCM
pcm_signal_inst = round((m_t_inst + 1) * (pcm_levels - 1) / 2);
grid on;

% Normaliza las señales para que estén en el mismo rango de amplitud (0 a 1)
m_t_norm = (m_t - min(m_t)) / (max(m_t) - min(m_t));
m_t_inst_norm = (m_t_inst - min(m_t_inst)) / (max(m_t_inst) - min(m_t_inst));
pcm_signal_inst_norm = (pcm_signal_inst - min(pcm_signal_inst)) / (max(pcm_signal_inst) - min(pcm_signal_inst));

% Calcular el error de cuantización para la señal PAM cuantificada (PCM)
quantization_error_inst = m_t_inst - ((2 * pcm_signal_inst / (pcm_levels - 1)) -1);

% Crear una figura para mostrar todas las señales en un mismo gráfico
figure;
plot(t, m_t_norm, 'b', 'LineWidth', 1.5); 
hold on;
plot(t, m_t_inst_norm, 'r', 'LineWidth', 1.5); 
stem(t, pcm_signal_inst_norm, 'g', 'Marker', 'o', 'LineWidth', 1.5);

xlabel('Tiempo (s)');
ylabel('Amplitud Normalizada');
title('Señal Original, Señal PAM Instantánea y Señal PAM Cuantificada (PCM)');
legend('Señal Original', 'Señal PAM Instantánea', 'Señal PAM Cuantificada (PCM)');
grid on;

% Crear una nueva figura para mostrar el error de cuantización
figure;
plot(t, quantization_error_inst, 'k--', 'LineWidth', 1.5);

xlabel('Tiempo (s)');
ylabel('Error de Cuantización');
title('Error de Cuantización para la Señal PAM Cuantificada (PCM)');
grid on;
