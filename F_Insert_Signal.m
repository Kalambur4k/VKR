function ST = F_Insert_Signal(t_zad,T_per,dt,dev_F,ampl,phi)
%% Программа вставки эталонной функции отклика цели (ФОЦ)
%  в период повторения импульсов (ППИ) с заданным временем запаздывания.
%  ФОЦ - функция sinc(), соответствует комплексной АКФ ЛЧМ сигнала.
%  Важно! Вставка сигнала выполняется с учетом дискретности оцифровки,
%         т.е. момент появления сигнала м.б. расположен между 
%         дискретными отсчетами времени.
%  t_zad  - время задержки сигнала (с)
%  T_per  - сетка времени одного ППИ
%  dt     - шаг дискретизации по времени (с)
%  dev_F  - ширина спектра сигнала (Гц)
%  phi    - начальная фаза сигнала (рад)
%       
%       ST - отсчеты в одном ППИ
t_imp = 1/dev_F;           % длительность сигнала
T     = -t_imp:dt:t_imp;   % сетка времени сигнала
nS    = length(T);
%% Вставка сигнала в период ПИ
nT_per = length(T_per);
n_zad  = floor(t_zad/dt);
dt_zad = t_zad - n_zad*dt;
S   = sinc(dev_F.*(T-dt_zad)).*exp(-1i*(pi*dev_F.*(T-dt_zad)+phi));% АКФ 
S   = S.*exp(-((1.2/t_imp)*(T-dt_zad)).^2);  % весовое окно
ST  = zeros(1,nT_per);
nS2 = floor(nS/2);
n1 = n_zad - nS2;
n2 = n_zad + nS2;
if (n1 > 1) && (n2 <= nT_per)
   ST(n_zad-nS2:n_zad+nS2) = ampl*S;
end
% figure(5); plot(T_per,real(ST),'b',T_per,imag(ST),'r'); grid on;
% figure(6); stem(real(ST)); grid on;

end