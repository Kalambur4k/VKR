clear all;close all;


%% Параметры РСА
lamb  = 0.1;
t_per = 0.000004;          % период повторения импульсов
dev_F = 50e+7;              % девиация частоты (ширина спектра)
t_imp = 1/dev_F;           % длительность сигнала
phi   = 0;                 % начальная фаза сигнала
df    = 5*dev_F;           % частота дискретизации
dt    = 1/df;              % шаг дискретизации по времени
T     = -t_imp:dt:t_imp;   % сетка времени сигнала
T_per = 0:dt:0.5*t_per-dt;     % сетка времени периода ПИ
t_dur = 1.5;                 % Время полета (с)
t_dur_1 = 1;

h_la = 100;                % высота полета ЛА
v_la = 100;                % Скорость ЛА (м/с)
% h_la = 100;   % образец    % высота полета ЛА
% v_la = 100;   % образец    % Скорость ЛА (м/с)
% f0   = 10e9;             % Несущая (Гц)
% La   = 2;                % длина реальной антенны (м)
x_c   = 2000;              % Расстояние до центра целевой области (м).
x_c2  = 200;               % Половина ширины целевой области (м)
teta  = 60;                % ширина ДНА
% x_c   = 2000; % образец  % Расстояние до центра целевой области (м).
% x_c2  = 200;  % образец  % Половина ширины целевой области (м)
% teta  = 60;   % образец  % ширина ДНА
c     = 3e8;               % скорость света

nT_per = length(T_per);
nS     = length(T);

% Define the coordinates of reflection points
C_trg = [
  % Those 10 signals are taken as etalons for 10 filters  
    0.75*t_dur_1*v_la, 7.5, 0; %7.5 25
    0.75*t_dur_1*v_la, 24.5, 0; %24.5 75
    0.75*t_dur_1*v_la, 37.5, 0; %37.5 125
    0.75*t_dur_1*v_la,  52.5, 0;
    0.75*t_dur_1*v_la, 67.5, 0;
    0.75*t_dur_1*v_la, 82.5, 0;
    0.75*t_dur_1*v_la,  97.5, 0;
    0.9*t_dur_1*v_la, 112.5, 0;
    0.9*t_dur_1*v_la, 127.5, 0;
    0.9*t_dur_1*v_la,  142.5, 0;
% Script takes first n signals as etalons for n filters, where n is number of filters (variable flt_num)

0.45*t_dur_1*v_la, 80, 0; % tail
0.47*t_dur_1*v_la, 80, 0;
0.49*t_dur_1*v_la, 80, 0;
0.51*t_dur_1*v_la, 80, 0;
0.53*t_dur_1*v_la, 80, 0; % start of wings
0.55*t_dur_1*v_la, 82, 0;
0.57*t_dur_1*v_la, 84, 0;
0.59*t_dur_1*v_la, 86, 0;
0.61*t_dur_1*v_la, 88, 0; % end of wings
0.61*t_dur_1*v_la, 86, 0;
0.61*t_dur_1*v_la, 84, 0;
0.61*t_dur_1*v_la, 82, 0;
0.63*t_dur_1*v_la, 80, 0; % body
0.65*t_dur_1*v_la, 80, 0;
0.67*t_dur_1*v_la, 80, 0;
0.69*t_dur_1*v_la, 80, 0;
0.71*t_dur_1*v_la, 80, 0; % start of wings
0.73*t_dur_1*v_la, 78, 0;
0.75*t_dur_1*v_la, 76, 0;
0.77*t_dur_1*v_la, 74, 0;
0.79*t_dur_1*v_la, 72, 0; % end of wings
0.79*t_dur_1*v_la, 74, 0;
0.79*t_dur_1*v_la, 76, 0;
0.79*t_dur_1*v_la, 78, 0;
0.81*t_dur_1*v_la, 80, 0; % nose
0.79*t_dur_1*v_la, 80, 0; % closing the shape
0.77*t_dur_1*v_la, 80, 0;
0.75*t_dur_1*v_la, 80, 0;
0.73*t_dur_1*v_la, 80, 0;
0.71*t_dur_1*v_la, 80, 0;
0.69*t_dur_1*v_la, 80, 0;
0.67*t_dur_1*v_la, 80, 0;
0.65*t_dur_1*v_la, 80, 0;
0.63*t_dur_1*v_la, 80, 0;
0.61*t_dur_1*v_la, 80, 0;
0.59*t_dur_1*v_la, 80, 0;
0.57*t_dur_1*v_la, 80, 0;
0.55*t_dur_1*v_la, 80, 0;
0.53*t_dur_1*v_la, 80, 0;
0.51*t_dur_1*v_la, 80, 0;
0.49*t_dur_1*v_la, 80, 0;
0.47*t_dur_1*v_la, 80, 0;
0.45*t_dur_1*v_la, 80, 0; % back to tail

0.12*t_dur_1*v_la, 145, 0;
0.78*t_dur_1*v_la, 67, 0;
0.34*t_dur_1*v_la, 89, 0;
0.56*t_dur_1*v_la, 123, 0;
0.91*t_dur_1*v_la, 56, 0;
0.67*t_dur_1*v_la, 78, 0;
0.29*t_dur_1*v_la, 134, 0;
0.85*t_dur_1*v_la, 90, 0;
0.47*t_dur_1*v_la, 115, 0;
0.73*t_dur_1*v_la, 45, 0;
0.58*t_dur_1*v_la, 67, 0;
0.39*t_dur_1*v_la, 123, 0;
0.96*t_dur_1*v_la, 78, 0;
0.52*t_dur_1*v_la, 90, 0;
0.81*t_dur_1*v_la, 112, 0;
0.63*t_dur_1*v_la, 145, 0;
0.37*t_dur_1*v_la, 67, 0;
0.94*t_dur_1*v_la, 123, 0;
0.51*t_dur_1*v_la, 149, 0;
0.72*t_dur_1*v_la, 78, 0;
0.86*t_dur_1*v_la, 134, 0;
0.43*t_dur_1*v_la, 90, 0;
0.59*t_dur_1*v_la, 112, 0;
0.77*t_dur_1*v_la, 145, 0;
0.33*t_dur_1*v_la, 67, 0;
0.88*t_dur_1*v_la, 123, 0;
0.46*t_dur_1*v_la, 26, 0;
0.71*t_dur_1*v_la, 98, 0;
0.35*t_dur_1*v_la, 134, 0;
0.92*t_dur_1*v_la, 120, 0;


];

flt_num = 10; % Number of filters in range [0;150] ( etalons in C_trg must be corrected according to number of filters)
% If not corrected, filtering will give wrong result


% Basic parameters
dx = 0.05;
dr = 3e8*dt*0.5;
nx = round(t_dur*v_la/dx);
ny = nT_per;

X_la = 0:dx:dx*(nx-1);
C_la = zeros(1, 3);
C_la(2) = 0;
C_la(3) = h_la;

% Initialize etalon signals and Sig_Trek
Sig_Etal = cell(1,flt_num); %
Sig_Trek = cell(1,flt_num); %
for i=1:flt_num %
    Sig_Etal{i} = [];
    Sig_Trek{i} = zeros(ny, nx);
end

% Initialize categorization lists
C_dist = cell(1,flt_num); %
for i=1:flt_num
    C_dist{i} = [];
end

% Categorize signals and apply to corresponding etalons
for trg_idx = 1:size(C_trg, 1)
    x0 = C_trg(trg_idx, 1);
    y0 = C_trg(trg_idx, 2);
    C_trg_single = [x0, y0, 0];
    
    % Categorize based on distance
  if y0 >= 1 && y0 <= 150
    dist_cat = ceil(y0 / (150 / flt_num)); %<- 150 / number of filters
else
    dist_cat = 0;


    end
    C_dist{dist_cat} = [C_dist{dist_cat}; trg_idx];

    % Loop through each
    indx_1 = 0;
    for ix = 1:nx
        C_la(1) = (ix - 1) * dx;
        r_trg = norm(C_trg_single - C_la);
        alf = atand(abs((C_la(1) - C_trg_single(1))) / C_trg_single(2));
        ku = exp(-0.7 * (2 * alf / 30)^2)^4;
        t_zad = 2 * r_trg / 3e8;
        phi = (4 * pi / lamb) * r_trg;
        ST = F_Insert_Signal(t_zad,T_per,dt,dev_F,ku,phi);
        
        % Add to corresponding Sig_Trek
        Sig_Trek{dist_cat}(:,ix) = Sig_Trek{dist_cat}(:,ix) + ST';
        
        % Add to corresponding Sig_Etal
        if ku > 0.1
            indx_1 = indx_1 + 1;
            ST = F_Insert_Signal(t_zad/20,T_per,dt,dev_F,ku,phi);
            nST = length(ST);
            n0  = round(nST/20);
            Sig_Etal{dist_cat}(1:n0,indx_1) = ST(1:n0)';
        end
    end
end

% Display Sig_Trek
for i=1:flt_num %
    figure(i+1); imshow(abs(Sig_Trek{i}),[]); 
end

Combined_Signal = zeros(ny, nx);

% Filter each Sig_Trek and add to Combined_Signal
for i=1:flt_num %
    Sig_Etalon = Sig_Etal{i};
    figure(3); imshow(abs(Sig_Etalon),[]); 

    A = Sig_Etalon;
    A = A / max(abs(A(:)));
    A = imrotate(A,180);
    A = conj(A);
    [mB, nB] = size(A);
    C1 = [A zeros(mB, nx-nB)];
    C2 = [C1; zeros(ny-mB, nx)];
    FS = fft2(Sig_Trek{i});
    FC = fft2(C2);
    FH = FS.* FC;
    H = ifft2(FH);

    Combined_Signal = Combined_Signal + H;
end


% Display the final result for the combined signal
figure(50); imshow(abs(Combined_Signal),[]);
title('Final Combined Filter Result');