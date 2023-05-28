clc
close all
clear
% Github Token Acesss ccccc

%% Definição de parâmetros (Etapa de escolha dos parâmetros a serem adotados no simulador.)
% Parâmetros Físicos (Valores usuais na física/engenharia)
c = physconst('LightSpeed'); % velocidade da luz

% Parametros do Rádio (Valores customizaveis do rádio escolhido)
freq = 2.4e9; % frequência em Hz
lambda = c/freq; % comprimento de onda em metros
Pt = 10; % potência de transmissão em dBm
B = 20e6; % largura de banda em Hz

% Parametros das Antenas (Valores customizaveis das antenas)
Gt = 0; % ganho da antena de transmissão em dBi
Gr = 0; % ganho da antena de recepção em dBis

% Parametros de ruído (Valores customizaveis de ruído)
N0 = -174 + 10*log10(B); % densidade espectral de ruído em dBm/Hz
Pn = N0 + 10*log10(B); % potência de ruído em dBm
% Parametros de canal (Valores customizaveis de canal)
channelName = 'TraditionalWirelessChannel'; % Opçãões:[TraditionalWirelessChannel,TRGRWirelessChannel,COST231WirelessChannel]
add_fast_fading = true; % Adicionar fast fading no canal.
add_shadowing =true; % Adicionar sombreamento no canal;
add_delay_profile =true; % Adicionar perfil de atraso no canal;
add_doppler_effect =true; % Adicionar atenuação por efeito dopple;
d0 = 1; % distância de referência em metros
alpha = 2; % expoente de perda de percurso []
sigma_sf = 2; % desvio padrão do sombreamento log-normal em dB [padrão:4, 2.49 para LTE basiado em Analysis and Modeling for Train-Ground Wireless Wideband Channel of LTE on High-Speed Railway]
sigma_ff = 3; % desvio padrão do fast fading em dB [padrão:8]
SF= 0; % Valor inicial do sombreamento em dB
FF= 0; % Valor inicial do fast fading em dB
DP= 0; % Valor inicial do atenuação de atraso de perfil (delay profile) em dB
% Especifico para COST-231
hte = 5; %Altura da antena Emissora (No caso o Ponto de Acesso) (em metros) Recomenda-se que esteja na gama [4,50] metros
hre = 3; %Altura efectiva da antena no terminal movél (No caso o trem) (em metros) Recomenda-se que esteja na gama [1,3] metros
ws = 12; %Largura da rua (expressa em metros), onde o terminal móvel se encontra; admite-se que este se localiza no centro da via
wb =10; % Distancia média entre edificios, marcada entre os pontos centrais dos mesmos, para efeitos do calculo da média
hb= 20; % Altura média dos edificios; neste parametros, inclui-se a contribuição dos telhados
zona =1; % Tipo de ambiente urbano em causa: 1 - Grandes Centros Urbanos 0 - Cidades de dimensões razoáveis ou areas Sub-urbanas.
% O seguinte artigo é um survay para canais de trem. Channel Measurements and Models for High-Speed Train Communication Systems: A Survey
% Parâmetros de simulação (Valores relacionados a caracteristicas do simulador)
T = 10; % duração da simulação em segundos
fs = 100; % frequência de amostragem em Hz
t = linspace(0, T, T*fs); % vetor de tempo.

% Parâmetros do Trem (Valores customizaveis do trem)
v = 180/3.6; % velocidade do trem em metros por segundo
ht = 5; % Altura da antena do trem (m)
% Parâmetros de infraestrutura (Valores relacioandos a infraestrutura do cenário)
AP_pos = [1,0;50,0;100,0]; % posições dos pontos de acesso em metros. Ex:[Xap1,Yap1;Xap2,Yap1], onde Xap1 e Yap1 corresponde a posição da Ponto de Acesso 1 no eixo X e Y, respectivamente.
hr = 1.5; % Altura da antena do AP (m)
%% Seção de Movimentação do Trem
pos_trem = v*t'; % Posição do trem em linha reta.
for i = 1:3
    d(i,:) = sqrt(sum((AP_pos(i,:) - reshape([pos_trem' zeros(1,length(v*t)) ], length(t), []) ).^2, 2))'; % distâncias dos pontos de acesso ao trem
end
%% Seção dos Canais
% Parte referente aos modelos de canais sem fio. (Obs: A saída do canal pode ser em dB, linear(watts) ou amplitude. Deve-se verificar bem isso para evitar confunsão adiante.)
switch channelName
    case 'TraditionalWirelessChannel'
        % O trecho de código apresentado é utilizado para modelar a perda de
        % percurso, o sombreamento log-normal e o desvanecimento rápido (fast fading) em um canal sem fio.
        % Ele é uma representação matemática dos efeitos que um sinal de comunicação sofre ao se propagar
        % através do ambiente sem fio, levando em consideração a atenuação do sinal devido à distância percorrida,
        % variações aleatórias devido ao sombreamento do ambiente e variações aleatórias devido ao desvanecimento rápido.
        PL = -10*alpha*log10(d./d0); % perda de percurso em dB
        if(add_shadowing)
            SF = sigma_sf*randn(3, length(t)); % sombreamento log-normal em dB
        end
        if(add_fast_fading)
            FF = sigma_ff*randn(3, length(t));  % fast fading em dB
        end
        if(add_delay_profile)
            DP = DelayProfile(c, lambda, d0, alpha, Gr , Gt);  % Atenuação por perfil de atraso
        end
   
        channel = PL + SF + FF + DP ; % Canal sem fio em dB
    case 'TRGRWirelessChannel'
        % Cálculo da perda de percurso usando o modelo Two-Ray Ground Reflection
        d0 = (ht - hr)^2 / (16 * hr);  % distância de referência em metros]
        [lin,col]=size(d);
        L = zeros(lin,col);
        for l = 1:lin
            for c=1:col
                if d(l,c) <= d0
                    L(l,c) = -20*log10((4*pi*d(l,c)*hr*ht) / lambda^2);
                else
                    L(l,c) = -20*log10((4*pi*d(l,c)*hr*ht) / lambda^2) - 10*alpha*log10(d(l,c)/d0);
                end
            end
        end
        if(add_shadowing)
            SF = sigma_sf*randn(3, length(t)); % sombreamento log-normal em dB
        end
        if(add_fast_fading)
            FF = sigma_ff*randn(3, length(t));  % fast fading em dB
        end
        if(add_delay_profile)
            DP = DelayProfile(c, lambda, d0, alpha, Gr , Gt);  % Atenuação por perfil de atraso
        end
    
        channel = PL + SF + FF + DP; % Canal sem fio em dB
    case 'COST231WirelessChannel'
        L = COST231WirelessChannel(d, freq, hte, hre, ws, wb, hb, zona);
        if(add_shadowing)
            SF = sigma_sf*randn(3, length(t)); % sombreamento log-normal em dB
        end
        if(add_fast_fading)
            FF = sigma_ff*randn(3, length(t));  % fast fading em dB
        end
        if(add_delay_profile)
            DP = DelayProfile(c, lambda, d0, alpha, Gr , Gt);  % Atenuação por perfil de atraso
        end
   
        channel = PL + SF + FF + DP; % Canal sem fio em dB
    otherwise
        disp('O canal a ser executado não exite. Verifique se foi chamado corretamente.')
end

%% Seção dos Sinais de Informação (Bits de Informação)
% Geração de bits aleatórios
%nbits = 1e6;         % Número total de bits transmitidos
%bits = randi([0 1], 3, nbits);
% Modulação dos bits
%M = 64; % Modulação 64-QAM Se encaixa com MCS 7 Abaixo da taxa de transmissão
%data = qammod(bits, M);
%% Seção da Recepção da Informação pelo AP (Bits de Informação)

Pr = zeros(3, length(t)); % potência de recepção em dBm
for i = 1:3
    Pr(i, :) = Pt + Gt + Gr + channel(i, :) + Pn;
end
%% Seção Metricas de Desempenho da Rede
% Cálculo da cobertura
threshold = -80; % threshold de potência em dBm
coverage = (Pr > threshold); % cobertura em porcentagem
% Cálculo da taxa de transmissão
MCS = 7; % esquema de modulação e codificação
Tb = 1/(MCS*1000*2); % tempo de bit em segundos
% Taxa de transmissão em Mbps
R = zeros(3, length(t));
SNRdB = Pr-Pn; % dB = dBm -dBm

dataload = load('Wifi_80211g.mat', 'SNRb','throughput_of_file');
throughputFile= dataload.('throughput_of_file');
SNRdBFile = dataload.('SNRb');
for l = 1:3
    %R(i, :) = abs(MCS*B*log2(1 + 10.^(Pr(i, :)/10)/N0)/1e6); % Analisar melhor
   % Encontrar os índices dos valores iguais ou próximos
for i = 1:length(SNRdB)
    [~, indice] = min(abs(SNRdBFile - SNRdB(l,i)));
    indicesEncontrados(i) = indice;
end
    R(l, :) = throughputFile(indicesEncontrados);
end
avg_R = mean(R, 2); % taxa de transmissão média em Mbps

% Cálculo da latência
latency = d/v + T/2;

% Cálculo do outage
%outage = sum(Pr < threshold, 2)/length(t); % outage em porcentagem
outage = (Pr < threshold); % outage em porcentagem
% Cálculo da QoE
%QoE = 1./(1 + exp(-0.1*(avg_R - 5)));
QoE = 1./(1 + exp(-0.1.*(R - 5)));
% Cálculo da QoS
%QoS = sum(coverage.*QoE);
QoS = coverage.*QoE;
% Cálculo da SNR
SNR = Pr - Pn;

%
% Transmissão dos dados
Pr_aux = Pr;
Pr_aux(Pr_aux == -inf) = 0; % Retira -Inf da potencia de recepção e substitui por 0.
Pr_aux(Pr_aux == inf) = 0; % Retira +Inf da potencia de recepção e substitui por 0.
%txsignal = awgn(abs(mean(sqrt(Pr_aux),2)).*data, Pn, 'measured');
% Recebimento dos dados
%rxsignal = awgn(txsignal, Pn, 'measured');
% Demodulação dos dados
%rxdata = qamdemod(rxsignal, M);
% Cálculo da BER
%ber = sum(xor(bits, rxdata),2)/nbits;
% for f=1:1:length(t)
%     for i=1:1:3
% txsignal = awgn(abs(sqrt(Pr_aux(i,f))).*data(i,:), Pn, 'measured');
% % Recebimento dos dados
% rxsignal = awgn(txsignal, Pn, 'measured');
% % Demodulação dos dados
% rxdata = qamdemod(rxsignal, M);
% % Cálculo da BER
% ber(i,f) = sum(xor(bits(1,:), rxdata),2)/nbits;
%     end
% end

for i=1:1:3
    ber(i,:) = 1/ (2^(6/5) - 1) * erfc(sqrt((10.^(SNR(i,:)./10))*5/6/10));
end

% Visualização dos resultados
disp('Cobertura:');
disp(coverage);
disp('Taxa de transmissão média (Mbps):');
disp(avg_R);
disp('Latência (s):');
disp(latency);
disp('Outage:');
disp(outage);
disp('QoE:');
disp(QoE);
disp('QoS:');
disp(QoS);
disp('SNR (dBm):');
disp(SNR);

% Visualização da potência de recepção ao longo do tempo
figure;
plot(t, Pr(1, :), 'r', t, Pr(2, :), 'g', t, Pr(3, :), 'b');
xlabel('Tempo (s)');
ylabel('Potência de recepção (dBm)');
legend('AP1', 'AP2', 'AP3');

% Visualização da cobertura ao longo do tempo
figure;
plot(t, coverage(1,:), 'r', t, coverage(2,:), 'g', t, coverage(3,:), 'b');
xlabel('Tempo (s)');
ylabel('Cobertura (%)');
legend('AP1', 'AP2', 'AP3');

% Visualização da taxa de transmissão ao longo do tempo
figure;
plot(t, R(1, :), 'r', t, R(2, :), 'g', t, R(3, :), 'b');
xlabel('Tempo (s)');
ylabel('Taxa de transmissão (Mbps)');
legend('AP1', 'AP2', 'AP3');

% Visualização da latência ao longo do tempo
figure;
plot(t, latency(1,:), 'r', t, latency(2,:), 'g', t, latency(3,:), 'b');
xlabel('Tempo (s)');
ylabel('Latência (s)');
legend('AP1', 'AP2', 'AP3');

% Visualização do outage ao longo do tempo
figure;
plot(t, 100*outage(1,:), 'r', t, 100*outage(2,:), 'g', t, 100*outage(3,:), 'b');
xlabel('Tempo (s)');
ylabel('Outage (%)');
legend('AP1', 'AP2', 'AP3');

% Visualização da QoE ao longo do tempo
figure;
plot(t, QoE(1,:), 'r', t, QoE(2,:), 'g', t, QoE(3,:), 'b');
xlabel('Tempo (s)');
ylabel('QoE');
legend('AP1', 'AP2', 'AP3');

% Visualização da QoS ao longo do tempo
figure;
plot(t, QoS(1,:), 'r', t, QoS(2,:), 'g', t, QoS(3,:), 'b');
xlabel('Tempo (s)');
ylabel('QoS');
legend('AP1', 'AP2', 'AP3');

% Visualização da SNR ao longo do tempo
figure;
plot(t, SNR(1, :), 'r', t, SNR(2, :), 'g', t, SNR(3, :), 'b');
xlabel('Tempo (s)');
ylabel('SNR (dBm)');
legend('AP1', 'AP2', 'AP3');

% Visualização da BER ao longo do tempo
figure;
plot(t, ber(1,:), 'r', t, ber(2, :), 'g', t, ber(3, :), 'b');
xlabel('Tempo (s)');
ylabel('BER ');
legend('AP1', 'AP2', 'AP3');

%% Plot da constelação e do sinal recebido
% figure;
% scatterplot(rxsignal(1, :));
% title('Constelação do sinal recebido');
% xlabel('I');
% ylabel('Q');
%
% figure;
% plot(real(rxsignal(1, :)));
% hold on;
% plot(imag(rxsignal(1, :)));
% title('Sinal recebido');
% xlabel('Tempo');
% ylabel('Amplitude');
% legend('Parte real','Parte imaginária');

%% Print da BER
%fprintf('Bit Error Rate (BER): %g\n',ber);
