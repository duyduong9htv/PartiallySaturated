

for ii = 2:10:1342
    time_series = []; 
    for jj = ii:(ii+9)
        load_command = ['load pool_test_nov_' num2str(jj) '.txt']; 
        eval(load_command); 
        assign_command = ['temp = pool_test_nov_' num2str(jj) ';']; 
        eval(assign_command); 
        time_series = [time_series; temp(:, 2)]; 
    end 
    [S, F, T, P] = spectrogram(time_series, 128, 56, 2048, 2000); 
    figure(1); clf; 
    imagesc(T, F, 10*log10(P)); 
    setFont(20);setFigureAuto; 
    colorbar; 
    axis xy; 
    ti = ['pool test nov 2 sample ' num2str(ii)];
    title(ti); 
    ylabel('Frequency (Hz)'); 
    xlabel('Time (seconds)'); 
    save_command = ['print -dtiff pool_test_nov_2_' num2str(ii) '.tif']; 
    eval(save_command);     
end

%% STFT to look at amplitude distribution 

Pow = []; 
Sf_all = []; 
for ii = 1:(0.1*2000):(length(time_series)-(0.1*2000))
    st = time_series(ii:(ii+200)); 
    Sf = fft(st); 
    Sf_all = [Sf_all; Sf(45)]; 
    Pow = [Pow 20*log10(abs(Sf(45)))];     
end


%% phase distribution 

%data check 
prefix_str = 'pool_test_nov_'; 
fileno = 2; 
s = []; 
for fileno = 500:1:1500
    s1 = load([prefix_str num2str(fileno) '.txt']); 
    t = s1(:, 1); %time vector 
    s1 = s1(:, 2); %signal
    s = [s; s1]; 
end 


figure; plot(1/2000*[1:1:length(s)], s); 


% 

t_start = 2.5; 
T = 2; 
step = 4; 
fs = 2000; 
fc = 550; 
ph = []; 
for t1 = t_start:step:54.5
    u = s(t1*fs:(t1+T)*fs); 
    uh = hilbert(u);
    Uf = fft(uh); 
    ph = [ph; mod(phase(Uf(1100)), 2*pi)]; 
    figure(10); hold on; 
    polar(mod(phase(Uf(1100)), 2*pi), abs(Uf(1100)), 'go'); 
    pause(1);     
end

%design matched filter to determine arrival time 
t = linspace(0, T, length(Uf)); 
w = [ones(1, fs) zeros(1, length(t)-fs)]; 
h = cos(2*pi*fc*t); 
h = h.*w; 
hu = hilbert(h); 
H = fft(hu); 

MF = H.*transpose(Uf); 


%% 
cos_phi = []
for ind = 2700:2000:(length(s_demod) - 2000)
    cos_phi = [cos_phi; 2*mean(s_demod(ind:ind+2000))]; 
end


%% testing the demodulating method 


s3  = []; 
T = 4;
fs = 2000; 
fc = 550; 
t = linspace(0, T, T*fs); 
w = [ones(1, fs) zeros(1, length(t)-fs)]; 
h = sin(2*pi*fc*t).*w; 

for k = 1:10
    s3 = [s3 h]; 
end


%reference signal 
t_ref = 1/2000*[1:1:length(s)]; 
fc = 550; 
s_ref = sin(2*pi*fc*t_ref); 
figure(3); hold on; 
plot(t_ref, 1/3*s_ref, '-r'); 

% check frequency content 

%data check 
prefix_str = 'pootest_nov_'; 
s = []; 
peaks = []; 
for fileno = 382:1:397
    s1 = load([prefix_str num2str(fileno) '.txt']); 
    N = length(s1(:, 2)); 
    disp(N);
    s_h = hilbert(s1(:, 2)); 
    Sf = fft(s_h);
    ave_energy = mean(abs(Sf).^2);
    max_energy = max(abs(Sf).^2); 
    if (10*log10(max_energy) - 10*log10(ave_energy) > 3)
        max_ind = find(abs(Sf(:)).^2 == max_energy);
        freq_peak = (max_ind/N)*2000; 
        peaks = [peaks freq_peak]; 
    end
    
    s = [s; s1(:, 2)]; 
end 







