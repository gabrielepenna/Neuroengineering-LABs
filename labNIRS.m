%LAB NIRS
clc
clear all

fs = 10;
fNy = fs/2;

%carico i pazienti da S01 a S08
for k=1:8
    formatSpec = 'S0%d.mat';
    A1 = k;
    str = sprintf(formatSpec,A1);
    
    load (str)
    n_run=length(data);

%carico i RUN
for j=1:n_run
channels = data{1,j}.X;
oxy_Hb = channels(:,1:52);
deoxy_Hb = channels(:,53:104);
total_Hb = channels(:,105:156);

%PUNTO a)--------------
%scelgo canale 48 per APFC
APFC_oxy = oxy_Hb(:,48);
APFC_deoxy = deoxy_Hb(:,48);
APFC_total = total_Hb(:,48);

%scelgo canale 29 per DLPFC
DLPFC_oxy = oxy_Hb(:,29);
DLPFC_deoxy = deoxy_Hb(:,29);
DLPFC_total = total_Hb(:,29);

figure (1)
subplot(2,1,1)
plot(APFC_oxy,'r') 
hold on
plot(APFC_deoxy,'b')
title('APFC')
subplot(2,1,2)
plot(DLPFC_oxy,'r') 
hold on
plot(DLPFC_deoxy,'b')
title('DLPFC')
close all

%PUNTO b)----------------
%low pass filter
ft = 0.09;
[b,a]=butter(4,ft/fNy,'low');  
%freqz(b,a,100,fs)

%high pass filter
ft = 0.01;
[bb,aa]=butter(4,ft/fNy,'high');  
%freqz(bb,aa,100,fs)
 
APFC_oxy_filt = filtfilt(b,a,APFC_oxy);
APFC_deoxy_filt = filtfilt(b,a,APFC_deoxy);
APFC_oxy_filt = filtfilt(bb,aa,APFC_oxy_filt);
APFC_deoxy_filt = filtfilt(bb,aa,APFC_deoxy_filt); 

DLPFC_oxy_filt = filtfilt(b,a,DLPFC_oxy);
DLPFC_deoxy_filt = filtfilt(b,a,DLPFC_deoxy);
DLPFC_oxy_filt = filtfilt(bb,aa,DLPFC_oxy_filt);
DLPFC_deoxy_filt = filtfilt(bb,aa,DLPFC_deoxy_filt); 

figure (2)
subplot(2,1,1)
plot(APFC_oxy_filt,'r') 
hold on
plot(APFC_deoxy_filt,'b')
title('APFC filtered')
subplot(2,1,2)
plot(DLPFC_oxy_filt,'r') 
hold on
plot(DLPFC_deoxy_filt,'b')
title('DLPFC fitered')
close all

%PUNTO c)------------------
%divido in epoche di 50s

for i=1:6
MA = data{1,1}.trial(i*2-1);
APFC_oxy_epochs(:,i) = APFC_oxy_filt(MA-10*fs:(MA+40*fs)-1,1);
APFC_deoxy_epochs(:,i) = APFC_deoxy_filt(MA-10*fs:(MA+40*fs)-1,1);
DLPFC_oxy_epochs(:,i) = DLPFC_oxy_filt(MA-10*fs:(MA+40*fs)-1,1);
DLPFC_deoxy_epochs(:,i) = DLPFC_deoxy_filt(MA-10*fs:(MA+40*fs)-1,1);
end

if j==1
APFC_oxy_epochs_tot = APFC_oxy_epochs;
APFC_deoxy_epochs_tot = APFC_deoxy_epochs;
DLPFC_oxy_epochs_tot = DLPFC_oxy_epochs;
DLPFC_deoxy_epochs_tot = DLPFC_deoxy_epochs;
end

if j>1
APFC_oxy_epochs_tot = [APFC_oxy_epochs_tot APFC_oxy_epochs];
APFC_deoxy_epochs_tot = [APFC_deoxy_epochs_tot APFC_deoxy_epochs];
DLPFC_oxy_epochs_tot = [DLPFC_oxy_epochs_tot DLPFC_oxy_epochs];
DLPFC_deoxy_epochs_tot = [DLPFC_deoxy_epochs_tot DLPFC_deoxy_epochs];
end

end

%PUNTO d)-----media del subject S0%d--------
dim = size(APFC_oxy_epochs_tot);
n_epoch = dim(2);

APFC_oxy_mean = mean(APFC_oxy_epochs_tot');
APFC_oxy_stde = std(APFC_oxy_epochs_tot')/sqrt(n_epoch);
APFC_deoxy_mean = mean(APFC_deoxy_epochs_tot');
APFC_deoxy_stde = std(APFC_deoxy_epochs_tot')/sqrt(n_epoch);
DLPFC_oxy_mean = mean(DLPFC_oxy_epochs_tot');
DLPFC_oxy_stde = std(DLPFC_oxy_epochs_tot')/sqrt(n_epoch);
DLPFC_deoxy_mean = mean(DLPFC_deoxy_epochs_tot');
DLPFC_deoxy_stde = std(DLPFC_deoxy_epochs_tot')/sqrt(n_epoch);

%plot APFC
plot(APFC_oxy_mean,'r')
hold on 
plot(APFC_oxy_mean+APFC_oxy_stde,'--')
hold on
plot(APFC_deoxy_mean,'b')
hold on
plot(APFC_deoxy_mean+APFC_deoxy_stde,'--')
title('S01 - HHB and O2HB signal, APFC channel')
close all

%plot DLPFC
plot(DLPFC_oxy_mean,'r')
hold on
plot(DLPFC_oxy_mean+APFC_oxy_stde,'--')
hold on 
plot(DLPFC_deoxy_mean,'b')
hold on
plot(DLPFC_deoxy_mean+APFC_deoxy_stde,'--')
title('S01 - HHB and O2HB signal, DLPFC channel')
close all
%pause

%calcolo la somma dei vari pazienti
if k==1
    APFC_oxy_totS = APFC_oxy_epochs_tot;
    APFC_deoxy_totS = APFC_deoxy_epochs_tot;
    DLPFC_oxy_totS = DLPFC_oxy_epochs_tot;
    DLPFC_deoxy_totS = DLPFC_deoxy_epochs_tot;
end
if k>1
   APFC_oxy_totS = [APFC_oxy_totS APFC_oxy_epochs_tot];
   APFC_deoxy_totS = [APFC_deoxy_totS APFC_deoxy_epochs_tot];
   DLPFC_oxy_totS = [DLPFC_oxy_totS DLPFC_oxy_epochs_tot];
   DLPFC_deoxy_totS = [DLPFC_deoxy_totS DLPFC_deoxy_epochs_tot];
end

end

%PUNTO f)----------media di tutti i subject--------------
dim = size(APFC_oxy_totS);
n_epoch = dim(2);

APFC_oxy_mean_Stot = mean(APFC_oxy_totS');
APFC_oxy_stde_Stot = std(APFC_oxy_totS')/sqrt(n_epoch);
APFC_deoxy_mean_Stot = mean(APFC_deoxy_totS');
APFC_deoxy_stde_Stot = std(APFC_deoxy_totS')/sqrt(n_epoch);
DLPFC_oxy_mean_Stot = mean(DLPFC_oxy_totS');
DLPFC_oxy_stde_Stot = std(DLPFC_oxy_totS')/sqrt(n_epoch);
DLPFC_deoxy_mean_Stot = mean(DLPFC_deoxy_totS');
DLPFC_deoxy_stde_Stot = std(DLPFC_deoxy_totS')/sqrt(n_epoch);

%plot APFC
plot(APFC_oxy_mean_Stot,'r')
hold on 
plot(APFC_oxy_mean_Stot+APFC_oxy_stde_Stot,'--')
hold on
plot(APFC_oxy_mean_Stot-APFC_oxy_stde_Stot,'--')
hold on
plot(APFC_deoxy_mean_Stot,'b')
hold on
plot(APFC_deoxy_mean_Stot+APFC_deoxy_stde_Stot,'--')
hold on
plot(APFC_deoxy_mean_Stot-APFC_deoxy_stde_Stot,'--')
title('ALL SUBJECTS - HHB and O2HB signal, APFC channel')
close all
