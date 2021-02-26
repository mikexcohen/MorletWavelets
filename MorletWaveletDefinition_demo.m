%%%
%    Morelet wavelet FWHM definition
%    
%    This code accompanies the paper:
%           "A better way to define and describe Morlet wavelets for time-frequency analysis"
%           by MX Cohen (NeuroImage, 2019)
% 
%    Questions? -> mikexcohen@gmail.com
%%%

%% Figure 1: The happy time-frequency man

clear

%%% generate signal
srate = 1000;
time  = (0:3*srate)/srate;
pnts  = length(time);

f1 = [3 20];
f2 = [20 3];

transgaus = exp( -(time-mean(time)).^2 / .2 );

signal = sin(2*pi*time.*linspace(f1(1),mean(f1),pnts)) + ...
         sin(2*pi*time.*linspace(f2(1),mean(f2),pnts)) + ...
         sin(2*pi*time*20).*transgaus;


%%% static power spectrum
hz = linspace(0,srate/2,floor(pnts/2)+1);
powr = (2*abs(fft(signal)/pnts)).^2;


%%% time-frequency analysis
nfrex = 40;
frex = linspace(2,25,nfrex);
fwhm = linspace(.8,.7,nfrex);

% setup wavelet and convolution parameters
wavet = -2:1/srate:2;
halfw = floor(length(wavet)/2)+1;
nConv = pnts + length(wavet) - 1;

% initialize time-frequency matrix
tf = zeros(length(frex),pnts);

% spectrum of data
dataX = fft(signal,nConv);

% loop over frequencies
for fi=1:length(frex)
    
    % create wavelet
    waveX = fft( exp(2*1i*pi*frex(fi)*wavet).*exp(-4*log(2)*wavet.^2/fwhm(fi).^2),nConv );
    waveX = waveX./max(waveX); % normalize
    
    % convolve
    as = ifft( waveX.*dataX );
    % trim and reshape
    as = as(halfw:end-halfw+1);
    
    % power
    tf(fi,:) = abs(as).^2;
end


%%% plotting

figure(1), clf
subplot(211)
plot(time,signal,'k','linew',2)
xlabel('Time (s)'), ylabel('Amplitude')
set(gca,'xlim',time([1 end]))
title('A) Time domain')


subplot(223)
plot(hz,powr(1:length(hz)),'k','linew',2)
set(gca,'xlim',[0 30])
xlabel('Frequency (Hz)'), ylabel('Power')
title('B) Frequency domain')


subplot(224)
contourf(time,frex,tf,40,'linecolor','none')
xlabel('Time (sec.)'), ylabel('Frequency (Hz)')
title('C) Time-frequency domain')

colormap(bluewhitered(64))


%% Figure 2: Big confusion from little wavelets

clear

% simulation parameters
freq1 = 7; % of the wavelet
freq2 = 12;
srate = 1000;
time  = -2:1/srate:2;
pnts  = length(time);

% define number of cycles
numcycles = [ 3 8 ];

% create the sine wave and two Gaussian windows
sinwave1 = cos(2*pi*freq1*time);
sinwave2 = cos(2*pi*freq2*time);
gauswin1 = exp( -time.^2 / (2* (numcycles(1)/(2*pi*freq1))^2 ) );
gauswin2 = exp( -time.^2 / (2* (numcycles(2)/(2*pi*freq1))^2 ) );
gauswin3 = exp( -time.^2 / (2* (numcycles(2)/(2*pi*freq2))^2 ) );


% create the three wavelets
morletwave1 = sinwave1 .* gauswin1;
morletwave2 = sinwave1 .* gauswin2;
morletwave3 = sinwave2 .* gauswin3;


% normalized power spectra of the wavelets
powspect1 = abs(fft(morletwave1)/pnts).^2;
powspect1 = powspect1 ./ max(powspect1);

powspect2 = abs(fft(morletwave2)/pnts).^2;
powspect2 = powspect2 ./ max(powspect2);

powspect3 = abs(fft(morletwave3)/pnts).^2;
powspect3 = powspect3 ./ max(powspect3);


hz = linspace(0,srate/2,floor(pnts/2)+1);

% compute the empirical temporal FWHM in seconds (later converted to ms)
midp = dsearchn(time',0);
fwhmT(1) = time(midp-1+dsearchn(gauswin1(midp:end)',.5)) - time(dsearchn(gauswin1(1:midp)',.5));
fwhmT(2) = time(midp-1+dsearchn(gauswin2(midp:end)',.5)) - time(dsearchn(gauswin2(1:midp)',.5));
fwhmT(3) = time(midp-1+dsearchn(gauswin3(midp:end)',.5)) - time(dsearchn(gauswin3(1:midp)',.5));


idx1 = dsearchn(hz',freq1);
idx2 = dsearchn(hz',freq2);
fwhmF(1) = hz(idx1-1+dsearchn(powspect1(idx1:end)',.5)) - hz(dsearchn(powspect1(1:idx1)',.5));
fwhmF(2) = hz(idx1-1+dsearchn(powspect2(idx1:end)',.5)) - hz(dsearchn(powspect2(1:idx1)',.5));
fwhmF(3) = hz(idx2-1+dsearchn(powspect3(idx2:end)',.5)) - hz(dsearchn(powspect3(1:idx2)',.5));






%%% plotting
figure(2), clf
subplot(231)
plot(time,morletwave1,'k')
xlabel('Time (s)'), axis square
title([ num2str(numcycles(1)) ' cycles, FWHM: ' num2str(fwhmT(1)*1000) ' ms' ])
set(gca,'xlim',[-1 1])

subplot(232)
plot(time,morletwave2,'k')
xlabel('Time (s)'), axis square
title([ num2str(numcycles(2)) ' cycles, FWHM: ' num2str(fwhmT(2)*1000) ' ms' ])
set(gca,'xlim',[-1 1])

subplot(233)
plot(time,morletwave3,'k')
xlabel('Time (s)'), axis square
title([ num2str(numcycles(2)) ' cycles, FWHM: ' num2str(fwhmT(3)*1000) ' ms' ])
set(gca,'xlim',[-1 1])




subplot(234)
plot(hz,powspect1(1:length(hz)),'k','linew',2)
set(gca,'xlim',[0 freq2*2])
xlabel('Frequency (Hz)'), axis square
title([ num2str(numcycles(1)) ' cycles, FWHM: ' num2str(fwhmF(1)) ' Hz' ])

subplot(235)
plot(hz,powspect2(1:length(hz)),'k','linew',2)
set(gca,'xlim',[0 freq2*2])
xlabel('Frequency (Hz)'), axis square
title([ num2str(numcycles(2)) ' cycles, FWHM: ' num2str(fwhmF(2)) ' Hz' ])

subplot(236)
plot(hz,powspect3(1:length(hz)),'k','linew',2)
set(gca,'xlim',[0 freq2*2])
xlabel('Frequency (Hz)'), axis square
title([ num2str(numcycles(2)) ' cycles, FWHM: ' num2str(fwhmF(3)) ' Hz' ])


%% Figure 3: FWHM vs. number-of-cycles

clear


% specify frequencies
frex = linspace(3,60,50);
fwhm = linspace(.4,.1,length(frex)); % in ms
fwhm = logspace(log10(.4),log10(.1),length(frex)); % in ms
ncyc = linspace(3.2,16,length(frex));



% parameters for complex Morlet wavelets
srate = 1024;
wavtime  = -2:1/srate:2;
midp = dsearchn(wavtime',0);


% outputs
empfwhm = zeros(length(frex),2);

% loop over frequencies
for fi=1:length(frex)
    
    % create the Gaussian using the FWHM formula (equation 3)
    gwin = exp( (-4*log(2)*wavtime.^2) ./ fwhm(fi)^2 );
    
    % measure the empirical fwhm
    empfwhm(fi,1) = wavtime(midp-1+dsearchn(gwin(midp:end)',.5)) - wavtime(dsearchn(gwin(1:midp)',.5));
    
    
    
    % create the Gaussian using the n-cycles formula (equations 1-2)
    s    = ncyc(fi) / (2*pi*frex(fi));
    gwin = exp( -wavtime.^2 ./ (2*s^2) );
    
    % empirical FWHM
    empfwhm(fi,2) = wavtime(midp-1+dsearchn(gwin(midp:end)',.5)) - wavtime(dsearchn(gwin(1:midp)',.5));
end

figure(3), clf

plot(frex,empfwhm*1000,'o-','markersize',8,'markerfacecolor','w','linew',2)
xlabel('Wavelet frequency (Hz)'), ylabel('FWHM (ms)')
legend({'Using FWHM';'Using n-cycles'})
set(gca,'xlim',[frex(1)-1 frex(end)+1])

%% Figure 4: Defining wavelets in the frequency domain

clear


% specify wavelet parameters
peakf = 11;
fwhm  = 5.2;

% specify simulation details
npnts = 8001;
srate = 1000;


% vector of frequencies
hz = linspace(0,srate,npnts);

% frequency-domain Gaussian
s  = fwhm*(2*pi-1)/(4*pi); % normalized width
x  = hz-peakf;             % shifted frequencies
fx = exp(-.5*(x/s).^2);    % gaussian

% empirical FWHM
idx = dsearchn(hz',peakf);
empfwhmF = hz(idx-1+dsearchn(fx(idx:end)',.5)) - hz(dsearchn(fx(1:idx)',.5));



% time-domain wavelet
morletwavelet = fftshift( ifft(fx) );
time = (-floor(npnts/2):floor(npnts/2))/srate;

% FWHM of wavelet in the time domain
midp = dsearchn(time',0);
mw_amp = abs(morletwavelet);
mw_amp = mw_amp./max(mw_amp);
empfwhmT = time(midp-1+dsearchn(mw_amp(midp:end)',.5)) - time(dsearchn(mw_amp(1:midp)',.5));



%%% plotting
figure(5), clf
subplot(211)
plot(hz,fx,'k','linew',2)
set(gca,'xlim',[0 peakf*3])
xlabel('Frequency (Hz)'), ylabel('Amplitude (gain)')
title([ 'FWHM specified: ' num2str(fwhm) ', obtained: ' num2str(empfwhmF) ' Hz' ])

subplot(212), hold on
plot(time,real(morletwavelet),'linew',2)
plot(time,imag(morletwavelet),'--','linew',2)
h = plot(time,abs(morletwavelet),'linew',2);
set(h,'color',[1 1 1]*.8)
set(gca,'xlim',[-1 1])
legend({'Real part';'Imag part';'Envelope'})
xlabel('Time (sec.)')
title([ 'FWHM: ' num2str(empfwhmT*1000) 'ms' ])

%% Figure 5a: The real deal (fixed FWHM)

% note: this is one channel from one dataset from 
% Cohen, M.X., 2015. Comparison of different spatial transformations applied to EEG data: A case study of error processing. Int. J. Psychophysiol. 97, 245–257

load MorletWaveletDefinition_data.mat

%%% time-frequency parameters
nfrex = 80;
frex  = linspace(2,40,nfrex);
fwhm  = [.1 .5 2];

% timimg parameters
bidx = dsearchn(EEG.times',[-500 -200]');
tidx = dsearchn(EEG.times',-500):dsearchn(EEG.times',1300);


% setup wavelet and convolution parameters
wavet = -5:1/EEG.srate:5;
halfw = floor(length(wavet)/2)+1;
nConv = EEG.pnts*EEG.trials + length(wavet) - 1;


% initialize time-frequency matrix
tf = zeros(length(frex),length(tidx),length(fwhm)+1);
empfwhm = zeros(length(fwhm),length(frex)+1);


% spectrum of data
dataX = fft(reshape(EEG.data,1,[]),nConv);


figure(5), clf
clim = [-2 2];


%%% loop over FWHM parameter settings
for fwhmi=1:length(fwhm)
    
    % loop over frequencies
    for fi=1:length(frex)
        
        % create wavelet
        waveX = fft( exp(2*1i*pi*frex(fi)*wavet).*exp(-4*log(2)*wavet.^2/fwhm(fwhmi).^2),nConv );
        waveX = waveX./max(waveX); % normalize
        
        % convolve
        as = ifft( waveX.*dataX );
        % trim and reshape
        as = reshape(as(halfw:end-halfw+1),EEG.pnts,EEG.trials);
        
        % power
        p = mean( abs(as).^2 ,2);
        tf(fi,:,fwhmi) = 10*log10( p(tidx)/mean(p(bidx(1):bidx(2))) );
        
        % empirical FWHM
        hz = linspace(0,EEG.srate,nConv);
        idx = dsearchn(hz',frex(fi));
        fx = abs(waveX);
        empfwhm(fwhmi,fi) = hz(idx-1+dsearchn(fx(idx:end)',.5)) - hz(dsearchn(fx(1:idx)',.5));

    end
    
    % plots.
    subplot(2,2,fwhmi)
    contourf(EEG.times(tidx),frex,squeeze(tf(:,:,fwhmi)),60,'linecolor','none')
    set(gca,'clim',clim,'YScale','log','YTick',logspace(log10(frex(1)),log10(frex(end)),10),'YTickLabel',round(logspace(log10(frex(1)),log10(frex(end)),10),2))
    axis square
    title([ num2str(fwhm(fwhmi)) ' sec' ])
end

colormap(bluewhitered(64))

%% Figure 5a: The real deal (variable FWHM)

fwhm = linspace(1,.2,length(frex));

% loop over frequencies
for fi=1:length(frex)
    
    % create wavelet
    waveX = fft( exp(2*1i*pi*frex(fi)*wavet).*exp(-4*log(2)*wavet.^2/fwhm(fi).^2),nConv );
    waveX = waveX./max(waveX); % normalize
    
    % convolve
    as = ifft( waveX.*dataX );
    % trim and reshape2 to 40 Hz in 
    as = reshape(as(halfw:end-halfw+1),EEG.pnts,EEG.trials);
    
    % power
    p = mean( abs(as).^2 ,2);
    tf(fi,:,4) = 10*log10( p(tidx)/mean(p(bidx(1):bidx(2))) );
    
    
    % empirical FWHM
    hz = linspace(0,EEG.srate,nConv);
    idx = dsearchn(hz',frex(fi));
    fx = abs(waveX);
    empfwhm(4,fi) = hz(idx-1+dsearchn(fx(idx:end)',.5)) - hz(dsearchn(fx(1:idx)',.5));

end

% plots.
subplot(2,2,4)
contourf(EEG.times(tidx),frex,squeeze(tf(:,:,4)),60,'linecolor','none')
set(gca,'clim',clim,'YScale','log','YTick',logspace(log10(frex(1)),log10(frex(end)),10),'YTickLabel',round(logspace(log10(frex(1)),log10(frex(end)),10),2))
axis square
title([ num2str(fwhm(1)) '-' num2str(fwhm(end)) ])

%% Figure 5b: the real deal, part deux

figure(6), clf

empfwhm(empfwhm>100) = nan;
plot(frex,empfwhm,'s','markerfacecolor','w','linew',2,'markersize',8)
legend({'100','500','2000','var'})
set(gca,'ylim',[0 18],'xlim',[frex(1)-1 frex(end)+1])
xlabel('Peak frequency (Hz)'), ylabel('FWHM (Hz)')
axis square

%% done.
