global KD KD2 HFold Hunfold

KD = 2000; % diffusion constant. It is given by D/(a^2), where D - diffusion coefficient, a lengh of DNA segment. 
KD2 = 15000; % second DNA diffusion constant. *60 is to convert this time to /minutes for longer calculation times
HFold = 15; % folding/unfolding constant
Hunfold = 15; 

N = 30;
Niter = 340000*5;
Lsave = zeros(Niter,N);

V = zeros(N,1);

figure
hold on
xlabel('Time (s)')
ylabel('Loop size (bp)')
for i=1:N
    L=Loop5(200);
    [L,TS,NS,HF,HB]=L.MMK(Niter);
    plot(TS,NS*5/0.34,'.-')
    p = polyfit(TS,NS*5/0.34,1);
    V(i) = p(1);
    pause(0.5)
    Lsave(:,i) = NS*5/0.34;         % in nanometers
    title(['i = ' num2str(i) ' of ' num2str(N)])
end
h1 = axis;
yyaxis right
axis([h1(1) h1(2) h1(3)*0.34/1000 h1(4)*0.34/1000])
ylabel('Loop size (um)')
yyaxis left
title(['LE rate: ' num2str(round(mean(V))) ' +/- ' num2str(round(std(V))) ' bps'])


figure
histogram(V,10)

disp(['LE rate: ' num2str(round(mean(V))) ' +/- ' num2str(round(std(V))) ' bps'])

z= dotplot(V);
figure
plot(z,V,'o')

