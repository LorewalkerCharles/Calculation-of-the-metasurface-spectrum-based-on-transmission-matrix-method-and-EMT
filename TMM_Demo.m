load nmms.mat
load nsi.mat
Freq = 0.1e12:2.4e9:2.5e12;
[r,t] = TransMatMthd(Freq,nmms,nsi,43e-6,244e-6);
figure(1)
plot(Freq',abs(r).^2,'-.')
hold on
plot(Freq',abs(t).^2)
plot(Freq',1-abs(r).^2-abs(t).^2,'--')
axis([1e12 1.75e12 0 1.05])
legend('r','t','absor')

