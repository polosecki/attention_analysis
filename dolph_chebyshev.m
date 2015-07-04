function [x] = dolph_chebyshev(fs,fc,as)

asa = 10.^(-as/20);
thc = 2*pi*fc/fs;

m = acosh(1./asa) ./ (2*acosh(1./cos(thc/2)));

mr = ceil(m);

n = 2*mr+1;

x = zeros(1,n);

x0 = 1./cos(thc/2);
asa = 1 ./ cosh(2*mr*acosh(x0));

for i=0:mr
    sm = 0;
    thi = 2*pi*i/n;
    for j=1:mr
        thj = 2*pi*j/n;
        sm = sm + cos(2*mr*acos(x0*cos(thj/2))).*cos(j*thi);
    end
    x(i+mr+1) = 1/n*(1+2*asa*sm);
    x(mr-i+1) = 1/n*(1+2*asa*sm);
end

% f = 0:fs/16384:fs/2;
% fx = abs(fft(x,16384));
% fx = fx(1:8193);
% db = 20*log10(fx);
% subplot(2,1,1);
% plot(f,db);
% grid on;
% subplot(2,1,2);
% plot(x);
% grid on;
