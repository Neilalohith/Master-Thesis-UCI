
function dpdw = pol(t,pw,delta,t2,t1,E0,f)

h = 6.62607015e-34 ;%* 1e-12;  %to ps, M1L1T-1
hbar=h/(2*pi);
%ft=linspace(0,1e-12,2000);
kappasq=(2*5.2*3.33564095e-30/hbar)^2;
efield=E0;
%efield = E0*exp(-1i*2*pi*f*ft) + E0*exp(1i*2*pi*f*ft);
dpdw = zeros(2,1);

%efield = interp1(ft,efield,t,'spline');

dpdw(1) = ((((1i*delta)-(1/t2))*pw(1)) - ((h/(8*pi))*1i.*efield*pw(2)*kappasq));
dpdw(2) = (((-pw(2))/t1) - ((8*pi/h)*imag(conj(pw(1)).*efield)));
end 