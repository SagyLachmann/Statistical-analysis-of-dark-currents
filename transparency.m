%the quantum transparency function as calculated by Nordheim.
%takes the energy range of the electron (Ex) the applied electric field (E)
%and the work function of the metal (work) to return the transparency of
%the barrier (D)
%both Ex and work are in eV and E0 is in V/m
function [D] = transparency(Ex,E0,work,Ef)
%     e=1;%electron charge in elementary units
%     h=4.136e-15;%planck's constant in eV*sec;
%     m=9.109e-31;%electron mass in kg
%     y=sqrt(e^3*E0)/work;%parameter for Nordheim's function
%     theta=1-0.142*y-0.855*y^2;
%     convert=1.602e-19;%the convertion from eV to joules
%     Ex=Ex*convert;
%     work=work*convert;
%     h=h*convert;
%     kappa=4*pi*sqrt(m)/h;
%     convert=1;%for some reason this convertion is unnecessary. i will look at the original paper to see why.
%     D=exp((8*pi*sqrt(2*m))/(3*e*h)*(Ex.^1.5/E0)/convert*theta);
%     D=exp(-4*kappa*(work-Ex).^1.5/(3*E0)).*(4*(Ex.*(work-Ex)).^0.5/work);
    
    %Sommerfeld and bethe
    y=3.79e-4*E0^0.5./(work+Ef-Ex);
    v=ellipticFNB(y);%eliptic function by nordheim Proc. Roy. Soc. A121, 626 (1928)
        %corrected by burgess kromer and huston phys rev 90 515 (1953)
    %for values of y greater than 1, f has an imagionary part but for our
    %purposes it's smaller than 1e-15 so i just take the real part.
    
    %used to be without the 1, corrected based on murphy-good and the error was for high fields
    D=(1+exp(6.83e7*(work+Ef-Ex).^1.5.*v/E0)).^(-1);
    D(y<0)=1;
end
%the transparency is the chance of an electron with energy Ex to pass
%through the barrier
%the convertion from eV to joules is for convenience in choosing the units
%of the mass, electric field and energy to use the charge as unity.

%this f function (denotion by dyke&dolan) is the FN v function (as denoted 
%by literally anyone else)
function v=ellipticFNB(y) 
%     a=(1-y.^2).^0.5;
%     ksq=(2*a).^0.5./(1+a).^0.5;
%     f=2^(-0.5)*(1+a).^0.5.*(ellipticE(ksq)-(1-a).*ellipticK(ksq));
%     ksq=2*(1-y.^2).^0.5./(1+(1-y.^2).^0.5);
%     f=2^(-0.5)*(1+(1-y.^2).^0.5).^0.5.*(ellipticE(ksq)-y.^2.*ellipticK(ksq)./(1+(1-y.^2).^0.5));
%     f=(f+conj(f))/2;
    v=1-y.^2+1/3*y.^2.*log(y);
    v(y==0)=1;%this is what the function approaches as y goes to 0 though the actual calculated value is NaN
end