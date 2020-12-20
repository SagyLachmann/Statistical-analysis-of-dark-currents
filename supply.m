%the supply function as calculated by Nordheim.
%takes the energy range of the electron (Ex) the temperature of the system
%(T) and the Fermi level of the metal to calculate the supply function (n)
%both Ex and Ef are in eV
function [n] = supply(Ex,T,Ef)
    e=1;%electron charge in elementary units
    h=4.136e-15;%planck's constant in eV*sec;
    m=9.109e-31;%electron mass in kg
    k=8.617e-5;%boltzmann's constant in eV/K
%     n=(4*pi*m*k*T)/(h^3)*log(1+exp((Ex-Ef)/(k*T)));%the supply function
    n=log(1+exp(-(Ex-Ef)/(k*T)));
%this modification was inserted for low temperatures where matlab thinks 
%n=inf for low energy - it calculates the exp which gets to infinity but
%the log of it doesn't go to infinity.
    i=(n==inf);
    n(i)=-(Ex(i)-Ef)/(k*T);
end
%the supply function is the flux of electrons hitting the barrier