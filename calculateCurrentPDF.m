function [dat]=calculateCurrentPDF(mat, E, T, S, dt)
%CALCULATECURRENTPDF this function calculates a current data set for given
%parameters. this method calculates the amount of charge measured, not the
%current and for that one needs to devide the returned matrix by 1e-6.

%mat - the material header. needs to be in the 'mats' file in the file './mats.mat'
%E - E field strength in V/cm
%T - temperature in kelvin
%S - The surface area of the emitter in cm^2
%dt - the time resolution of the measurement in sec (default 1e-6)

    if(~exist('dt','var'))
        dt=1e-6;
    end
    
    %load from mats file
    temp=load('./mats.mat');
    mats=temp.mats;
    mat=getfield(mats,mat); %#ok<GFLD>
    work=mat.work;
    Ef=mat.Ef;
    
    e=1.6e-19;%electron charge in coloumbs
    c=3e10;%light speed cm/sec
    h=4.135668e-15;%planck's constant in eV*sec;
    m=0.51e6/c^2;%electron mass in kg
    k=8.617e-5;%boltzmann's constant in eV/K

    fac=e*4*pi*m*k*T/h^3;%normalization factor
    
    %this piece of code takes the entire spectrum and returns the part that
    %has the highest contribution to 
    x=0:0.1:0.99*(work+Ef);
    n_temp=supply(x,T,Ef)*dt;
    D_temp=transparency(x,E,work,Ef);
    D_temp([false,(D_temp(2:end)-D_temp(1:end-1)<0)])=1;
    temp=n_temp.*D_temp.*fac;
    max_val=max(temp);
    ind=find(temp>max_val/10000);
    xstart=x(ind(1));
    xfin=x(ind(end));

    dx=0.001;%ev
    energy=xstart:dx:xfin;%relevant energy specrtum in eV
    n=supply(energy,T,Ef)*S;
    D=transparency(energy,E,work,Ef);
    D([false,(D(2:end)-D(1:end-1)<0)])=1;
    j=n.*D.*fac;
    mu=j*dx;%the mean current of each energy value 
    sigma_sq=e*(1-D).*mu*dt^(-1);%the variance of each energy value
    
    %this part generates a random sample of currents for each energy band
    %that can be later summed for the "measured" current
    N=1e4;
    dat=randn(N,length(mu)).*repmat(sqrt(sigma_sq),N,1)+repmat(mu,N,1);
    
end

