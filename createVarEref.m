%this code creates a var E reference data set. based on the varEsim2 code
%which is why i still use the rep_fac variable and the way the code is
%structured

function data=createVarEref(type, params, beta_range, E0_range)
%type - a string of wether to use Gaussian ('gauss'), sinusoidal ('sine') 
%   or saw-tooth ('saw') statistics
%params - a vector of the parameters of the distribution:
%   for Gaussian: the STD of E0 in V/m (constant in field)
%   for sinusoidal: magnitude of sine in V/m
%   for saw-tooth: magnitude of variation in V/m, in the range 
%   [-params(1),+params(1)).
%   FOR ALL - the last two parameters are the number of points used and
%   integration time in seconds (respectively).
%   NOTE - for the sine and saw cases the points are destributed in the
%   cycle so the the next point to have been measured is the first one
%   again.
%beta_range - the expected beta range in the system.
%E0_range - range of the applied fields to run in V/m

%data - the reference data. the data is also saved during the run of the
%   code.

    switch type
        case 'gauss'
            fun_type=1;
        case 'sine'
            fun_type=2;
        case 'saw'
            fun_type=3;
    end

    n_b=7;%number of beta values to use in the range
    rep_fac=10;%number of E0 samples
    N=5*rep_fac;%5 repeating measurements for each unique E0 values
    %sigmaE is the 'magnitude' of variation - in the gaussian case it's the
    %std and in the sine and saw cases it's the amplitude of the
    %wave-function (from -sigmaE to sigmaE)
    sigmaE=ones(n_b,N)*params(1)/100;%variation amplitude in V/cm
    E0_range=E0_range/100;%range in V/cm
    sigE=zeros(n_b,N);%std of simulated field
    sigI=zeros(n_b,N);%std of simulated current
    E_used=zeros(n_b,N);%the E values used after the gaussian random
    beta=zeros(n_b,N);%This is also a remnant of the varEsim2 code, the actual range is 1-3*1e8V/cm (100-300MV/cm)
    meanI=zeros(n_b,N);%mean I value
    meanE=zeros(n_b,N);%mean E value
    beta_r=linspace(beta_range(1),beta_range(2),n_b);%the range of beta values to use
    prog=0;
    wb=waitbar(prog,['E var refrence generation at ',num2str(prog*100),'%']);

    T=293;
    S=(1)^2;%tunneling surface in cm^2. This value does not affect results.
    
    m=params(end-1);%number of points to use for estimation
    dt=params(end);%integration time
    
    tic;
    for k=1:length(beta_r)
        for i=1:N
            beta(k,i)=beta_r(k);%beta
            E_used(k,i)=E0_range(1)+floor(i/(N+1)*rep_fac)*diff(E0_range)/(rep_fac-1);
            %applied field in V/cm
            E=(double(fun_type==1)*randn(1,m)+...%gaussian case
                double(fun_type==2)*(sin(rand*2*pi+[0:m-1]/m*2*pi))+...%sinusoidal case
                double(fun_type==3)*(rand*1/(m/2)+[0:m-1]/(m/2)-1))*...%saw-tooth case
                sigmaE(k,i)+E_used(k,i);%the amplitude and bias
            %additional cases can be added here to in the form of:
            %'+double(fun_type==n)*<expression>' with 'n' the number
            %representing the new case and '<expression>' being an
            %expression that produces a (1,m) sized vector.
            
            %this part simulates the current. it's generally the part that takes time
            %so it's the one with the waitbar
            I=zeros(1,m);
            for j=1:m
                temp=calculateCurrentPDF('Cu', E(j)*beta(k,i), T, S,dt);%returns 2.5e4 points
                temp_I=sum(temp,2);
                
                I(j)=mean(temp_I);
                
                prog=(j+m*(i-1)+(k-1)*m*N)/(m*N*n_b);
                t=toc;
                waitbar(prog,wb,['E var refrence generation at ',num2str(prog*100,'%0.2f')...
                    ,'% estimated time: ',num2str(t/prog/60-t/60,'%0.2f'),' minutes']);
                clear temp temp_I
            end
            
            
            sigI(k,i)=std(I);
            sigE(k,i)=std(E);
            meanI(k,i)=mean(I);
            meanE(k,i)=mean(E);
            
        end
    end
    
    close(wb)
    
    data.type=type;
    data.m=m;
    data.dt=dt;
    data.E0_range=E0_range*100;
    data.beta_range=beta_range;
    data.E0_range=E0_range*100;
    data.meanE=meanE;
    data.meanI=meanI;
    data.sigE=sigE;
    data.sigI=sigI;
    data.beta=beta;
    data.sigmaE=sigmaE;
    data.beta_range=beta_range;
    
    E0=data.meanE(:)/1e6;
    b=data.beta(:);
    y=data.sigI(:)./(data.sigE(:).*data.meanI(:));
    
    [f,g]= fit([E0,b],y,'poly53');
    data.f=f;
    data.g=g;
    
end
