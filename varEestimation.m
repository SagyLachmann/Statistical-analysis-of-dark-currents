function [beta_est,beta_err]=varEestimation(sigI,sigE,E0,I0,l,ref_data)
%this function takes in a refrence data matrix (saved by the function
%'createVarEref') and the measured parameters to return an estimation of
%beta and the error in beta.
%the fit is assumed to be a power of 2 in beta so that two solutions exist
%for any pair of eta and E0. The power of E0 in the fit to eta isn't
%assumed to be any specific value.

%sigI - the std measured in the current (in A)
%sigE - the std measured in the applied field in V/m (std(V)/gap)
%E0 - the mean applied field in V/m
%I0 - the mean measured current in A
%l - the number of points used for the integration of the current
%   measurement. i.e. the integration time for a single current measurement
%   divided by the time resolution of the measurement. used for error
%   estimation.
%ref_data - a structure variable of the refrence data, returned by the
%   function 'createVarEref'. this can be either for Gaussian, sinusoidal
%   or saw-tooth statistics.
%NOTE - both values of variation (sigI and sigE) are the STD of their
%   respective variable, regardless of the statistics. The code is written
%   so that they all use the STD. It might be that for some statistics
%   other ways of identifying the variation (e.g. for sine to use the mean
%   or the initial phase but that was not explored in this work.

%beta_est - the estimated value of beta
%beta_err - the estimated error of the estimated beta

    %transfering to the proper units
    sigE=sigE/100;
    E0=abs(E0)/100/1e6;%the ref data mean field is normalized to be MV/cm in the fit. abs to account for negative bias
    
    f=ref_data.f;
    coef_names=coeffnames(f);
    coef_values=coeffvalues(f);
    
    eta_ex=abs(sigI/(sigE*I0));%experimentally measured value. abs for negative bias cases.
    
    syms b x;
    eq=0;%building the equation of the fit, assumed polynomial and that the "y" parameter is beta
    for i=1:length(coef_names)
        eq=eq+coef_values(i)*E0^str2double(coef_names{i}(2))*b^str2double(coef_names{i}(3));
    end
    
    beta_est=double(vpasolve(eq==eta_ex));%this returns beta
    %this beta estimation has two values, one of which should be outside
    %the beta range that was used for the reference data. this next part
    %chooses only one value and returns warnings if both are outside the
    %range (in which case the one closest is chosen).
    beta_est=reshape(beta_est,[],1);
    beta_est(imag(beta_est)~=0)=[];%remove imaginary
    if(length(beta_est)>1)
        inside=(beta_est>=ref_data.beta_range(1)).*(beta_est<=ref_data.beta_range(2));
        if(sum(inside)==0)%no one inside
            dis=min(abs(beta_est-ones(2,1)*ref_data.beta_range(1)),...
                abs(beta_est-ones(2,1)*ref_data.beta_range(2)));
            beta_est=beta_est(dis==min(dis));%return the value closest to the range 
            warning('both beta solutions are outside the desired range. returning the solution closest to the range.');
            %TODO throw warning
        elseif(sum(inside)==1)
            beta_est=beta_est(logical(inside));%return the value inside the range
        else%both are inside
            warning('both beta solutions are inside the desired range. returning both solutions.');
            %TODO add case
        end
    end
    
    %error estimation
    err=estimateErr(l,sigI,sigE,I0);%estimated error in eta_ex
    if(length(beta_est)==1)
        beta_p=double(vpasolve(eq==eta_ex+err,b,beta_est));
        beta_m=double(vpasolve(eq==eta_ex-err,b,beta_est));
        beta_p=beta_p(min(abs(beta_p-beta_est))==abs(beta_p-beta_est));
        beta_m=beta_m(min(abs(beta_m-beta_est))==abs(beta_m-beta_est));
        beta_err=abs(beta_p-beta_m)/2;
        %just as there are 2 beta solutions, there are two solutions for
        %the beta error.
    else
        beta_p1=double(vpasolve(eq==eta_ex+err,b,beta_est(1)));
        beta_m1=double(vpasolve(eq==eta_ex-err,b,beta_est(1)));
        beta_p1=beta_p1(min(abs(beta_p1-beta_est))==abs(beta_p1-beta_est));
        beta_m1=beta_m1(min(abs(beta_m1-beta_est))==abs(beta_m1-beta_est));
        beta_p2=double(vpasolve(eq==eta_ex+err,b,beta_est(2)));
        beta_m2=double(vpasolve(eq==eta_ex-err,b,beta_est(2)));
        beta_p2=beta_p2(min(abs(beta_p2-beta_est))==abs(beta_p2-beta_est));
        beta_m2=beta_m2(min(abs(beta_m2-beta_est))==abs(beta_m2-beta_est));
        beta_err=[abs(beta_p1-beta_m1)/2;abs(beta_p2-beta_m2)/2];
    end
end

%this function estimates the error of sigI/(sigE*I0) based on the supplied parameters
%and a pre-calculated matrix of errors for Gaussian parameter estimation
%error at '../data/mean and std relative error.mat'
%NOTE - the assumed gaussian distribution isn't that of the E0 wavefunction
%but of the signals measured in the scope for the "instanteneous" I and E0
%values that are assumed to have Gaussian noise. these errors dI0, dI, dE
%are the errors in evaluating the mean and std of a Gaussian variable given
%the number of points used in the evaluation.
function err=estimateErr(l,sigI,sigE,I0)

    load('./mean and std relative error.mat'); %#ok<LOAD>
    x=log10(l);
    dI0=10^(fm.p1*x+fm.p2)*I0;
    dI=10^(fs.p1*x+fs.p2)*sigI;
    dE=10^(fs.p1*x+fs.p2)*sigE;
    err=sqrt((1/sigE/I0)^2*dI^2+(sigI/sigE^2/I0)^2*dE^2+(sigI/sigE/I0^2)^2*dI0^2);

end
