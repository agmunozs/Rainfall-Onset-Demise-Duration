% Matlab script to compute onset/demise (and duration)
% AG Munoz (IRI Columbia U) - agmunoz@iri.columbia.edu
% Project: Predictability of the onset, duration and demise of the North and South American Monsoon Systems: the role of cross-equatorial interactions
% First edition: Aug 4, 2017 
% Last edition: Jan 30, 2018
% Notes: 
% + Subset of the full code written for the N/S American Monsoon
%   project, to show the Bombardi and Carvalho approach. 
% + For illustration purposes, only one gridbox/station is used.
% + The method has issues when the year is too dry (as one should expect).


%%%%%START OF USER-MODIFIABLE SECTION%%%%%%%%%%%%

disp('Start...');
% set working directory
clear all
% set working directory
cd /Users/agmunoz/Documents/Angel/Weather_within_climate
%parpool('local')

%Read data via OpenDAP?
down=1;   %1=yes; 0=no (0 assumes the data is locally available in *.mat files; *not* NetCDF!)

%Define temporal parameters:
seasons='Jan';
seasone='Dec';
yeari=1981; %first year (MUST BE >=1995 AND <2005!)
yeare=2014; %last year  (MUST BE >2006!)    %Note: for chi there're missing values in 2009!!!!

%Define gridbox coordinates:
slat=18.3810;
slon=-67.1570;

%%%END OF USER-MODIFIABLE SECTION (DO NOT MODIFY ANYTHING BELOW THIS LINE)%%%%%
%%
%Onset/Demise/Duration
%Pentads:
if down==1  
 clear pentad
 %The following URL uses the entire Jan-Dec period!!!! 
 iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.RETRO/.rain/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.REALTIME/.rain/appendstream/T/(0000%201%20Jan%20' num2str(yeari) ')/(0000%201%20Jan%202006)/RANGEEDGES/pentadAverage/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.REALTIME/.rain/T/(0000%201%20Jan%202006)/(0000%201%20Jan%20' num2str(yeare+1) ')/RANGEEDGES/pentadAverage/appendstream/T/73/splitstreamgrid/X/' num2str(slon) '/VALUE/Y/' num2str(slat) '/VALUE/dods'];
 pentad = squeeze(double(ncread(iridl,'rain')));  %pentad 1 is 1-5 Jan 1981. There're 73 pentads per "year", and a total of 34 "years" (until 2014)
 save -v7.3 pentad.mat pentad 
else
 load pentad.mat pentad 
end
 clear penmean onset demise dur
 penmean=squeeze(nanmean(pentad,2));
 ns=size(pentad);nyear=ns(2);
%%
% %Computing onset/demise a la Bombardi and Carvalho 2009:
clear onset demise inipent inipent2


if down==1
%Find initial pentad for onset
display('Computing initial pentad... (this takes some time)')

%inipent/inipent2 is t0 in Bombardi and Carvalho 2009, initial pentad for
%onset and demise, respectively. Here wrt Jan 1 every year

        clear f1 f2
            if(~(isnan(penmean(:))))
                %Compute first harmonic of the climatological rainfall for
                %the location:
                f1 = fit([1:73]',squeeze(penmean(1:73)),'fourier2');
                %inipent=round(fminbnd(f1,1,73)); %min, for onset
                inipent=round(interp1(feval(f1,[1:73]'),1:73,min(feval(f1,1:73))));%min, for onset
                %f2 = fit([1:73]',squeeze(-penmean(1:73)),'fourier2');  %now negative for maximum
                %inipent2=round(fminbnd(f2,1,73));%max, for demise
                inipent2=round(interp1(feval(f1,[1:73]'),1:73,max(feval(f1,1:73))));%max, for demise
            else
                inipent=NaN;
                inipent2=NaN;
            end
            
inipent(inipent==0)=NaN;
inipent2(inipent2==0)=NaN;
display('t0 Done!')
figure(1);clf
plot(f1);hold on; plot(1:73,mean(pentad,2),inipent,f1(inipent),'*k',inipent2,f1(inipent2),'*k','markers',12)

%Compute S
clear S S0
display('Computing S...')

        for ity=1:nyear-1 
            S0=0;
            S02=0;
            %if(isnan(penmean) | isnan(inipent) | isnan(inipent2))
            %    S(ity)=NaN;
            %    S2(ity)=NaN;   %S2 for demise, using inipent2 
            %else
                for itp=1:73
                    itk=inipent+itp;
                    itk2=inipent2+itp;  %for demise, using inipent2 
                    if itp == 1
                        S(ity,itp)=squeeze(pentad(inipent,ity)-penmean(itp,1));
                        S2(ity,itp)=pentad(inipent2,ity)-penmean(itp,1);
                    else
                        if itk<=73
                            S(ity,itp)=sum(pentad(inipent:itk,ity)-penmean(itp,1),1);
                        else %use next year's data  
                            S0=sum(pentad(inipent:73,ity)-penmean(itp,1),1); %accumulated until the end of previous year
                            S(ity,itp)= S0 + ( sum(pentad(1:itk-73,ity+1)-penmean(itp,1),1) ); %total, considering demise on following year
                         end
                        if itk2<=73  %demise
                            S2(ity,itp)=sum(pentad(inipent2:itk2,ity)-penmean(itp,1),1);
                        else %use next year's data  
                            S02=sum(pentad(inipent2:73,ity)-penmean(itp,1),1); %accumulated until the end of previous year
                            S2(ity,itp)= S0 + ( sum(pentad(1:itk2-73,ity+1)-penmean(itp,1),1) ); %total, considering demise on following year
                        end
                    end
                end
            %end
        end

S(S==0)=NaN;
S2(S2==0)=NaN;
display('S Done!')

%Compute Sb
clear Sb
display('Smoothing with a Savitzky-Golay 2nd-degree filter...')

        for ity=1:nyear-1
            if( isnan(S(ity,:)) | isnan(S2(ity,:)) )
                Sb(ity,:)=NaN;
                Sb2(ity,:)=NaN;
            else
                %for pentads, R. Bombardi suggests to smooth S ~20 times;
                %for daily data ~150 times !!! NOT USING THAT, USING SGOLAY
                 Sb(ity,1:73)=smooth(S(ity,1:73),30,'sgolay');
                 Sb2(ity,1:73)=smooth(S2(ity,1:73),30,'sgolay');
            end
        end

display('Sb Done!')

%Compute dSb/dt
clear dSb onset
%zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);    %function to find zero-crossings   --NOT USED ANYMORE

        for ity=1:nyear-1 
            if(isnan(Sb(ity,:)))
                dSb(ity,:)=NaN;
            else
                dSb(ity,:)=diff(Sb(ity,:));   
            end
            if(isnan(Sb2(ity,:)))
                dSb2(ity,:)=NaN;
            else
                dSb2(ity,:)=diff(Sb2(ity,:));   
            end
        end

display('dSb Done!')

%Compute onset and demise
clear onset demise

        for ity=1:nyear-1 
            o=0;
            d=0;
            if(isnan(dSb(ity,:)) | isnan(dSb2(ity,:)))
                onset(ity)=NaN;
                demise(ity)=NaN;
            else
                for itp=2:71  %pentad relative to t0!!!
                    if o==1
                    else
                        if ( dSb(ity,itp-1) < 0 && dSb(ity,itp) > 0 && dSb(ity,itp+1) > 0 )
                        onset(ity)=itp%+inipent;  %relative to Jan 1st!
                        onset4anom=onset;
                          if onset(ity)>73
                              onset(ity)=onset(ity)-73;
                          end
                        o=1;
                        else
                          onset(ity)=NaN;  
                        end
                    end
                    if d==1
                    else
                        if ( dSb2(ity,itp-1) > 0 && dSb2(ity,itp) < 0 && dSb2(ity,itp+1) < 0 )
                        demise(ity)=itp+inipent2;  %relative to Jan 1st!
                        demise4anom=demise;
                        %if ( dSb(ilo,ila,ity,itp-1) > 0 && dSb(ilo,ila,ity,itp) < 0  && dSb(ilo,ila,ity,itp+1) < 0 )
                        %demise(ilo,ila,ity)=itp+inipent(ilo,ila);  %relative to Jan 1st!
                          if demise(ity)>73
                              demise(ity)=demise(ity)-73;
                          end
                        d=1;
                        else
                          demise(ity)=NaN;
                        end
                    end
                end
            end
        end

display('Onset and Demise Done!')
save -v7.3 onsetdemise.mat onset onset4anom demise demise4anom inipent inipent2 S Sb dSb pentad 
else
load onsetdemise.mat onset onset4anom demise demise4anom inipent inipent2 S Sb dSb pentad 
end

%%
%TEST PLOTS
%Select a year:
iyr=4;
%
figure(2);clf %S,Sb,dSb for onset
pfig=plot(1:72,squeeze(S(iyr,1:72)),1:72,squeeze(Sb(iyr,1:72)),1:72,squeeze(dSb(iyr,1:72)),'LineWidth',3);hold on
bfig=bar(1:72,squeeze(pentad(1:72,iyr)'))   
plot(1:73,pentad(:,iyr),onset(iyr),pentad(onset(iyr)),'*k',demise(iyr),pentad(demise(iyr)),'*k','markers',12)
title('S,Sb,dSb for onset')
xlabel('pentad')
legend('S','Sb','dSb/dt','rainfall','zero line','onset','demise')
legend('Location','southwest')
legend('boxoff')

figure(3);clf %S,Sb,dSb for demise
pfig=plot(1:72,squeeze(S2(iyr,1:72)),1:72,squeeze(Sb2(iyr,1:72)),1:72,squeeze(dSb2(iyr,1:72)),'LineWidth',3);hold on
bfig=bar(1:72,squeeze(pentad(1:72,iyr)'))   
plot(1:73,pentad(:,iyr),onset(iyr),pentad(onset(iyr)),'*k',demise(iyr),pentad(demise(iyr)),'*k','markers',12)
title('S,Sb,dSb for demise')
xlabel('pentad')
legend('S','Sb','dSb/dt','rainfall','zero line','onset','demise')
legend('Location','southwest')
legend('boxoff')