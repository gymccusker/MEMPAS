%%
%
%       SIZE DISTRIBUTION FUNCTION
%

function [OUTPUT,zz]=sem_sizedist(OUTPUT,flag)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %     INITIALISE DATA
% % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


area_cm=OUTPUT.INITIAL.AREA_mm./100;
D_mm=OUTPUT.INITIAL.FILTERDIAM_mm-10;
D_cm=D_mm/10;        %   exposed filter diameter in cm
A_mm2=pi*(D_mm/2)^2;
A_cm2=pi*(D_cm/2)^2;

% init_flow=OUTPUT.INPUT.VOLUME_L*1e3;          %   total flow in cm3

% STP_T=273.15;        % STP temp in kelvin
% STP_P=100000;        % STP pressure in Pa
% STP_rho=STP_P./(287.05.*STP_T); % STP density in kg/m3
% mean_T=254.6961;            % filter averaged mean temp in K
% mean_P=983.9006.*100;     % filter averaged mean pressure in Pa
% mean_rho=mean_P./(287.05.*mean_T); 
% tot_flow=(init_flow./mean_rho).*STP_rho;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %     SIZE BINS
% % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


OUTPUT.All.SIZEBINS=[0 unique(OUTPUT.All.Raw(:,3)')]';
OUTPUT.SIZEDIST.SIZEBINS=OUTPUT.All.SIZEBINS;
n_tot=histc(OUTPUT.All.Raw(:,3),OUTPUT.All.SIZEBINS);
part_denstot=n_tot./OUTPUT.INITIAL.AREA_mm;
part_tot=part_denstot.*A_cm2;
% part_tot=nansum(tot,2)./nansum(tot~=0,2);               % for more than 1 scan, sum up number in each bin, ignoring zeros
OUTPUT.SIZEDIST.NUMCONC_cm3=part_tot./(OUTPUT.INITIAL.VOLUME_L.*1e3);       % number concentration (#/cc) in each size bin

% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % %  Calculate dN/dlogD
% % % % % % % % % % % % % % % % % % 

% OUTPUT.SIZEDIST.DNDLOGD_cm3=(OUTPUT.SIZEDIST.NUMCONC_cm3.*log(10).*OUTPUT.All.SIZEBINS)./(0.01);
edge_upp=OUTPUT.SIZEDIST.SIZEBINS(1:end-1)+(OUTPUT.SIZEDIST.SIZEBINS(2:end)-OUTPUT.SIZEDIST.SIZEBINS(1:end-1))/2;
edge_low=[0; OUTPUT.SIZEDIST.SIZEBINS(2:end-1)-(OUTPUT.SIZEDIST.SIZEBINS(3:end)-OUTPUT.SIZEDIST.SIZEBINS(2:end-1))/2];
OUTPUT.SIZEDIST.DNDLOGD_cm3=OUTPUT.SIZEDIST.NUMCONC_cm3(1:end-1,:)./(real(log10(edge_upp))-real(log10(edge_low)));
bin_width=edge_upp-edge_low;

% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % %  AVERAGE/SMOOTH
% % % % % % % % % % % % % % % % % % 

k=length(OUTPUT.SIZEDIST.DNDLOGD_cm3);   
for i=1:round(k/2)-1,
        if i<round(k/2)-1,
            OUTPUT.SIZEDIST.SIZEBINS_AV(i)=nanmean([OUTPUT.SIZEDIST.SIZEBINS(2*i-1) OUTPUT.SIZEDIST.SIZEBINS(2*i)]);
            OUTPUT.SIZEDIST.DNDLOGD_AV_cm3(i)=nanmean([OUTPUT.SIZEDIST.DNDLOGD_cm3(2*i-1) OUTPUT.SIZEDIST.DNDLOGD_cm3(2*i)]);
            OUTPUT.SIZEDIST.NUMCONC_AV_cm3(i)=nanmean([OUTPUT.SIZEDIST.NUMCONC_cm3(2*i-1) OUTPUT.SIZEDIST.NUMCONC_cm3(2*i)]);
        else
            OUTPUT.SIZEDIST.SIZEBINS_AV(i)=OUTPUT.SIZEDIST.SIZEBINS(2*i);
            OUTPUT.SIZEDIST.DNDLOGD_AV_cm3(i)=OUTPUT.SIZEDIST.DNDLOGD_cm3(2*i);
            OUTPUT.SIZEDIST.NUMCONC_AV_cm3(i)=OUTPUT.SIZEDIST.NUMCONC_cm3(2*i);
        end
end
% k=length(OUTPUT.SIZEDIST.DNDLOGD_cm3);   
% for i=1:round(k/15)-1,
%         if i<round(k/15),
%             OUTPUT.SIZEDIST.SIZEBINS_AV2(i)=nanmean([OUTPUT.SIZEDIST.SIZEBINS(15*i-14); OUTPUT.SIZEDIST.SIZEBINS(15*i-13);...
%                 OUTPUT.SIZEDIST.SIZEBINS(15*i-12); OUTPUT.SIZEDIST.SIZEBINS(15*i-11);...
%                 OUTPUT.SIZEDIST.SIZEBINS(15*i-10); OUTPUT.SIZEDIST.SIZEBINS(15*i-9); ...
%                 OUTPUT.SIZEDIST.SIZEBINS(15*i-8); OUTPUT.SIZEDIST.SIZEBINS(15*i-7); OUTPUT.SIZEDIST.SIZEBINS(15*i-6);...
%                 OUTPUT.SIZEDIST.SIZEBINS(15*i-5); OUTPUT.SIZEDIST.SIZEBINS(15*i-4); OUTPUT.SIZEDIST.SIZEBINS(15*i-3);...
%                 OUTPUT.SIZEDIST.SIZEBINS(15*i-2); OUTPUT.SIZEDIST.SIZEBINS(15*i-1); OUTPUT.SIZEDIST.SIZEBINS(15*i)]);
%             OUTPUT.SIZEDIST.DNDLOGD_AV2_cm3(i)=nanmean([OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-14); OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-13);...
%                 OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-12); OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-11);...
%                 OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-10); OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-9); OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-8);...
%                 OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-7); OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-6);...
%                 OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-5); OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-4); OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-3);...
%                 OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-2); OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i-1); OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i)]);
%         else
%             OUTPUT.SIZEDIST.SIZEBINS_AV2(i)=OUTPUT.SIZEDIST.SIZEBINS(15*i);
%             OUTPUT.SIZEDIST.DNDLOGD_AV2_cm3(i)=OUTPUT.SIZEDIST.DNDLOGD_cm3(15*i);
%         end
% end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % %  CONC OF CLASSIFICATIONS >0.5UM
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

if flag==1,
    classes={'Silicate','MixedSil','CaRich','FreshCl','AgedCl','Sulph',...
        'Gypsum','Carbon','Bio','Metal','Phos','Biomass','Other'};
elseif flag==2
    classes={'Silicate','Quartz','Carbonates','Chlorides','Oxides',...
        'Sulphates','Phosphates','Mixtures','Secondary','Other'};
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % %  SILICATES
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%         'OUTPUT.',char(classes(i)),'.DNDLOGD_cm3=OUTPUT.',char(classes(i)),'.NumberConcentration_cm3(1:end-1,:)./(real(log10(edge_upp))-real(log10(edge_low)));',...

for i=1:length(classes),
     eval(['if length(OUTPUT.',char(classes(i)),'.All_Indices)>1,',...
        'OUTPUT.',char(classes(i)),'.Number_on_filter=histc(OUTPUT.',char(classes(i)),'.Raw(:,3),OUTPUT.All.SIZEBINS);',...
        'OUTPUT.',char(classes(i)),'.ParticleDensity=OUTPUT.',char(classes(i)),'.Number_on_filter./area_cm;',...
        'OUTPUT.',char(classes(i)),'.ParticleNumber=OUTPUT.',char(classes(i)),'.ParticleDensity.*A_cm2;',...
        'OUTPUT.',char(classes(i)),'.NumberConcentration_cm3=OUTPUT.',char(classes(i)),'.ParticleNumber./(OUTPUT.INITIAL.VOLUME_L.*1e3);',...
            'for j=1:length(bin_width);',...
                'OUTPUT.',char(classes(i)),'.DNDLOGD_cm3=(OUTPUT.',char(classes(i)),'.NumberConcentration_cm3.*log(10).*OUTPUT.All.SIZEBINS)./(bin_width(j));',...       
            'end;',...
    'elseif length(OUTPUT.',char(classes(i)),'.All_Indices)==1,',...
        'OUTPUT.',char(classes(i)),'.Number_on_filter=histc(OUTPUT.',char(classes(i)),'.Raw(:,3),OUTPUT.All.SIZEBINS);',...
        'OUTPUT.',char(classes(i)),'.ParticleDensity=OUTPUT.',char(classes(i)),'.Number_on_filter./area_cm;',...
        'OUTPUT.',char(classes(i)),'.ParticleNumber=OUTPUT.',char(classes(i)),'.ParticleDensity.*A_cm2;',...
        'OUTPUT.',char(classes(i)),'.NumberConcentration_cm3=OUTPUT.',char(classes(i)),'.ParticleNumber./(OUTPUT.INITIAL.VOLUME_L.*1e3);',...    
     'end;'])
end


zz=figure(); hold on;
h1=plot(OUTPUT.SIZEDIST.SIZEBINS_AV,OUTPUT.SIZEDIST.DNDLOGD_AV_cm3,'x','linewidth',2,'markersize',6);
xlabel('Size [\mum]','fontsize',16)
ylabel('dN/dlog_{10}D_P [cm^{-3}]','fontsize',16)
set(gca,'xscale','log','yscale','log')
set(gca,'fontsize',16)
% plot(OUTPUT.SIZEDIST.SIZEBINS_AV2,OUTPUT.SIZEDIST.DNDLOGD_AV2_cm3,'k-','linewidth',2);
% legend('Raw (averaged 1x)','Averaged (averaged 15x)')


end



