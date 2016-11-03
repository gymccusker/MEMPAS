%%
%
%       CLASSIFICATION FUNCTION
%

function [OUTPUT]=sem_classifications(input,txt,final_entry)

% field=input(:,1);
% area=input(:,2);
partsize=input(:,3);
% volume=input(:,4);
n=length(input(:,1));
input(isnan(input))=0;

if length(find(strcmp('B', txt)==1))==1,  B=input(:,find(strcmp('B', txt)==1));   B(isnan(B))=0;   OUTPUT.RAW.B=B;   else B=zeros(n,1);  end
if length(find(strcmp('C', txt)==1))==1,  C=input(:,find(strcmp('C', txt)==1));   C(isnan(C))=0;   OUTPUT.RAW.C=C;   end
if length(find(strcmp('N', txt)==1))==1,  N=input(:,find(strcmp('N', txt)==1));   N(isnan(N))=0;   OUTPUT.RAW.N=N;   else N=zeros(n,1);  end
if length(find(strcmp('O', txt)==1))==1,  O=input(:,find(strcmp('O', txt)==1));   O(isnan(O))=0;   OUTPUT.RAW.O=O;   end
if length(find(strcmp('F', txt)==1))==1,  F=input(:,find(strcmp('F', txt)==1));   F(isnan(F))=0;   OUTPUT.RAW.F=F;   else F=zeros(n,1);  end
if length(find(strcmp('Na', txt)==1))==1, Na=input(:,find(strcmp('Na', txt)==1)); Na(isnan(Na))=0; OUTPUT.RAW.Na=Na; else Na=zeros(n,1); end
if length(find(strcmp('Mg', txt)==1))==1, Mg=input(:,find(strcmp('Mg', txt)==1)); Mg(isnan(Mg))=0; OUTPUT.RAW.Mg=Mg; else Mg=zeros(n,1); end
if length(find(strcmp('Al', txt)==1))==1, Al=input(:,find(strcmp('Al', txt)==1)); Al(isnan(Al))=0; OUTPUT.RAW.Al=Al; else Al=zeros(n,1); end
if length(find(strcmp('Si', txt)==1))==1, Si=input(:,find(strcmp('Si', txt)==1)); Si(isnan(Si))=0; OUTPUT.RAW.Si=Si; else Si=zeros(n,1); end
if length(find(strcmp('P', txt)==1))==1,  P=input(:,find(strcmp('P', txt)==1));   P(isnan(P))=0;   OUTPUT.RAW.P=P;   else P=zeros(n,1);  end
if length(find(strcmp('S', txt)==1))==1,  S=input(:,find(strcmp('S', txt)==1));   S(isnan(S))=0;   OUTPUT.RAW.S=S;   else S=zeros(n,1);  end
if length(find(strcmp('Cl', txt)==1))==1, Cl=input(:,find(strcmp('Cl', txt)==1)); Cl(isnan(Cl))=0; OUTPUT.RAW.Cl=Cl; else Cl=zeros(n,1); end
if length(find(strcmp('K', txt)==1))==1,  K=input(:,find(strcmp('K', txt)==1));   K(isnan(K))=0;   OUTPUT.RAW.K=K;   else K=zeros(n,1);  end
if length(find(strcmp('Ca', txt)==1))==1, Ca=input(:,find(strcmp('Ca', txt)==1)); Ca(isnan(Ca))=0; OUTPUT.RAW.Ca=Ca; else Ca=zeros(n,1); end
if length(find(strcmp('Sc', txt)==1))==1, Sc=input(:,find(strcmp('Sc', txt)==1)); Sc(isnan(Sc))=0; OUTPUT.RAW.Sc=Sc; else Sc=zeros(n,1); end
if length(find(strcmp('Ti', txt)==1))==1, Ti=input(:,find(strcmp('Ti', txt)==1)); Ti(isnan(Ti))=0; OUTPUT.RAW.Ti=Ti; else Ti=zeros(n,1); end
if length(find(strcmp('V', txt)==1))==1,  V=input(:,find(strcmp('V', txt)==1));   V(isnan(V))=0;   OUTPUT.RAW.V=V;   else V=zeros(n,1); end
if length(find(strcmp('Cr', txt)==1))==1, Cr=input(:,find(strcmp('Cr', txt)==1)); Cr(isnan(Cr))=0; OUTPUT.RAW.Cr=Cr; else Cr=zeros(n,1); end
if length(find(strcmp('Mn', txt)==1))==1, Mn=input(:,find(strcmp('Mn', txt)==1)); Mn(isnan(Mn))=0; OUTPUT.RAW.Mn=Mn; else Mn=zeros(n,1); end
if length(find(strcmp('Fe', txt)==1))==1, Fe=input(:,find(strcmp('Fe', txt)==1)); Fe(isnan(Fe))=0; OUTPUT.RAW.Fe=Fe; else Fe=zeros(n,1); end
if length(find(strcmp('Ni', txt)==1))==1, Ni=input(:,find(strcmp('Ni', txt)==1)); Ni(isnan(Ni))=0; OUTPUT.RAW.Ni=Ni; else Ni=zeros(n,1); end
if length(find(strcmp('Cu', txt)==1))==1, Cu=input(:,find(strcmp('Cu', txt)==1)); Cu(isnan(Cu))=0; OUTPUT.RAW.Cu=Cu; else Cu=zeros(n,1); end
if length(find(strcmp('Zn', txt)==1))==1, Zn=input(:,find(strcmp('Zn', txt)==1)); Zn(isnan(Zn))=0; OUTPUT.RAW.Zn=Zn; else Zn=zeros(n,1); end
if length(find(strcmp('Br', txt)==1))==1, Br=input(:,find(strcmp('Br', txt)==1)); Br(isnan(Br))=0; OUTPUT.RAW.Br=Br; else Br=zeros(n,1); end
if length(find(strcmp('Sr', txt)==1))==1, Sr=input(:,find(strcmp('Sr', txt)==1)); Sr(isnan(Sr))=0; OUTPUT.RAW.Sr=Sr; else Sr=zeros(n,1); end
if length(find(strcmp('Mo', txt)==1))==1, Mo=input(:,find(strcmp('Mo', txt)==1)); Mo(isnan(Mo))=0; OUTPUT.RAW.Mo=Mo; else Mo=zeros(n,1); end
if length(find(strcmp('Ru', txt)==1))==1, Ru=input(:,find(strcmp('Ru', txt)==1)); Ru(isnan(Ru))=0; OUTPUT.RAW.Ru=Ru; else Ru=zeros(n,1); end
if length(find(strcmp('Rh', txt)==1))==1, Rh=input(:,find(strcmp('Rh', txt)==1)); Rh(isnan(Rh))=0; OUTPUT.RAW.Rh=Rh; else Rh=zeros(n,1); end
if length(find(strcmp('Pd', txt)==1))==1, Pd=input(:,find(strcmp('Pd', txt)==1)); Pd(isnan(Pd))=0; OUTPUT.RAW.Pd=Pd; else Pd=zeros(n,1); end
if length(find(strcmp('Ba', txt)==1))==1, Ba=input(:,find(strcmp('Ba', txt)==1)); Ba(isnan(Ba))=0; OUTPUT.RAW.Ba=Ba; else Ba=zeros(n,1); end
if length(find(strcmp('W', txt)==1))==1,  W=input(:,find(strcmp('W', txt)==1));   W(isnan(W))=0;   OUTPUT.RAW.W=W;   else W=zeros(n,1);  end
if length(find(strcmp('Pb', txt)==1))==1, Pb=input(:,find(strcmp('Pb', txt)==1)); Pb(isnan(Pb))=0; OUTPUT.RAW.Pb=Pb; else Pb=zeros(n,1); end
  
elms=txt(7:final_entry);
OUTPUT.RAW.Elements=elms;
OUTPUT.RAW.Headers=txt;
unclassed=input; 

total=zeros(n,1);
for i=1:n,
   total(i)=sum(input(i,9:final_entry)); 
   all(i)=sum(input(i,7:final_entry)); 
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % %     BEGIN CLASSIFICATION SCHEME
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
 
%
%
%       Carbonaceous Particles
%
%

organic=zeros(n,length(txt));
for i=1:n,
    if unclassed(i,7:final_entry)<=0.2*total(i)
        if C(i) + O(i) > 0.92*all(i);
            organic(i,:)=input(i,:);
        else
            organic(i,:)=0;
        end
    else
        organic(i,:)=0;
    end
    
    if C(i)+O(i)>= 0.9*all(i);   % 90% carbon+O from Geng et al 2010
        if unclassed(i,7:final_entry)<=0.2*total(i)
            if Mg(i) >= 0.1*total(i)              %   Mg
                organic(i,:)=input(i,:);
                
            end
            if Na(i) >= 0.1*total(i)              %   Na
                organic(i,:)=input(i,:);
                
            end
            if S(i) >= 0.1*total(i)             %   S
                organic(i,:)=input(i,:);
            end
        end
        
    end
end
organic_index=find(organic(:,3)>0);

if isempty(organic_index)
    disp('No Carbonaceous Particles Detected')
else
    unclassed(organic_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(organic_index)=0;']); end 
end



% %
% %
% %           Secondary Particles
% %
% %

Narich=zeros(n,length(txt));
ammsulph=zeros(n,length(txt));
for i=1:n,
       
    %
    %           Na Rich
    %

    if Na(i)/total(i)>=0.2
        if Na(i)/total(i)<=1.1
            if Cl(i)/total(i)<0.02499
                if Mg(i)/Na(i)<1.1
                    if Al(i)/Na(i)<0.75
                        if Si(i)/Na(i)<0.25
                            if P(i)/Na(i)<0.1
                                if S(i)/Na(i)<0.1
                                    if K(i)/Na(i)<0.5
                                        if Ca(i)/Na(i)<0.5
                                            if Ti(i)/Na(i)<0.05
                                                if Cr(i)/Na(i)<0.05
                                                    if Fe(i)/Na(i)<0.1
                                                        Narich(i,:)=input(i,:); 
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

end
Narich_index=find(Narich(:,3)>0);
if isempty(Narich_index)
    disp('No "Na-Rich" Particles Detected')
else
    unclassed(Narich_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Narich_index)=0;']); end       
end

for i=1:n,
    
    %
    %       Ammonium Sulphate
    %                    
    
    if S(i)/total(i)>=0.3
        if S(i)/total(i)<=1.1
            if Na(i)/S(i)<0.1
                if Mg(i)/S(i)<0.1
                    if Al(i)/S(i)<0.2
                        if Si(i)/S(i)<0.25
                            if P(i)/S(i)<0.1
                                if Cl(i)/S(i)<0.1
                                    if K(i)/S(i)<0.1
                                        if Ca(i)/S(i)<0.1
                                            if Ti(i)/S(i)<0.05
                                                if Cr(i)/S(i)<0.05
                                                    if Fe(i)/S(i)<0.1
                                                        ammsulph(i,:)=input(i,:);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
ammsulph_index=find(ammsulph(:,3)>0);
if isempty(ammsulph_index)
    disp('No "Ammonium Sulphate" Particles Detected')
else
    unclassed(ammsulph_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(ammsulph_index)=0;']); end         
end

%
%
%       Biogenic Particles
%               e.g. pollen, bacteria
%                       Traces of N, K, P, S and/or Cl (Geng et al 2010)
%

bio=zeros(n,length(txt));

for i=1:n,
    
    % Kandler 2011
    if partsize(i)>0.3
        if (K(i)+Na(i)+S(i)+P(i)+Ca(i))/total(i) >= 0.4
            if (K(i)+Na(i)+S(i)+P(i)+Ca(i))/total(i) <= 1.1
                if P(i)/total(i)>=0.05
                    if P(i)/total(i)<=0.8
                        if Na(i)/total(i) >=0.05
                            if Na(i)/total(i) <=0.8
                                if Ca(i)/total(i)>=0.05
                                    if Ca(i)/total(i)<=1.1
                                        if K(i)/total(i)>=0.025
                                            if K(i)/total(i)<=0.8
                                                if S(i)/total(i)>=0.025
                                                    if S(i)/total(i)<=0.8
                                                        if Mg(i)/(Na(i)+P(i)+S(i)+Ca(i))<0.1
                                                            if Al(i)/(Na(i)+P(i)+S(i)+Ca(i))<0.05
                                                                if Si(i)/(Na(i)+P(i)+S(i)+Ca(i))<0.1
                                                                    if Cl(i)/(Na(i)+P(i)+S(i)+Ca(i))<0.05
                                                                        if Ti(i)/(Na(i)+P(i)+S(i)+Ca(i))<0.05
                                                                            if Cr(i)/(Na(i)+P(i)+S(i)+Ca(i))<0.05
                                                                                if Fe(i)/(Na(i)+P(i)+S(i)+Ca(i))<0.1
                                                                                    bio(i,:)=input(i,:);
                                                                                else bio(i,:)=0;
                                                                                end
                                                                            else bio(i,:)=0;
                                                                            end
                                                                        else bio(i,:)=0;
                                                                        end
                                                                    else bio(i,:)=0;
                                                                    end
                                                                else bio(i,:)=0;
                                                                end
                                                            else bio(i,:)=0;
                                                            end
                                                        else bio(i,:)=0;
                                                        end
                                                    else bio(i,:)=0;
                                                    end
                                                else bio(i,:)=0;
                                                end
                                            else bio(i,:)=0;
                                            end
                                        else bio(i,:)=0;
                                        end
                                    else bio(i,:)=0;
                                    end
                                    
                                else bio(i,:)=0;
                                end
                            else bio(i,:)=0;
                            end
                        else bio(i,:)=0;
                        end
                    else bio(i,:)=0;
                    end
                else bio(i,:)=0;
                end
            else bio(i,:)=0;
            end
        else bio(i,:)=0;
        end
    else bio(i,:)=0;
    end
    
    thresh=0.9*all(i);            % 90% C+O from Geng et al 2010
    if C(i)+O(i) >= thresh
        if Si(i)/Cl(i) < 0.2
            if Na(i)/Cl(i)<0.3           % Make sure not sea salt
                if Na(i)/S(i)<0.3           % Make sure not sulphate
                    if K(i) >= 0.2*total(i)         % K
                        bio(i,:)=input(i,:);
                    end
                    if P(i) >= 0.2*total(i)         % P
                        bio(i,:)=input(i,:);
                    end
                    if S(i) >= 0.2*total(i)         % S
                        bio(i,:)=input(i,:);
                    end
                    if Cl(i) >= 0.2*total(i)         % Cl
                        bio(i,:)=input(i,:);
                    end
                    if Ca(i) + K(i)>=0.3*total(i)           % Ca + K >30%        % Krejci et al 2005
                        bio(i,:)=input(i,:);
                    end
                    if Na(i)+P(i)+K(i)>= 0.3*total(i)    % Na + P + K >30%
                        bio(i,:)=input(i,:);
                    end
                    if Na(i)+Mg(i)+Zn(i)>=0.3*total(i)        % Na + Mg + Zn>30%
                        bio(i,:)=input(i,:);
                    end
                end
            end
        end
    end

end


bio_index=find(bio(:,3)>0);

if isempty(bio_index)
    disp('No Biogenic Particles Detected')
else
    unclassed(bio_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(bio_index)=0;']); end        
end


%
%
%           Sulphates
%
%

NaS=zeros(n,length(txt));
for i=1:n,

    if Na(i)/S(i)>=0.101
        if Na(i)/S(i)<=10
            if (Na(i)+S(i))/total(i)>=0.1
                if (Na(i)+S(i))/total(i)<=1.1
                    if S(i)/total(i)>=0.025
                        if S(i)/total(i)<=1.1
                            if Mg(i)/(Na(i)+S(i))<0.5
                                if Al(i)/(Na(i)+S(i))<0.1
                                    if Si(i)/(Na(i)+S(i))<0.15
                                        if P(i)/(Na(i)+S(i))<0.5
                                            if Cl(i)/(Na(i)+S(i))<0.1
                                                if K(i)/(Na(i)+S(i))<0.1
                                                    if Ca(i)/(Na(i)+S(i))<0.05
                                                        if Ti(i)/(Na(i)+S(i))<0.05
                                                            if Cr(i)/(Na(i)+S(i))<0.05
                                                                if Fe(i)/(Na(i)+S(i))<0.1
                                                                    NaS(i,:)=input(i,:);
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
NaS_index=find(NaS(:,3)>0);
if isempty(NaS_index)
    disp('No "Na Sulphate" Particles Detected')
else
    unclassed(NaS_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(NaS_index)=0;']); end      
end  


%
%       Gypsum Particles
%

gypsum=zeros(n,length(txt));

for i=1:n,

    thresh=0.5*total(i);
    if Ca(i) + S(i) > thresh
        if Ca(i)/S(i) <= 4
            if Ca(i)/S(i) >= 0.25
                if Si(i)/Ca(i) < 0.5
                    gypsum(i,:)=input(i,:);
                else
                    gypsum(i,:)=0;
                end
            else
                gypsum(i,:)=0;
            end
        else
            gypsum(i,:)=0;
        end
    else
        gypsum(i,:)=0;
    end
    
    
    %   OR
    
   if Ca(i)/S(i)>=0.2
       if Ca(i)/S(i)<=10
           if (Ca(i)+S(i))/total(i)>=0.2
               if (Ca(i)+S(i))/total(i)<=1.1
                   if Na(i)/(Ca(i)+S(i))<0.1
                       if Mg(i)/(Ca(i)+S(i))<0.35
                           if Al(i)/(Ca(i)+S(i))<0.1
                               if Si(i)/(Ca(i)+S(i))<0.1
                                   if P(i)/(Ca(i)+S(i))<0.1
                                       if Cl(i)/(Ca(i)+S(i))<0.1
                                           if K(i)/(Ca(i)+S(i))<0.1
                                               if Ti(i)/(Ca(i)+S(i))<0.05
                                                   if Cr(i)/(Ca(i)+S(i))<0.05
                                                       if Fe(i)/(Ca(i)+S(i))<0.1
                                                           gypsum(i,:)=input(i,:);
                                                       end
                                                   end
                                               end
                                           end
                                       end
                                   end
                               end
                           end
                       end
                   end
               end
           end
       end
   end
end


gypsum_index=find(gypsum(:,3)>0);
if isempty(gypsum_index)
    disp('No Gypsum Particles Detected')
else
    unclassed(gypsum_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(gypsum_index)=0;']); end               
end

%
%           Sulphur-Rich (part of Sulphate group)
%

Srich=zeros(n,length(txt));

for i=1:n,

    thresh=0.4*total(i);
    if S(i)>thresh
        if  Si(i)/S(i)<0.5
            if S(i)>Na(i)
                if S(i)>Mg(i)
                    if S(i)>Al(i)
                        if S(i)>Si(i)
                            if S(i)>P(i)
                                if S(i)>Cl(i)
                                    if S(i)>K(i)
                                        if S(i)>Ca(i)
                                            if S(i)>Ti(i)
                                                if S(i)>Cr(i)
                                                    if S(i)>Mn(i)
                                                        if S(i)>Fe(i)
                                                            if S(i)>Ni(i)
                                                                if S(i)>Cu(i)
                                                                    if S(i)>Zn(i)
                                                                        Srich(i,:)=input(i,:);
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

Srich_index=find(Srich(:,3)>0);
if isempty(Srich_index)
    disp('No "S-Rich" Particles Detected')
else
    unclassed(Srich_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Srich_index)=0;']); end          
end 

CaNaS=zeros(n,length(txt));
for i=1:n,
 
   %
   %            CaNa Sulphates
   %
   
   if (Na(i)+S(i)+Ca(i))/total(i)>=0.15
       if (Na(i)+S(i)+Ca(i))/total(i)<=1.1
           if Na(i)/total(i)>=0.025
               if Na(i)/total(i)<=1.1
                   if S(i)/total(i)>=0.025
                       if S(i)/total(i)<=1.1
                           if Ca(i)/total(i)>=0.025
                               if Ca(i)/total(i)<=1.1
                                   if Na(i)/Ca(i)>=0.1
                                       if Na(i)/Ca(i)<=10
                                           if Mg(i)/(Na(i)+S(i)+Ca(i))<0.5
                                               if Al(i)/(Na(i)+S(i)+Ca(i))<0.05
                                                   if Si(i)/(Na(i)+S(i)+Ca(i))<0.05
                                                       if P(i)/(Na(i)+S(i)+Ca(i))<0.2
                                                           if Cl(i)/(Na(i)+S(i)+Ca(i))<0.1
                                                               if K(i)/(Na(i)+S(i)+Ca(i))<0.1
                                                                   if Ti(i)/(Na(i)+S(i)+Ca(i))<0.05
                                                                       if Cr(i)/(Na(i)+S(i)+Ca(i))<0.1
                                                                           if Fe(i)/(Na(i)+S(i)+Ca(i))<0.1
                                                                               if Ca(i)/(Na(i)+S(i))>=0.1001
                                                                                   if Ca(i)/(Na(i)+S(i))<=10
                                                                                       CaNaS(i,:)=input(i,:);
                                                                                   end
                                                                               end
                                                                           end
                                                                       end
                                                                   end
                                                               end
                                                           end
                                                       end
                                                   end
                                               end
                                           end
                                       end
                                   end
                               end
                           end
                       end
                   end
               end
           end
       end
   end
   
                                                        
end
CaNaS_index=find(CaNaS(:,3)>0);
if isempty(CaNaS_index)
    disp('No "CaNa Sulphate" Particles Detected')
else
    unclassed(CaNaS_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(CaNaS_index)=0;']); end      
end   



%
%
%           Calcium-Rich
%
%

%
%
%           Ca-Rich Particles
%
%

Carich=zeros(n,length(txt));

for i=1:n,

    if Ca(i) > 0.5*total(i)                  % Ca > 60%, little Si or Al
        if Si(i)/Ca(i) < 0.5
            if Al(i)/Ca(i) < 0.5
                if Ca(i)>Na(i)
                    if Ca(i)>Mg(i)
                        if Ca(i)>Al(i)
                            if Ca(i)>Si(i)
                                if Ca(i)>P(i)
                                    if Ca(i)>Cl(i)
                                        if Ca(i)>K(i)
                                            if Ca(i)>S(i)
                                                if Ca(i)>Ti(i)
                                                    if Ca(i)>Cr(i)
                                                        if Ca(i)>Mn(i)
                                                            if Ca(i)>Fe(i)
                                                                if Ca(i)>Ni(i)
                                                                    if Ca(i)>Cu(i)
                                                                        if Ca(i)>Zn(i)
                                                                            Carich(i,:)=input(i,:);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

Ca_index=find(Carich(:,3)>0);
if isempty(Ca_index)
    disp('No Calcium-Rich Particles Detected')
else
    unclassed(Ca_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Ca_index)=0;']); end          
end

%
%
%       Carbonate Particles
%
%

Cacarb=zeros(n,length(txt));
Camgcarb=zeros(n,length(txt));

for i=1:n,

    if Ca(i)/total(i)>=0.2
        if Ca(i)/total(i)<=1.1
            if Na(i)/Ca(i)<0.11
                if Mg(i)/Ca(i)<0.5
                    if Al(i)/Ca(i)<0.151
                        if Si(i)/Ca(i)<0.11
                            if P(i)/Ca(i)<0.1
                                if S(i)/Ca(i)<0.1
                                    if Cl(i)/Ca(i)<0.1
                                        if K(i)/Ca(i)<0.1
                                            if Ti(i)/Ca(i)<0.1
                                                if Cr(i)/Ca(i)<0.05
                                                    if Fe(i)/Ca(i)<0.1
                                                        Cacarb(i,:)=input(i,:);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
Cacarb_index=find(Cacarb(:,3)>0);
if isempty(Cacarb_index)
    disp('No "Ca Carbonate" Particles Detected')
else
    unclassed(Cacarb_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Cacarb_index)=0;']); end         
end  


for i=1:n,

    %
    %
    %       Ca Mg Carbonate
    %
    
    if (Ca(i)+Mg(i))/total(i)>=0.4
        if (Ca(i)+Mg(i))/total(i)<=1.1
            if Mg(i)/Ca(i)>=0.501
                if Mg(i)/Ca(i)<=2
                    if Na(i)/(Mg(i)+Ca(i))<0.5
                        if Al(i)/(Mg(i)+Ca(i))<0.1
                            if Si(i)/(Mg(i)+Ca(i))<0.2
                                if P(i)/(Mg(i)+Ca(i))<0.1
                                    if S(i)/(Mg(i)+Ca(i))<0.1
                                        if Cl(i)/(Mg(i)+Ca(i))<0.1
                                            if Ti(i)/(Mg(i)+Ca(i))<0.1
                                                if Cr(i)/(Mg(i)+Ca(i))<0.05
                                                    if Fe(i)/(Mg(i)+Ca(i))<0.1
                                                        Camgcarb(i,:)=input(i,:);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
end  
Camgcarb_index=find(Camgcarb(:,3)>0);
if isempty(Camgcarb_index)
    disp('No "CaMg Carbonate" Particles Detected')
else
    unclassed(Camgcarb_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Camgcarb_index)=0;']); end              
end  


carbonate=zeros(n,length(txt));

for i=1:n,

    thresh=0.5*total(i);
    if Ca(i)+Mg(i) > thresh
        if Si(i)/Ca(i) < 0.5
            if S(i)/Ca(i) < 0.25
                if P(i)/Ca(i) < 0.15;
                    if Mg(i)/Ca(i) <= 3
                        if Mg(i)/Ca(i) >= 0.33
                            carbonate(i,:)=input(i,:);
                        else
                            carbonate(i,:)=0;
                        end
                    else
                        carbonate(i,:)=0;
                    end
                else
                    carbonate(i,:)=0;
                end
            else
                carbonate(i,:)=0;
            end
        else
            carbonate(i,:)=0;
        end
    else
        carbonate(i,:)=0;
    end
    
end


carbonate_index=find(carbonate(:,3)>0);
if isempty(carbonate_index)
    disp('No Carbonate Particles Detected')
else
    unclassed(carbonate_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(carbonate_index)=0;']); end        
end
%
%
%               Phosphates
%
%


phos=zeros(n,length(txt));

for i=1:n

    if unclassed(i,4)>0,
        if P(i)/total(i)>=0.05
            if P(i)/total(i)<=1.1
                if Al(i)/(Ca(i)+P(i))<0.2
                    if Si(i)/(Ca(i)+P(i))<0.1
                        phos(i,:)=input(i,:);
                    end
                end
            end
        end
    end
    
    if unclassed(i,4)>0,
        if P(i)>=0.3*total(i)
            phos(i,:)=input(i,:);
        end
    end
end

phos_index=find(phos(:,3)>0);
if isempty(phos_index)
    disp('No "Phosphate" Particles Detected')
else
    unclassed(phos_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(phos_index)=0;']); end          
end  


%
%
%      Fresh Chlorides Particles
%
%

freshss=zeros(n,length(txt));

for i=1:n,

    thresh=0.5*total(i);
    if Na(i) + Cl(i) > thresh
        if S(i)/Na(i) < 0.375
            if S(i)/Cl(i) < 0.5
                if Si(i)/Cl(i) < 0.2
                    if Fe(i)/Cl(i) < 0.5
                        freshss(i,:)=input(i,:);
                    else
                        freshss(i,:)=0;
                    end
                else
                    freshss(i,:)=0;
                end
            else
                freshss(i,:)=0;
            end
        else
            freshss(i,:)=0;
        end
    else
        freshss(i,:)=0;
    end
    
    if Na(i)/Cl(i) <= 1.5                  %   Na/Cl ~ 1
        if Na(i)/Cl(i) >= 0.5
            if S(i)/Cl(i) < 0.5
                if Si(i)/Cl(i) < 0.2
                    if S(i)/Na(i) < 0.375
                        if Fe(i)/Cl(i) < 0.5
                            freshss(i,:)=input(i,:);
                        end
                    end
                end
            end
        end
    end
    
end

freshss_index=find(freshss(:,3)>0);
if isempty(freshss_index)
    disp('No Fresh Chlorides Particles Detected')
else
    unclassed(freshss_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(freshss_index)=0;']); end          
end

%
%           Chlorides from Kandler2011
%
%
NaCl=zeros(n,length(txt));
KCl=zeros(n,length(txt));
oCl=zeros(n,length(txt));


for i=1:n,

    if (Na(i)+Cl(i))/total(i)>=0.22
        if (Na(i)+Cl(i))/total(i)<=1.1
            if Na(i)/total(i)>=0.01
                if Na(i)/total(i)<=1.1
                    if Cl(i)/total(i)>=0.01
                        if Cl(i)/total(i)<=1.1
                            if Si(i)/total(i)<0.0499
                                if Al(i)/total(i)<0.0299
                                    if Mg(i)/(Na(i)+Cl(i))<2
                                        if P(i)/(Na(i)+Cl(i))<0.2
                                            if S(i)/(Na(i)+Cl(i))<0.25
                                                if K(i)/(Na(i)+Cl(i))<0.15
                                                    if Ti(i)/(Na(i)+Cl(i))<0.25
                                                        if Cr(i)/(Na(i)+Cl(i))<0.25
                                                            if Fe(i)/(Na(i)+Cl(i))<0.25
                                                                NaCl(i,:)=input(i,:);
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
NaCl_index=find(NaCl(:,3)>0);
if isempty(NaCl_index)
    disp('No "NaCl" Particles Detected')
else
    unclassed(NaCl_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(NaCl_index)=0;']); end        
end  

for i=1:n,

    
    %
    %       Potassium Chloride
    %
    
    if (K(i)+Cl(i))/total(i)>=0.3
        if (K(i)+Cl(i))/total(i)<=1.1
            if Na(i)/total(i)>=0.01
                if Na(i)/total(i)<=1.1
                    if Cl(i)/total(i)>=0.01
                        if Cl(i)/total(i)<=1.1
                            if Na(i)/(K(i)+Cl(i))<0.15
                                if Mg(i)/(K(i)+Cl(i))<0.1
                                    if Al(i)/(K(i)+Cl(i))<2
                                        if Si(i)/(K(i)+Cl(i))<0.25
                                            if P(i)/(K(i)+Cl(i))<0.2
                                                if S(i)/(K(i)+Cl(i))<0.25
                                                    if Ca(i)/(K(i)+Cl(i))<0.5
                                                        if Ti(i)/(K(i)+Cl(i))<0.25
                                                            if Cr(i)/(K(i)+Cl(i))<0.25
                                                                if Fe(i)/(Na(i)+Cl(i))<0.25
                                                                    KCl(i,:)=input(i,:);
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
KCl_index=find(KCl(:,3)>0);
if isempty(KCl_index)
    disp('No "KCl" Particles Detected')
else
    unclassed(KCl_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(KCl_index)=0;']); end    
end 

for i=1:n,

    %
    %           Other Chlorides
    %
    %
        
    if Cl(i)/total(i)>=0.25
        if Cl(i)/total(i)<=1.1
            if Si(i)/total(i)<0.0699
                if Al(i)/total(i)<0.0099
                    if Na(i)/Cl(i)<2
                        if Mg(i)/Cl(i)<2
                            if P(i)/Cl(i)<0.2
                                if S(i)/Cl(i)<0.2
                                    if K(i)/Cl(i)<2
                                        if Ca(i)/Cl(i)<2
                                            if Ti(i)/Cl(i)<0.1
                                                if Cr(i)/Cl(i)<0.1
                                                    if Fe(i)/Cl(i)<10
                                                        oCl(i,:)=input(i,:);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end 
oCl_index=find(oCl(:,3)>0);
if isempty(oCl_index)
    disp('No "Other Cl" Particles Detected')
else
    unclassed(oCl_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(oCl_index)=0;']); end
end  



%
%
%       Aged Chlorides Particles
%
%

agedss=zeros(n,length(txt));
mixClS=zeros(n,length(txt));

for i=1:n,

    thresh=0.4*total(i);
    if Na(i) + Cl(i) > thresh
        if S(i)/Na(i) < 5
            if S(i)/Cl(i) < 5
                if Si(i)/Cl(i) < 0.5
                    agedss(i,:)=input(i,:);
                else
                    agedss(i,:)=0;
                end
            else
                agedss(i,:)=0;
            end
        else
            agedss(i,:)=0;
        end
    else
        agedss(i,:)=0;
    end
    
   
    %
    %           Other Chlorides (adapted from K11 to include greater S
    %           fraction)
    %
    %
        
    if Cl(i)/total(i)>=0.1
        if Cl(i)/total(i)<=1.1
            if Si(i)/total(i)<=0.0699
                if Al(i)/total(i)<=0.0099
                    if Na(i)/Cl(i)<2
                        if Mg(i)/Cl(i)<2
                            if P(i)/Cl(i)<0.2
%                                 if S(i)/Cl(i)>0.1
                                    if K(i)/Cl(i)<2
                                        if Ca(i)/Cl(i)<2
                                            if Ti(i)/Cl(i)<0.1
                                                if Cr(i)/Cl(i)<0.1
                                                    if Fe(i)/Cl(i)<10
                                                        agedss(i,:)=input(i,:);
                                                    end
                                                end
                                            end
                                        end
                                    end
%                                 end
                            end
                        end
                    end
                end
            end
        end
    end
end    

agedss_index=find(agedss(:,3)>0);
if isempty(agedss_index)
    disp('No Aged Chlorides Particles Detected')
else
    unclassed(agedss_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(agedss_index)=0;']); end    
end

    %
    %           mixClS
    %
for i=1:n,   
    if Cl(i)/S(i)>=0.201
        if Cl(i)/S(i)<=10
            if (S(i)+Cl(i))/total(i)>=0.2
                if (S(i)+Cl(i))/total(i)<=1.1
                    if Cl(i)/total(i)>=0.025
                        if Cl(i)/total(i)<=1.1
                            if S(i)/total(i)>=0.025
                                if S(i)/total(i)<=1.1
                                    if S(i)/(Na(i)+Cl(i))>=0.1
                                        if S(i)/(Na(i)+Cl(i))<=20
                                            if Na(i)/(S(i)+Cl(i))<3
                                                if Mg(i)/(S(i)+Cl(i))<3
                                                    if Al(i)/(S(i)+Cl(i))<0.2
                                                        if Si(i)/(S(i)+Cl(i))<0.25
                                                            if P(i)/(S(i)+Cl(i))<0.25
                                                                if K(i)/(S(i)+Cl(i))<3
                                                                    if Fe(i)/(S(i)+Cl(i))<2
                                                                        mixClS(i,:)=input(i,:);
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
mixClS_index=find(mixClS(:,3)>0);
if isempty(mixClS_index)
    disp('No "ClS Mixtures" Particles Detected')
else
    unclassed(mixClS_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(mixClS_index)=0;']); end    
end 
 



%
%
%       Other Sulphate Particles
%
%

sulph=zeros(n,length(txt));

for i=1:n,

    if Cl(i)/S(i) < 0.5
        if Si(i)/S(i) < 0.5
            if Ti(i)/S(i) < 0.2
                if Cr(i)/S(i) < 0.2
                    if Fe(i)/S(i) < 0.5
                        if Ni(i)/S(i) < 0.2
                            if Cu(i)/S(i) < 0.2
                                if Zn(i)/S(i) < 0.2
                                    sulph(i,:)=input(i,:);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    

    %
    %           OR
    %
    
   %
   %            Other Sulphates (K11)
   %
   %
   
   if S(i)/total(i)>=0.2
       if S(i)/total(i)<=1.1
           if Na(i)/S(i)<2
               if Mg(i)/S(i)<2
                   if Al(i)/S(i)<2.5
                       if Si(i)/S(i)<0.25
                           if P(i)/S(i)<0.2
                               if Cl(i)/S(i)<0.2
                                   if K(i)/S(i)<10
                                       if Ca(i)/S(i)<2
                                           if Ti(i)/S(i)<0.5
                                               if Cr(i)/S(i)<0.5
                                                   if Fe(i)/S(i)<2
                                                       sulph(i,:)=input(i,:);
                                                   end
                                               end
                                           end
                                       end
                                   end
                               end
                           end
                       end
                   end
               end
           end
       end
   end
   
end
sulph_index=find(sulph(:,3)>0);
if isempty(sulph_index)
    disp('No Sulphate Particles Detected')
else
    unclassed(sulph_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(sulph_index)=0;']); end    
end

%
%
%               Oxides
%
%


FeO=zeros(n,length(txt));
TiO=zeros(n,length(txt));
FeTiO=zeros(n,length(txt));
AlO=zeros(n,length(txt));

for i=1:n,

    if Fe(i)/total(i)>=0.25
        if Fe(i)/total(i)<=1.1
            if Na(i)/Fe(i)<0.1
                if Mg(i)/Fe(i)<0.25
                    if Al(i)/Fe(i)<0.2
                        if Si(i)/Fe(i)<0.25
                            if P(i)/Fe(i)<0.2
                                if S(i)/Fe(i)<0.2
                                    if Cl(i)/Fe(i)<0.1
                                        if K(i)/Fe(i)<0.1
                                            if Ca(i)/Fe(i)<0.1
                                                if Ti(i)/Fe(i)<0.25
                                                    if Cr(i)/Fe(i)<0.05
                                                        FeO(i,:)=input(i,:);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
end
FeO_index=find(FeO(:,3)>0);
if isempty(FeO_index)
    disp('No "Iron Oxide" Particles Detected')
else
    unclassed(FeO_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(FeO_index)=0;']); end   
end  

for i=1:n,

    %
    %           Ti Oxide
    %

    if Ti(i)/total(i)>=0.25
        if Ti(i)/total(i)<=1.1
            if Na(i)/Ti(i)<0.18
                if Mg(i)/Ti(i)<0.1
                    if Al(i)/Ti(i)<0.2
                        if Si(i)/Ti(i)<0.25
                            if P(i)/Ti(i)<0.2
                                if S(i)/Ti(i)<0.2
                                    if Cl(i)/Ti(i)<0.1
                                        if K(i)/Ti(i)<0.1
                                            if Ca(i)/Ti(i)<0.1
                                                if Cr(i)/Ti(i)<0.05
                                                    if Fe(i)/Ti(i)<0.25
                                                        TiO(i,:)=input(i,:);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
TiO_index=find(TiO(:,3)>0);
if isempty(TiO_index)
    disp('No "Titanium Oxide" Particles Detected')
else
    unclassed(TiO_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(TiO_index)=0;']); end   
end  


for i=1:n,

    %
    %           Fe-Ti Oxides
    %
    
    if Ti(i)/Fe(i)>=0.2501
        if Ti(i)/Fe(i)<=4
            if (Ti(i)+Fe(i))/total(i)>=0.25
                if (Ti(i)+Fe(i))/total(i)<=1.1
                    if Na(i)/(Ti(i)+Fe(i))<0.2
                        if Mg(i)/(Ti(i)+Fe(i))<0.1
                            if Al(i)/(Ti(i)+Fe(i))<0.2
                                if Si(i)/(Ti(i)+Fe(i))<0.25
                                    if P(i)/(Ti(i)+Fe(i))<0.2
                                        if S(i)/(Ti(i)+Fe(i))<0.2
                                            if Cl(i)/(Ti(i)+Fe(i))<0.1
                                                if K(i)/(Ti(i)+Fe(i))<0.1
                                                    if Ca(i)/(Ti(i)+Fe(i))<0.1
                                                        if Cr(i)/(Ti(i)+Fe(i))<0.05
                                                            FeTiO(i,:)=input(i,:);
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
FeTiO_index=find(FeTiO(:,3)>0);
if isempty(FeTiO_index)
    disp('No "Fe-Ti Oxide" Particles Detected')
else
    unclassed(FeTiO_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(FeTiO_index)=0;']); end    
end  

for i=1:n,

    %
    %           Al Oxide
    %
    
    if Al(i)/total(i)>=0.2
        if Al(i)/total(i)<=1.1
            if Na(i)/Al(i)<0.2
                if Mg(i)/Al(i)<0.1
                    if Si(i)/Al(i)<0.2499
                        if P(i)/Al(i)<0.2
                            if S(i)/Al(i)<0.2
                                if Cl(i)/Al(i)<0.1
                                    if K(i)/Al(i)<0.1
                                        if Ca(i)/Al(i)<0.1
                                            if Ti(i)/Al(i)<0.1
                                                if Fe(i)/Al(i)<1
                                                    AlO(i,:)=input(i,:);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
AlO_index=find(AlO(:,3)>0);
if isempty(AlO_index)
    disp('No "Al Oxide" Particles Detected')
else
    unclassed(AlO_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(AlO_index)=0;']); end    
end  

%
%
%       Iron-Rich Particles
%
%

Fe=zeros(n,length(txt));

for i=1:n,

    thresh=0.3*total(i);
    if Fe(i) > thresh
        if Si(i)/Fe(i) < 0.2
            if Al(i)/Fe(i) < 0.2
                if Cl(i)/Fe(i) < 0.2
                    if Mg(i)/Fe(i) < 0.2
                        if Ca(i)/Fe(i) < 0.2
                            if Ti(i)/Fe(i) < 1.33
                                Fe(i,:)=input(i,:);
                            else
                                Fe(i)=0;
                            end
                        else
                            Fe(i)=0;
                        end
                    else
                        Fe(i)=0;
                    end
                else
                    Fe(i)=0;
                end
            else
                Fe(i)=0;
            end
        else
            Fe(i)=0;
        end
    end
    
end


Fe_index=find(Fe(:,3)>0);
if isempty(Fe_index)
    disp('No Iron-Rich Particles Detected')
else
    unclassed(Fe_index,:)=0;
     for j=1:length(elms), eval([char(elms(j)),'(Fe_index)=0;']); end   
end


%
%
%       Titanium-Rich Particles
%
%

Tirich=zeros(n,length(txt));

for i=1:n,

    thresh=0.3*total(i);
    if Ti(i) > thresh
        if Na(i)/Ti(i) < 1
            if Mg(i)/Ti(i) < 1
                if Al(i)/Ti(i) < 0.2
                    if Si(i)/Ti(i) < 0.2
                        if S(i)/Ti(i) < 1
                            if Fe(i)/Ti(i) < 1
                                Tirich(i,:)=input(i,:);
                            end
                        end
                    end
                end
            end
        end
    end

end


Ti_index=find(Tirich(:,3)>0);
if isempty(Ti_index)
    disp('No Titanium-Rich Particles Detected')
else
    unclassed(Ti_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Ti_index)=0;']); end   
end



%
%
%       MINERAL DUSTS
%
%
%

%
%
%       Quartz Particles
%
%

quartz=zeros(n,length(txt));

for i=1:n

    if partsize(i)>0.2
        if Si(i)/total(i)>=0.4
            if Si(i)/total(i)<=1.1
                if Al(i)/Si(i)<0.2
                    if Na(i)/Si(i)<0.1
                        if Mg(i)/Si(i)<0.1
                            if P(i)/Si(i)<0.2
                                if S(i)/Si(i)<0.2
                                    if Cl(i)/Si(i)<0.05
                                        if K(i)/Si(i)<0.1
                                            if Ca(i)/Si(i)<0.05
                                                if Ti(i)/Si(i)<0.1
                                                    if Cr(i)/Si(i)<0.05
                                                        if Fe(i)/Si(i)<0.1
                                                            quartz(i,:)=input(i,:);
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end

quartz_index=find(quartz(:,3)>0);
if isempty(quartz_index)
    disp('No "Quartz" Particles Detected')
else
    unclassed(quartz_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(quartz_index)=0;']); end    
end  


%
%
%               Silicates (K2011 criteria)
%
%

SiAl=zeros(n,length(txt));
SiAlK=zeros(n,length(txt));
SiAlNa=zeros(n,length(txt));
SiAlNaCa=zeros(n,length(txt));
SiAlNaK=zeros(n,length(txt));
SiAlCaFeMg=zeros(n,length(txt));
SiAlKFeMg=zeros(n,length(txt));
SiAlFeMg=zeros(n,length(txt));
SiMgFe=zeros(n,length(txt));
SiMg=zeros(n,length(txt));
SiCaTi=zeros(n,length(txt));

for i=1:n,

    if Al(i)/Si(i)>=0.201
        if Al(i)/Si(i)<=4
            if (Al(i)+Si(i))/total(i)>=0.4
                if (Al(i)+Si(i))/total(i)<=1.1
                    if Al(i)/total(i)>=0.05
                        if Al(i)/total(i)<=1.1
                            if Na(i)/(Al(i)+Si(i))<0.05
                                if Mg(i)/(Al(i)+Si(i))<0.05
                                    if P(i)/(Al(i)+Si(i))<0.2
                                        if S(i)/(Al(i)+Si(i))<0.2
                                            if Cl(i)/(Al(i)+Si(i))<0.1
                                                if K(i)/(Al(i)+Si(i))<0.05
                                                    if Ca(i)/(Al(i)+Si(i))<0.05
                                                        if Ti(i)/(Al(i)+Si(i))<0.1
                                                            if Cr(i)/(Al(i)+Si(i))<0.1
                                                                if Fe(i)/(Al(i)+Si(i))<0.1
                                                                    SiAl(i,:)=input(i,:);
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
SiAl_index=find(SiAl(:,3)>0);
if isempty(SiAl_index)
    disp('No "SiAl" Particles Detected')
else
    unclassed(SiAl_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiAl_index)=0;']); end    
end  


for i=1:n,

    %
    %       SiAlK
    %

    if K(i)/(Al(i)+Si(i))>=0.101
        if K(i)/(Al(i)+Si(i))<=3
            if Al(i)/Si(i)>=0.2
                if Al(i)/Si(i)<=2
                    if (Al(i)+Si(i)+K(i))/total(i)>=0.4
                        if (Al(i)+Si(i)+K(i))/total(i)<=1.1
                            if K(i)/total(i)>=0.0025
                                if K(i)/total(i)<=1.1
                                    if Na(i)/(Al(i)+Si(i)+K(i))<0.05
                                        if Mg(i)/(Al(i)+Si(i)+K(i))<0.08
                                            if P(i)/(Al(i)+Si(i)+K(i))<0.2
                                                if S(i)/(Al(i)+Si(i)+K(i))<0.1
                                                    if Cl(i)/(Al(i)+Si(i)+K(i))<0.1
                                                        if Ca(i)/(Al(i)+Si(i)+K(i))<0.1
                                                            if Ti(i)/(Al(i)+Si(i)+K(i))<0.05
                                                                if Cr(i)/(Al(i)+Si(i)+K(i))<0.05
                                                                    if Fe(i)/(Al(i)+Si(i)+K(i))<0.05
                                                                        SiAlK(i,:)=input(i,:);
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
SiAlK_index=find(SiAlK(:,3)>0);
if isempty(SiAlK_index)
    disp('No "SiAlK" Particles Detected')
else
    unclassed(SiAlK_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiAlK_index)=0;']); end    
end  

for i=1:n,

    %
    %  SiAlNa
    %
    if Na(i)/(Al(i)+Si(i))>=0.101
        if Na(i)/(Al(i)+Si(i))<=3
            if Al(i)/Si(i)>=0.2
                if Al(i)/Si(i)<=2
                    if Ca(i)/Na(i)<0.25
                        if (Al(i)+Si(i)+Na(i))/total(i)>=0.4
                            if (Al(i)+Si(i)+Na(i))/total(i)<=1.1
                                if Mg(i)/(Al(i)+Si(i)+Na(i))<0.15
                                    if P(i)/(Al(i)+Si(i)+Na(i))<0.2
                                        if S(i)/(Al(i)+Si(i)+Na(i))<0.1
                                            if Cl(i)/(Al(i)+Si(i)+Na(i))<0.05
                                                if K(i)/(Al(i)+Si(i)+Na(i))<0.05
                                                    if Ca(i)/(Al(i)+Si(i)+Na(i))<0.05
                                                        if Ti(i)/(Al(i)+Si(i)+Na(i))<0.05
                                                            if Cr(i)/(Al(i)+Si(i)+Na(i))<0.05
                                                                if Fe(i)/(Al(i)+Si(i)+Na(i))<0.15
                                                                    SiAlNa(i,:)=input(i,:);
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
SiAlNa_index=find(SiAlNa(:,3)>0);
if isempty(SiAlNa_index)
    disp('No "SiAlNa" Particles Detected')
else
    unclassed(SiAlNa_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiAlNa_index)=0;']); end   
end  

for i=1:n,

    %
    %       SiAlNaCa
    %
    
    if (Na(i)+Ca(i))/(Al(i)+Si(i))>=0.101
        if (Na(i)+Ca(i))/(Al(i)+Si(i))<=3
            if Al(i)/Si(i)>=0.2
                if Al(i)/Si(i)<=2
                    if Ca(i)/(Al(i)+Si(i))>=0.101
                        if Ca(i)/(Al(i)+Si(i))<=3
                            if Ca(i)/Na(i)>=0.2501
                                if Ca(i)/Na(i)<=5.5
                                    if (Al(i)+Ca(i)+Si(i)+Na(i))/total(i)>=0.4
                                        if (Al(i)+Ca(i)+Si(i)+Na(i))/total(i)<=1.1
                                            if Mg(i)/(Al(i)+Si(i)+Ca(i)+Na(i))<0.1
                                                if P(i)/(Al(i)+Si(i)+Na(i)+Ca(i))<0.2
                                                    if S(i)/(Al(i)+Si(i)+Na(i)+Ca(i))<0.2
                                                        if Cl(i)/(Al(i)+Si(i)+Na(i)+Ca(i))<0.05
                                                            if K(i)/(Al(i)+Si(i)+Na(i)+Ca(i))<0.1
                                                                if Ti(i)/(Al(i)+Si(i)+Na(i)+Ca(i))<0.05
                                                                    if Cr(i)/(Al(i)+Si(i)+Na(i)+Ca(i))<0.05
                                                                        if Fe(i)/(Al(i)+Si(i)+Na(i)+Ca(i))<0.1
                                                                            SiAlNaCa(i,:)=input(i,:);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
SiAlNaCa_index=find(SiAlNaCa(:,3)>0);
if isempty(SiAlNaCa_index)
    disp('No "SiAlNaCa" Particles Detected')
else
    unclassed(SiAlNaCa_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiAlNaCa_index)=0;']); end     
end  


for i=1:n,

    %
    %       SiAlNaK
    %
    
    if (Na(i)+K(i))/(Al(i)+Si(i))>=0.101
        if (Na(i)+K(i))/(Al(i)+Si(i))<=3
            if Al(i)/Si(i)>=0.2
                if Al(i)/Si(i)<=2
                    if K(i)/Na(i)>=0.25
                        if K(i)/Na(i)<=4
                            if (Al(i)+K(i)+Si(i)+Na(i))/total(i)>=0.2501
                                if (Al(i)+K(i)+Si(i)+Na(i))/total(i)<=1.1
                                    if Na(i)/total(i)>=0.05
                                        if Na(i)/total(i)<=1.1
                                            if K(i)/total(i)>=0.05
                                                if K(i)/total(i)<=1.1
                                                    if Mg(i)/(Al(i)+Si(i)+K(i)+Na(i))<0.05
                                                        if P(i)/(Al(i)+Si(i)+Na(i)+K(i))<0.2
                                                            if S(i)/(Al(i)+Si(i)+Na(i)+K(i))<0.2
                                                                if Cl(i)/(Al(i)+Si(i)+Na(i)+K(i))<0.05
                                                                    if Ca(i)/(Al(i)+Si(i)+Na(i)+K(i))<0.1
                                                                        if Ti(i)/(Al(i)+Si(i)+Na(i)+K(i))<0.05
                                                                            if Cr(i)/(Al(i)+Si(i)+Na(i)+K(i))<0.05
                                                                                if Fe(i)/(Al(i)+Si(i)+Na(i)+K(i))<0.05
                                                                                    SiAlNaK(i,:)=input(i,:);
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
SiAlNaK_index=find(SiAlNaK(:,3)>0);
if isempty(SiAlNaK_index)
    disp('No "SiAlNaK" Particles Detected')
else
    unclassed(SiAlNaK_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiAlNaK_index)=0;']); end     
end  


for i=1:n,

    %
    %           SiAlCaFeMg
    %
    
    if (Ca(i)+Fe(i)+Mg(i))/(Al(i)+Si(i))>=0.101
        if (Ca(i)+Fe(i)+Mg(i))/(Al(i)+Si(i))<=3
            if Al(i)/Si(i)>=0.2
                if Al(i)/Si(i)<=2
                    if Ca(i)/(Fe(i)+Mg(i))>=0.25
                        if Ca(i)/(Fe(i)+Mg(i))<=10
                            if (Ca(i)+Fe(i)+Mg(i)+Al(i)+Si(i))/total(i)>=0.4
                                if (Ca(i)+Fe(i)+Mg(i)+Al(i)+Si(i))/total(i)<=1.1
                                    if Ca(i)/total(i)>=0.05
                                        if Ca(i)/total(i)<=1.1
                                            if Fe(i)/total(i)>=0.025
                                                if Fe(i)/total(i)<=1.1
                                                    if Mg(i)/total(i)>=0.025
                                                        if Mg(i)/total(i)<=1.1
                                                            if Ca(i)/(Al(i)+Si(i))<0.5
                                                                if Na(i)/(Si(i)+Al(i)+Ca(i)+Fe(i)+Mg(i))<0.05
                                                                    if P(i)/(Si(i)+Al(i)+Ca(i)+Fe(i)+Mg(i))<0.2
                                                                        if S(i)/(Si(i)+Al(i)+Ca(i)+Fe(i)+Mg(i))<0.2
                                                                            if Cl(i)/(Si(i)+Al(i)+Ca(i)+Fe(i)+Mg(i))<0.1
                                                                                if K(i)/(Si(i)+Al(i)+Ca(i)+Fe(i)+Mg(i))<0.05
                                                                                    if Ti(i)/(Si(i)+Al(i)+Ca(i)+Fe(i)+Mg(i))<0.05
                                                                                        if Cr(i)/(Si(i)+Al(i)+Ca(i)+Fe(i)+Mg(i))<0.05
                                                                                            SiAlCaFeMg(i,:)=input(i,:);
                                                                                        end
                                                                                    end
                                                                                end
                                                                            end
                                                                            
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                            
                                        end
                                    end
                                end
                                
                            end
                        end
                    end
                    
                end
            end
            
        end
                                                                                        
    end

end
SiAlCaFeMg_index=find(SiAlCaFeMg(:,3)>0);
if isempty(SiAlCaFeMg_index)
    disp('No "SiAlCaFeMg" Particles Detected')
else
    unclassed(SiAlCaFeMg_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiAlCaFeMg_index)=0;']); end     
end  


for i=1:n,

    %
    %           SiAlKFeMg
    %
    
    if (K(i)+Fe(i)+Mg(i))/(Al(i)+Si(i))>=0.101
        if (K(i)+Fe(i)+Mg(i))/(Al(i)+Si(i))<=3
            if K(i)/(Al(i)+Si(i))>=0.101
                if K(i)/(Al(i)+Si(i))<=3
                    if (Fe(i)+Mg(i))/(Al(i)+Si(i))>=0.101
                        if (Fe(i)+Mg(i))/(Al(i)+Si(i))<=3
                            if K(i)/(Fe(i)+Mg(i))>=0.25
                                if K(i)/(Fe(i)+Mg(i))<=4
                                    if (K(i)+Fe(i)+Mg(i)+Al(i)+Si(i))/total(i)>=0.4
                                        if (K(i)+Fe(i)+Mg(i)+Al(i)+Si(i))/total(i)<=1.1
                                            if Ca(i)/total(i)<0.05
                                                if Na(i)/(Si(i)+Al(i)+K(i)+Fe(i)+Mg(i))<0.1
                                                    if P(i)/(Si(i)+Al(i)+K(i)+Fe(i)+Mg(i))<0.2
                                                        if S(i)/(Si(i)+Al(i)+K(i)+Fe(i)+Mg(i))<0.2
                                                            if Cl(i)/(Si(i)+Al(i)+K(i)+Fe(i)+Mg(i))<0.1
                                                                if Ca(i)/(Si(i)+Al(i)+K(i)+Fe(i)+Mg(i))<0.05
                                                                    if Ti(i)/(Si(i)+Al(i)+K(i)+Fe(i)+Mg(i))<0.05
                                                                        if Cr(i)/(Si(i)+Al(i)+K(i)+Fe(i)+Mg(i))<0.05
                                                                            SiAlKFeMg(i,:)=input(i,:);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                            
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                                
                            end
                        end
                    end
                    
                end
            end
        end
    end
    
end
SiAlKFeMg_index=find(SiAlKFeMg(:,3)>0);
if isempty(SiAlKFeMg_index)
    disp('No "SiAlKFeMg" Particles Detected')
else
    unclassed(SiAlKFeMg_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiAlKFeMg_index)=0;']); end     
end  


for i=1:n,

    %
    %       SiAlFeMg
    %
    
    if Al(i)/total(i)>=0.1
        if Al(i)/total(i)<=0.8
            if Fe(i)/total(i)>=0.05
                if Fe(i)/total(i)<=0.8
                    if Mg(i)/total(i)>=0.05
                        if Mg(i)/total(i)<=0.8
                            if Ca(i)/total(i)<0.05
                                if (Fe(i)+Mg(i))/(Al(i)+Si(i))>=0.101
                                    if (Fe(i)+Mg(i))/(Al(i)+Si(i))<=3
                                        if Al(i)/Si(i)>=0.201
                                            if Al(i)/Si(i)<=2
                                                if K(i)/(Al(i)+Si(i))<0.1
                                                    if (Al(i)+Si(i)+Fe(i)+Mg(i))/total(i)>=0.5
                                                        if (Al(i)+Si(i)+Fe(i)+Mg(i))/total(i)<=1.1
                                                            if Na(i)/(Si(i)+Al(i)+Fe(i)+Mg(i))<0.05
                                                                if P(i)/(Si(i)+Al(i)+Fe(i)+Mg(i))<0.2
                                                                    if S(i)/(Si(i)+Al(i)+Fe(i)+Mg(i))<0.2
                                                                        if Cl(i)/(Si(i)+Al(i)+Fe(i)+Mg(i))<0.05
                                                                            if K(i)/(Si(i)+Al(i)+Fe(i)+Mg(i))<0.1
                                                                                if Ca(i)/(Si(i)+Al(i)+Fe(i)+Mg(i))<0.1
                                                                                    if Ti(i)/(Si(i)+Al(i)+Fe(i)+Mg(i))<0.05
                                                                                        if Cr(i)/(Si(i)+Al(i)+Fe(i)+Mg(i))<0.05
                                                                                            SiAlFeMg(i,:)=input(i,:);
                                                                                        end
                                                                                    end
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
SiAlFeMg_index=find(SiAlFeMg(:,3)>0);
if isempty(SiAlFeMg_index)
    disp('No "SiAlFeMg" Particles Detected')
else
    unclassed(SiAlFeMg_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiAlFeMg_index)=0;']); end     
end  


for i=1:n,

    %
    %       SiMgFe
    %
    
    if Fe(i)/(Si(i)+Mg(i))>=0.201
        if Fe(i)/(Si(i)+Mg(i))<=10
            if (Mg(i)+Fe(i))/Si(i)>=0.25
                if (Mg(i)+Fe(i))/Si(i)<=4
                    if Al(i)/Si(i)<0.2
                        if (Fe(i)+Mg(i)+Si(i))/total(i)>=0.4
                            if (Fe(i)+Mg(i)+Si(i))/total(i)<=1.1
                                if Na(i)/(Si(i)+Fe(i)+Mg(i))<0.1
                                    if Al(i)/(Si(i)+Fe(i)+Mg(i))<0.05
                                        if P(i)/(Si(i)+Fe(i)+Mg(i))<0.2
                                            if S(i)/(Si(i)+Fe(i)+Mg(i))<0.2
                                                if Cl(i)/(Si(i)+Fe(i)+Mg(i))<0.1
                                                   if K(i)/(Si(i)+Fe(i)+Mg(i))<0.1
                                                       if Ca(i)/(Si(i)+Fe(i)+Mg(i))<0.1
                                                           if Ti(i)/(Si(i)+Fe(i)+Mg(i))<0.05
                                                               if Cr(i)/(Si(i)+Fe(i)+Mg(i))<0.05
                                                                   SiMgFe(i,:)=input(i,:);
                                                               end
                                                           end
                                                       end
                                                   end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
SiMgFe_index=find(SiMgFe(:,3)>0);
if isempty(SiMgFe_index)
    disp('No "SiMgFe" Particles Detected')
else
    unclassed(SiMgFe_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiMgFe_index)=0;']); end     
end  

for i=1:n,

    %
    %           SiMg
    %
    
    
    if Mg(i)/Si(i)>=0.25
        if Mg(i)/Si(i)<=4
            if Al(i)/Si(i)<0.2
                if (Mg(i)+Si(i))/total(i)>=0.4
                    if (Mg(i)+Si(i))/total(i)<=1.1
                        if Na(i)/(Si(i)+Mg(i))<0.1
                            if Al(i)/(Si(i)+Mg(i))<0.1
                                if P(i)/(Si(i)+Mg(i))<0.2
                                    if S(i)/(Si(i)+Mg(i))<0.2
                                        if Cl(i)/(Si(i)+Mg(i))<0.1
                                            if K(i)/(Si(i)+Mg(i))<0.1
                                                if Ca(i)/(Si(i)+Mg(i))<0.1
                                                    if Ti(i)/(Si(i)+Mg(i))<0.05
                                                       if Cr(i)/(Si(i)+Mg(i))<0.05
                                                           if Fe(i)/(Si(i)+Mg(i))<0.2
                                                               SiMg(i,:)=input(i,:);
                                                           end
                                                       end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
SiMg_index=find(SiMg(:,3)>0);
if isempty(SiMg_index)
    disp('No "SiMg" Particles Detected')
else
    unclassed(SiMg_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiMg_index)=0;']); end    
end  


for i=1:n,

    %
    %               SiCaTi
    %
                  
   if Ca(i)/Ti(i)>=0.25
       if Ca(i)/Ti(i)<=4
           if Al(i)/Si(i)<0.2
               if (Si(i)+Ca(i)+Ti(i))/total(i)>=0.4
                   if (Si(i)+Ca(i)+Ti(i))/total(i)<=1.1
                       if Ca(i)/Si(i)>=0.101
                           if Ca(i)/Si(i)<=10
                               if Ti(i)/Si(i)>=0.101
                                   if Ti(i)/Si(i)<=10
                                       if Na(i)/(Si(i)+Ca(i)+Ti(i))<0.1
                                           if Mg(i)/(Si(i)+Ca(i)+Ti(i))<0.1
                                               if P(i)/(Si(i)+Ca(i)+Ti(i))<0.2
                                                   if S(i)/(Si(i)+Ca(i)+Ti(i))<0.2
                                                       if Cl(i)/(Si(i)+Ca(i)+Ti(i))<0.1
                                                           if K(i)/(Si(i)+Ca(i)+Ti(i))<0.1
                                                              if Cr(i)/(Si(i)+Ca(i)+Ti(i))<0.05
                                                                  if Fe(i)/(Si(i)+Ca(i)+Ti(i))<0.2
                                                                      SiCaTi(i,:)=input(i,:);
                                                                  end
                                                              end
                                                           end
                                                       end
                                                   end
                                               end
                                           end
                                       end
                                   end
                               end
                           end
                       end
                   end
               end
           end
       end
   end
  
end
SiCaTi_index=find(SiCaTi(:,3)>0);
if isempty(SiCaTi_index)
    disp('No "SiCaTi" Particles Detected')
else
    unclassed(SiCaTi_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(SiCaTi_index)=0;']); end     
end 

%
%
%       Silicate Particles (Classifications from K7 and others)
%
%

silicateK7=zeros(n,length(txt));

for i=1:n,
    thresh=0.2*total(i);
    if Si(i) > thresh
        if Na(i)/Si(i) < 0.7
            if Mg(i)/Si(i) < 1.33
                if Al(i)/Si(i) < 1.33
                    if K(i)/Si(i) < 0.5
                        if Ca(i)/Si(i) < 0.5
                            if Ti(i)/Si(i) < 0.5
                                if Fe(i)/Si(i) < 0.5
                                    if P(i) + S(i) + Cl(i) < thresh
                                        silicateK7(i,:)=input(i,:);
                                    else
                                        silicateK7(i,:)=0;
                                    end
                                else
                                    silicateK7(i,:)=0;
                                end
                            else
                                silicateK7(i,:)=0;
                            end
                        else
                            silicateK7(i,:)=0;
                        end
                    else
                        silicateK7(i,:)=0;
                    end
                else
                    silicateK7(i,:)=0;
                end
            else
                silicateK7(i,:)=0;
            end
        else
            silicateK7(i,:)=0;
        end
    else
        silicateK7(i,:)=0;
    end
    
    if Si(i) >= 0.2*total(i)         % Si > 20% for silicate classification
        if Si(i) >= 0.6*total(i)         % Si > 60%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    silicateK7(i,:)=input(i,:);
                end
            end
        end
        
        if Al(i)+Si(i) >= 0.6*total(i)        % Al+Si > 60%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    silicateK7(i,:)=input(i,:);
                end
            end
        end
        
        if Si(i)+Fe(i) >= 0.6*total(i)        % Si+Fe > 60%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    silicateK7(i,:)=input(i,:);
                end
            end
        end
        
        if Al(i)+Si(i)+Fe(i) >= 0.5*total(i)        % Al+Si+Fe > 50%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    silicateK7(i,:)=input(i,:);
                end
            end
        end
        
        if Al(i)+Si(i)+Na(i) >= 0.5*total(i)        % Al+Si+Na > 50%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    silicateK7(i,:)=input(i,:);
                end
            end
        end
        
        if Al(i)+Si(i)+Mg(i) >= 0.5*total(i)        % Al+Si+Mg > 50%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    silicateK7(i,:)=input(i,:);
                end
            end
        end
        
        if Al(i)+Si(i)+K(i) >= 0.5*total(i)        % Al+Si+K > 50%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    silicateK7(i,:)=input(i,:);
                end
            end
        end
        
        if Al(i)+Si(i)+Ca(i) >= 0.5*total(i)        % Al+Si+Ca > 50%
            if Cl(i)/Si(i) < 0.2
                if S(i)/Si(i) < 0.2
                    silicateK7(i,:)=input(i,:);
                end
            end
        end
        
        if Al(i)+Si(i)+Ti(i) >= 0.5*total(i)        % Al+Si+Ti > 50%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    silicateK7(i,:)=input(i,:);
                end
            end
        end
        
        
        if Si(i) >= 0.5*total(i)        % Si > 50%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    if Mg(i)>=0.1*total(i)                % Mg > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                    if K(i)>=0.1*total(i)               % K > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                    if Ca(i)>=0.1*total(i)               % Ca > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                end
            end
        end
        
        
        if Al(i)+Si(i) >= 0.5*total(i)        % Si+Al > 50%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    if Mg(i)>=0.1*total(i)                % Mg > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                    if K(i)>=0.1*total(i)               % K > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                    if Ca(i)>=0.1*total(i)               % Ca > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                end
            end
        end
        
        
        if Si(i)+Fe(i) >= 0.5*total(i)        % Si+Fe > 50%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    if Mg(i)>=0.1*total(i)                % Mg > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                    if K(i)>=0.1*total(i)               % K > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                    if Ca(i)>=0.1*total(i)               % Ca > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                end
            end
        end
        
        if Al(i)+Si(i)+Fe(i) >= 0.5*total(i)        % Si+Al+Fe > 50%
            if S(i)/Si(i) < 0.2
                if Cl(i)/Si(i) < 0.2
                    if Mg(i)>=0.1*total(i)                % Mg > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                    if K(i)>=0.1*total(i)               % K > 10%
                        silicateK7(i,:)=input(i,:);
                        
                    end
                    if Ca(i)>=0.1*total(i)               % Ca > 10%
                        silicateK7(i,:)=input(i,:);
                    end
                end
            end
        end
    end
    
end



silicateK7_index=find(silicateK7(:,3)>0);
if isempty(silicateK7_index)
    disp('No Silicate (K7) Particles Detected')
else
    unclassed(silicateK7_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(silicateK7_index)=0;']); end     
end




%
%
%       Aged Silicatess
%
%

agedmd=zeros(n,length(txt));

for i=1:n,

    thresh=0.7*total(i);
    if Na(i) + S(i) + Mg(i) + Al(i) + Si(i) + K(i) + Ca(i) > thresh
        if S(i)/Si(i) <= 2
            if S(i)/Si(i) >= 0.6
                agedmd(i,:)=input(i,:);
            else
                agedmd(i,:)=0;
            end
        else
            agedmd(i,:)=0;
        end
    else
        agedmd(i,:)=0;
    end
    
    if Al(i)+Si(i) >= 0.6*total(i)
        if S(i)/Si(i) > 0.2
            agedmd(i,:)=input(i,:);
        end
        
    end
    
    thresh=0.2*total(i);
    if Si(i) >= thresh
        if Na(i)/Si(i) < 0.7
            if Mg(i)/Si(i) < 1.33
                if Al(i)/Si(i) < 1.33
                    if K(i)/Si(i) < 0.5
                        if Ca(i)/Si(i) < 0.5
                            if Ti(i)/Si(i) < 0.5
                                if Fe(i)/Si(i) < 0.5
                                    if P(i) + Cl(i) < thresh
                                        if S(i) > thresh
                                            agedmd(i,:)=input(i,:);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    if Si(i) >= 0.1*total(i)             % Si > 10% for Mixed Silicate Classification
        if Si(i) + S(i) >= 0.5*total(i)          % Si + S > 50%
            if S(i)/Si(i)>0.2
                agedmd(i,:)=input(i,:);
            end
        end
        
        if Al(i)+Si(i)+S(i) >= 0.5*total(i)        % Si+Al+S > 50%
            if S(i)/Si(i)>0.2
                agedmd(i,:)=input(i,:);
            end
        end
        
        if Si(i)+Fe(i)+S(i) >= 0.5*total(i)       % Si + S + Fe > 50%
            if S(i)/Si(i)>0.2
                agedmd(i,:)=input(i,:);
            end
        end
        
        if Al(i)+Si(i)+Fe(i)+S(i) >= 0.5*total(i)   % Si+Al+Fe+S > 50%
            if S(i)/Si(i)>0.2
                agedmd(i,:)=input(i,:);
            end
        end
        
        
        if Si(i) + S(i) >= 0.4*total(i)          % Si + S > 40%
            if S(i)/Si(i)>0.2
                if Mg(i)>=0.1*total(i)                        % & Mg > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
                if K(i)>=0.1*total(i)                        % & K > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
                if Ca(i)>=0.1*total(i)                        % & Ca > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
            end
        end
        
        
        if Al(i)+Si(i)+S(i) >= 0.4*total(i)   % Si + Al + S > 40%
            if S(i)/Si(i)>0.2
                if Mg(i)>=0.1*total(i)                        % & Mg > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
                if K(i)>=0.1*total(i)                        % & K > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
                if Ca(i)>=0.1*total(i)                        % & Ca > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
            end
        end
        
        
        
        if Si(i)+Fe(i)+S(i) >= 0.4*total(i)  % Si + Fe + S > 40%
            if S(i)/Si(i)>0.2
                if Mg(i)>=0.1*total(i)                        % & Mg > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
                if K(i)>=0.1*total(i)                        % & K > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
                if Ca(i)>=0.1*total(i)                        % & Ca > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
            end
        end
        
        if Al(i)+Si(i)+Fe(i)+S(i) >= 0.4*total(i)   % Si+Al+Fe+S>40%
            if S(i)/Si(i)>0.2
                if Mg(i)>=0.1*total(i)                        % & Mg > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
                if K(i)>=0.1*total(i)                        % & K > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
                if Ca(i)>=0.1*total(i)                        % & Ca > 10%
                    agedmd(i,:)=input(i,:);
                    
                end
            end
        end
    end
%     
%     %
% %
% %       Fe-S Particles
% %
% %

    thresh=0.15*total(i);
    if Fe(i) > thresh
        if Si(i)/Fe(i) < 1
            if Ti(i)/Fe(i) < 1.33
                if S(i)+Fe(i)>0.4*total(i)
                agedmd(i,:)=input(i,:);
                
                end
                
            end
        end
    end 

    
% %
% %
% %       Ti+S Particles
% %
% %

    thresh=0.3*total(i);
    if Ti(i) > thresh
        if Na(i)/Ti(i) < 1
            if Mg(i)/Ti(i) < 1
                if Al(i)/Ti(i) < 1
                    if Si(i)/Ti(i) < 1
                        if Fe(i)/Ti(i) < 1
                            if Ti(i)+S(i)>0.4*total(i)
                                agedmd(i,:)=input(i,:);
                            end
                        end
                    end
                end
            end
        end  
    end

    
    if Ti(i)+S(i)>0.5*total(i);
        agedmd(i,:)=input(i,:);
    end

end

agedmd_index=find(agedmd(:,3)>0);
if isempty(agedmd_index)
    disp('No Aged Silicates Particles Detected')
else
    unclassed(agedmd_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(agedmd_index)=0;']); end 
end

%
%
%               Mixtures
%
%


mixSiS=zeros(n,length(txt));
mixAlSiS=zeros(n,length(txt));
mixNaClSi=zeros(n,length(txt));
mixNaClAlSi=zeros(n,length(txt));
mixCaSi=zeros(n,length(txt));
mixCaAlSi=zeros(n,length(txt));

for i=1:n,
    
    %
    %           MixSiS
    %
    
    if partsize(i)>0.2
        if Al(i)/total(i)<0.05
            if S(i)/total(i)>=0.05
                if S(i)/total(i)<=0.9
                    if S(i)/Si(i)>=0.5
                        if S(i)/Si(i)<=4
                            if Al(i)/Si(i)<0.2
                                if (Si(i)+S(i))/total(i)>=0.3
                                    if (Si(i)+S(i))/total(i)<=1.1
                                        if Na(i)/(Si(i)+S(i))<2
                                            if Mg(i)/(Si(i)+S(i))<2
                                                if Al(i)/(Si(i)+S(i))<0.2
                                                    if P(i)/(Si(i)+S(i))<0.2
                                                        if Cl(i)/(Si(i)+S(i))<0.05
                                                            if K(i)/(Si(i)+S(i))<2
                                                                if Ca(i)/(Si(i)+S(i))<2
                                                                    if Ti(i)/(Si(i)+S(i))<0.2
                                                                        if Cr(i)/(Si(i)+S(i))<0.2
                                                                            if Fe(i)/(Si(i)+S(i))<5
                                                                                mixSiS(i,:)=input(i,:);
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
mixSiS_index=find(mixSiS(:,3)>0);
if isempty(mixSiS_index)
    disp('No "SiS Mixtures" Particles Detected')
else
    unclassed(mixSiS_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(mixSiS_index)=0;']); end     
end 

for i=1:n,

    %
    %       mixAlSiS
    %
    
    
    if Al(i)/total(i)>=0.05
        if Al(i)/total(i)<=0.9
            if Si(i)/total(i)>=0.1
                if Si(i)/total(i)<=0.9
                    if S(i)/total(i)>=0.1
                        if S(i)/total(i)<=0.9
                            if S(i)/Si(i)>=0.5
                                if S(i)/Si(i)<=10
                                    if Al(i)/Si(i)>=0.201
                                        if Al(i)/Si(i)<=5
                                            if (Si(i)+Al(i)+S(i))/total(i)>=0.3
                                                if (Si(i)+Al(i)+S(i))/total(i)<=1.1
                                                    if Na(i)/(Si(i)+Al(i)+S(i))<5
                                                        if Mg(i)/(Si(i)+Al(i)+S(i))<5
                                                            if P(i)/(Si(i)+Al(i)+S(i))<0.2
                                                                if Cl(i)/(Si(i)+Al(i)+S(i))<0.05
                                                                    if K(i)/(Si(i)+Al(i)+S(i))<5
                                                                        if Ca(i)/(Si(i)+Al(i)+S(i))<5
                                                                            if Ti(i)/(Si(i)+Al(i)+S(i))<0.2
                                                                                if Cr(i)/(Si(i)+Al(i)+S(i))<0.2
                                                                                    if Fe(i)/(Si(i)+Al(i)+S(i))<5
                                                                                        mixAlSiS(i,:)=input(i,:);
                                                                                    end
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
mixAlSiS_index=find(mixAlSiS(:,3)>0);
if isempty(mixAlSiS_index)
    disp('No "AlSiS Mixtures" Particles Detected')
else
    unclassed(mixAlSiS_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(mixAlSiS_index)=0;']); end     
end



for i=1:n,

    %
    %       mixNaClSi
    %
    
    if Si(i)/(Na(i)+Cl(i))>=0.5
        if Si(i)/(Na(i)+Cl(i))<=100
            if Al(i)/Si(i)<0.2
                if (Si(i)+Na(i)+Cl(i))/total(i)>=0.2
                    if (Si(i)+Na(i)+Cl(i))/total(i)<=1.1
                        if Cl(i)/total(i)>=0.05
                            if Cl(i)/total(i)<=1.1
                                if Na(i)/total(i)>=0.05
                                    if Na(i)/total(i)<=1.1
                                        if Si(i)/total(i)>=0.01
                                            if Si(i)/total(i)<=1.1
                                                mixNaClSi(i,:)=input(i,:);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
mixNaClSi_index=find(mixNaClSi(:,3)>0);
if isempty(mixNaClSi_index)
    disp('No "NaClSi Mixtures" Particles Detected')
else
    unclassed(mixNaClSi_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(mixNaClSi_index)=0;']); end     
end 

for i=1:n,

    %
    %   MixNaClAlSi
    %
    
    if (Si(i)+Al(i))/(Na(i)+Cl(i))>=0.075
        if (Si(i)+Al(i))/(Na(i)+Cl(i))<=100
            if Al(i)/Si(i)>=0.201
                if Al(i)/Si(i)<=100
                    if (Si(i)+Na(i)+Cl(i))/total(i)>=0.2
                        if (Si(i)+Na(i)+Cl(i))/total(i)<=1.1
                            if Cl(i)/total(i)>=0.05
                                if Cl(i)/total(i)<=1.1
                                    if Na(i)/total(i)>=0.05
                                        if Na(i)/total(i)<=1.1
                                            if Si(i)/total(i)>=0.025
                                                if Si(i)/total(i)<=1.1
                                                   if Al(i)/total(i)>=0.01
                                                       if Al(i)/total(i)<=1.1
                                                           mixNaClAlSi(i,:)=input(i,:);
                                                       end
                                                   end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
mixNaClAlSi_index=find(mixNaClAlSi(:,3)>0);
if isempty(mixNaClAlSi_index)
    disp('No "NaClAlSi Mixtures" Particles Detected')
else
    unclassed(mixNaClAlSi_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(mixNaClAlSi_index)=0;']); end     
end 


for i=1:n,

    %
    %           MixCaSi
    %
    
    
    if Si(i)/Ca(i)>=0.1001
        if Si(i)/Ca(i)<=100
            if Al(i)/Si(i)<0.2
                if (Ca(i)+Si(i))/total(i)>=0.2
                    if (Ca(i)+Si(i))/total(i)<=1.1
                        if Si(i)/total(i)>=0.01
                            if Si(i)/total(i)<=1.1
                                if Ca(i)/total(i)>=0.05
                                    if Ca(i)/total(i)<=1.1
                                        if Mg(i)/(Si(i)+Ca(i))<0.1
                                            if Al(i)/(Si(i)+Ca(i))<0.2
                                                if P(i)/(Si(i)+Ca(i))<0.2
                                                    if S(i)/(Si(i)+Ca(i))<0.2
                                                        mixCaSi(i,:)=input(i,:);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
mixCaSi_index=find(mixCaSi(:,3)>0);
if isempty(mixCaSi_index)
    disp('No "CaSi Mixtures" Particles Detected')
else
    unclassed(mixCaSi_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(mixCaSi_index)=0;']); end     
end 

for i=1:n,
    
    %
    %           MixCaAlSi
    %
    
    
    if Al(i)/Si(i)>=0.201
        if Al(i)/Si(i)<=20
            if (Ca(i)+Si(i)+Al(i))/total(i)>=0.2
                if (Ca(i)+Si(i)+Al(i))/total(i)<=1.1
                    if Al(i)/total(i)>=0.01
                        if Al(i)/total(i)<=1.1
                            if Si(i)/total(i)>=0.01
                                if Si(i)/total(i)<=1.1
                                    if Ca(i)/total(i)>=0.05
                                        if Ca(i)/total(i)<=1.1
                                            if Si(i)/Ca(i)>=0.1001
                                                if Si(i)/Ca(i)<=100
                                                    if Na(i)/(Si(i)+Ca(i))<0.2
                                                        if Mg(i)/(Si(i)+Ca(i))<2
                                                            if P(i)/(Si(i)+Ca(i))<0.2
                                                                if S(i)/(Si(i)+Ca(i))<0.2
                                                                    if Cl(i)/(Si(i)+Ca(i))<0.05
                                                                        if K(i)/(Si(i)+Ca(i))<1
                                                                            mixCaAlSi(i,:)=input(i,:);
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end                                                                                                                 
end
mixCaAlSi_index=find(mixCaAlSi(:,3)>0);
if isempty(mixCaAlSi_index)
    disp('No "CaAlSi Mixtures" Particles Detected')
else
    unclassed(mixCaAlSi_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(mixCaAlSi_index)=0;']); end     
end 



%
%
%       Fresh Chlorides with minor Ca
%
%

ssca=zeros(n,length(txt));

for i=1:n,

    if Na(i)+Cl(i)+Ca(i) >= 0.5*total(i)    % Na+Cl+Ca >50%
        if Na(i)/Cl(i)>=0.2
            if Na(i)/Cl(i)<=1.1
                if Si(i)/Cl(i)<0.2
                    if S(i)/Cl(i)<0.2
                        ssca(i,:)=input(i,:);
                    end
                end
            end
        end
    end
end
ssca_index=find(ssca(:,3)>0);
if isempty(ssca_index)
    disp('No Chlorides-Ca Mixed Particles Detected')
else
    unclassed(ssca_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(ssca_index)=0;']); end    
end


%
%
%       Aged Chlorides with minor Ca
%
%

agssca=zeros(n,length(txt));

for i=1:n,

    if Na(i)+Cl(i)+Ca(i)+S(i) >= 0.5*total(i)    % Na+Cl+Ca+S >50%
        if Na(i)/Cl(i)>=0.1
            if Na(i)/Cl(i)<=1.1
                if Si(i)/Cl(i)<0.2
                    if S(i)/Cl(i)>0.2
                        agssca(i,:)=input(i,:);
                    end
                end
            end
        end
    end
end

agssca_index=find(agssca(:,3)>0);
if isempty(agssca_index)
    disp('No Aged Chlorides-Ca Mixed Particles Detected')
else
    unclassed(agssca_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(agssca_index)=0;']); end     
end

%
%
%       Mixed Chlorides (low Si!)
%
%

mixcl=zeros(n,length(txt));

for i=1:n,

    if Cl(i)/total(i)>=0.1
        if Cl(i)/total(i)<=1.1
            if Si(i)/Cl(i)<0.1
                if S(i)/Cl(i)>0.2
                    if Cr(i)/Cl(i)<0.1
                    mixcl(i,:)=input(i,:);
                    end
                end
            end
        end
        
    end
end
mixcl_index=find(mixcl(:,3)>0);
if isempty(mixcl_index)
    disp('No Aged Chlorides-Ca Mixed Particles Detected')
else
    unclassed(mixcl_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(mixcl_index)=0;']); end     
end

%
%
%       Metallic
%               e.g. Metal-based fly-ashes or Metal Oxides
%                       Based on Fe, minor Ni, Cr, Zn, Cu
%

industry=zeros(n,length(txt));
NiRich=zeros(n,length(txt));
CrRich=zeros(n,length(txt));
ZnRich=zeros(n,length(txt));
CuRich=zeros(n,length(txt));

for i=1:n,

    if unclassed(:,3)>0
        if Fe(i)+Ni(i)+Cr(i)+Cu(i)+Zn(i)>= 0.5*total(i)  % Fe+Ni+Cr+Cu+Zn>50%
            %         if Si(i)/(Fe(i)+Ni(i)+Cu(i)+Zn(i))<0.05
            industry(i,:)=input(i,:);
            %         end
        end
    end
end
industry_index=find(industry(:,3)>0);

if isempty(industry_index)
    disp('No Potential Industry Tracer Particles Detected')
else
    unclassed(industry_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(industry_index)=0;']); end 
end
for i=1:n,

    if Zn(i)/total(i)>=0.2               % Zn rich
        if Zn(i)/total(i)<=1.1
            %                 if Si(i)/Zn(i)<0.2
            ZnRich(i,:)=input(i,:);
            %                 end
        end
    end
end
ZnRich_index=find(ZnRich(:,3)>0);
if isempty(ZnRich_index)
    disp('No Zn-Rich Particles Detected')
else
    unclassed(ZnRich_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(ZnRich_index)=0;']); end     
end
for i=1:n,

    if Ni(i)/total(i)>=0.2                      % Ni Rich
        if Ni(i)/total(i)<=1.1
            %             if Si(i)/Ni(i)<0.2
            NiRich(i,:)=input(i,:);
            %             end
        end
    end
end
NiRich_index=find(NiRich(:,3)>0);

if isempty(NiRich_index)
    disp('No Ni-Rich Particles Detected')
else
    unclassed(NiRich_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(NiRich_index)=0;']); end     
end
for i=1:n,

    if Cu(i)/total(i)>=0.2                      % Cu Rich
        if Cu(i)/total(i)>=1.1
            %             if Si(i)/Cu(i)<0.2
            CuRich(i,:)=input(i,:);
            %             end
        end
    end
end
CuRich_index=find(CuRich(:,3)>0);
if isempty(CuRich_index)
    disp('No Cu-Rich Particles Detected')
else
    unclassed(CuRich_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(CuRich_index)=0;']); end     
end
for i=1:n,

    if Cr(i)/total(i)>=0.2                      % Cr Rich
        if Cr(i)/total(i)>=1.1
            %             if Si(i)/Cr(i)<0.2
            CrRich(i,:)=input(i,:);
            %             end
        end
    end
end
CrRich_index=find(CrRich(:,3)>0);
if isempty(CrRich_index)
    disp('No Cr-Rich Particles Detected')
else
    unclassed(CrRich_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(CrRich_index)=0;']); end     
end



%
%
%           Leftover particles
%
%

Cadom=zeros(n,length(txt));
Aldom=zeros(n,length(txt));
AlSidom=zeros(n,length(txt));
Cldom=zeros(n,length(txt));
Pdom=zeros(n,length(txt));
Cudom=zeros(n,length(txt));
Sidom=zeros(n,length(txt));
Kdom=zeros(n,length(txt));
Mgdom=zeros(n,length(txt));

for i=1:n,
    tt=sum(unclassed(i,7:final_entry));
    if Si(i)>=0.1
    if Mg(i)/total(i)>=0.35
        if Mg(i)/total(i)<=1.1
            Mgdom(i,:)=input(i,:);
        else Mgdom(i,:)=0;
        end
    else Mgdom(i,:)=0;
    end
    end
end 
Mgdom_index=find(Mgdom(:,3)>0);
if isempty(Mgdom_index)
    disp('No Mg rich Particles Detected')
else
    unclassed(Mgdom_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Mgdom_index)=0;']); end     
end

for i=1:n,
    if P(i)/total(i)>=0.1
        if P(i)/total(i)<=1.1;
            if P(i)>Al(i)
                if P(i)>Si(i)
                    if P(i)>S(i)
                        if P(i)>Cl(i)
                            if P(i)>K(i)
                                if P(i)>Ca(i)
                                    if P(i)>Ti(i)
                                        if P(i)>Cr(i)
                                            if P(i)>Fe(i)
                                                if P(i)>Ni(i)
                                                    if P(i)>Cu(i)
                                                        if P(i)>Zn(i)
                                                            Pdom(i,:)=input(i,:);
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
            
end
Pdom_index=find(Pdom(:,3)>0);
if isempty(Pdom_index)
    disp('No P-dominant Particles Detected')
else
    unclassed(Pdom_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Pdom_index)=0;']); end     
end

for i=1:n,
    tt=sum(unclassed(i,7:final_entry));
    if Si(i)/total(i)>=0.1
        if Si(i)/total(i)<=1.1
            Sidom(i,:)=input(i,:);
        end
    end
end       
Sidom_index=find(Sidom(:,3)>0);
if isempty(Sidom_index)
    disp('No "Other Si-rich" Particles Detected')
else
    unclassed(Sidom_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Sidom_index)=0;']); end     
end 

for i=1:n,
    if Si(i)>=0.1
        if Al(i)/total(i)>=0.1
            if Al(i)/total(i)<=1.1
                Aldom(i,:)=input(i,:);
            end
        end
    end
end
Aldom_index=find(Aldom(:,3)>0);
if isempty(Aldom_index)
    disp('No Al-rich Particles Detected')
else
    unclassed(Aldom_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Aldom_index)=0;']); end     
end
for i=1:n,
    if Si(i)>=0.1
        if (Al(i)+Si(i))/total(i)>=0.2
            if (Al(i)+Si(i))/total(i)<=1.1
                AlSidom(i,:)=input(i,:);
            end
        end
    end
end
AlSidom_index=find(AlSidom(:,3)>0);
if isempty(AlSidom_index)
    disp('No AlSi-rich Particles Detected')
else
    unclassed(AlSidom_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(AlSidom_index)=0;']); end     
end

for i=1:n,
    tt=sum(unclassed(i,7:final_entry));
    if Cl(i)/total(i)>=0.1
        if Cl(i)/total(i)>=1.1;
                    Cldom(i,:)=input(i,:);
                end
            end
end
Cldom_index=find(Cldom(:,3)>0);
if isempty(Cldom_index)
    disp('No Cl-dominant Particles Detected')
else
    unclassed(Cldom_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Cldom_index)=0;']); end 
end

for i=1:n,
    tt=sum(unclassed(i,7:final_entry));
    if K(i)/total(i)>=0.25
        if K(i)/total(i)<=1.1
            Kdom(i,:)=input(i,:);
        end
    end
end
Kdom_index=find(Kdom(:,3)>0);
if isempty(Kdom_index)
    disp('No "K-dom" Particles Detected')
else
    unclassed(Kdom_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Kdom_index)=0;']); end    
end

for i=1:n,
    if Ca(i)/total(i)>=0.2
        if Ca(i)/total(i)<=1.1
            Cadom(i,:)=input(i,:);
        else Cadom(i,:)=0;
        end
    else Cadom(i,:)=0;
    end  
end
Cadom_index=find(Cadom(:,3)>0);
if isempty(Cadom_index)
    disp('No Ca-dominant Particles Detected')
else
    unclassed(Cadom_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Cadom_index)=0;']); end     
end

for i=1:n,
    tt=sum(unclassed(i,7:final_entry));
    if Cu(i)>=0.2*tt
        if Cu(i)>unclassed(i,7:find(strcmp('Cu', txt)==1))
            if Cu(i)>Zn(i)
                Cudom(i,:)=input(i,:);
            end
        end
    end
end
Cudom_index=find(Cudom(:,3)>0);
if isempty(Cudom_index)
    disp('No Cu-dominant Particles Detected')
else
    unclassed(Cudom_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(Cudom_index)=0;']); end     
end


%
%
%       Other Particles
%
%


other_index=find(unclassed(:,3)>0);

if isempty(other_index)
    disp('No Mixed Particles Detected')
else
    unclassed(other_index,:)=0;
    for j=1:length(elms), eval([char(elms(j)),'(other_index)=0;']); end     
end

carb_index=[ammsulph_index;Narich_index;organic_index];
carbon_index=unique(carb_index);
carbon=input(carbon_index,:);

silicate_index=[SiAl_index;SiAlK_index;SiAlNa_index;SiAlNaCa_index;...
    SiAlNaK_index;SiAlCaFeMg_index;SiAlKFeMg_index;SiAlFeMg_index;...
    SiMgFe_index;SiMg_index;SiCaTi_index;quartz_index;Sidom_index];
silicate=input(silicate_index,:);

mindust_index=[silicateK7_index;silicate_index;AlSidom_index;...
    Mgdom_index];
mineral_dust_index=unique(mindust_index);
mineral_dust=input(mineral_dust_index,:);

Ca_total=[carbonate_index;Ca_index;Cadom_index;Camgcarb_index;Cacarb_index];
Ca_rich_index=unique(Ca_total);
Ca_rich=input(Ca_rich_index,:);

ss=[freshss_index;Cldom_index;ssca_index;NaCl_index;KCl_index;oCl_index];
seasalt_index=unique(ss);

agss=[agedss_index;agssca_index;mixClS_index;mixcl_index];
agseasalt_index=unique(agss);

sul_index=[Srich_index;NaS_index;sulph_index;CaNaS_index];
sulphate_index=unique(sul_index);
sulphate=input(sulphate_index,:);

ox_index=[FeO_index;TiO_index;FeTiO_index;AlO_index;industry_index;Aldom_index;Cudom_index;...
    Ti_index;Fe_index;CrRich_index;CuRich_index;ZnRich_index;NiRich_index];
oxide_index=unique(ox_index);
oxide=input(oxide_index,:);

m_index=[mixSiS_index;mixAlSiS_index;mixNaClSi_index;mixNaClAlSi_index;...
    mixCaSi_index;mixCaAlSi_index;agedmd_index];
mix_index=unique(m_index);
mix=input(mix_index,:);

p_index=[phos_index;Pdom_index];
phosphate_index=unique(p_index);
phosphate=input(phosphate_index,:);

% clear carbon
% carbon=input(carbon_index,:);
clear bio
bio=input(bio_index,:);
clear freshss
freshss=input(seasalt_index,:);
clear agedss
agedss=input(agseasalt_index,:);
clear sulph
sulph=input(sulph_index,:);
clear other
other=input(other_index,:);
clear Narich
Narich=input(Narich_index,:);
clear ammsulph
ammsulph=input(ammsulph_index,:);
clear phos
phos=input(phos_index,:);
clear carbonate
carbonate=input(carbonate_index,:);
clear gypsum
gypsum=input(gypsum_index,:);
clear NaCl
NaCl=input(NaCl_index,:);
clear KCl
KCl=input(KCl_index,:);
clear oCl
oCl=input(oCl_index,:);
clear Ca
Ca=input(Ca_index,:);
clear FeO
FeO=input(FeO_index,:);
clear TiO
TiO=input(TiO_index,:);
clear FeTiO
FeTiO=input(FeTiO_index,:);
clear AlO
AlO=input(AlO_index,:);
clear quartz
quartz=input(quartz_index,:);
clear SiAl
SiAl=input(SiAl_index,:);
clear SiAlK
SiAlK=input(SiAlK_index,:);
clear SiAlNa
SiAlNa=input(SiAlNa_index,:);
clear SiAlNaCa
SiAlNaCa=input(SiAlNaCa_index,:);
clear SiAlNaK
SiAlNaK=input(SiAlNaK_index,:);
clear SiAlCaFeMg
SiAlCaFeMg=input(SiAlCaFeMg_index,:);
clear SiAlKFeMg
SiAlKFeMg=input(SiAlKFeMg_index,:);
clear SiAlFeMg
SiAlFeMg=input(SiAlFeMg_index,:);
clear SiMgFe
SiMgFe=input(SiMgFe_index,:);
clear SiMg
SiMg=input(SiMg_index,:);
clear SiCaTi
SiCaTi=input(SiCaTi_index,:);
clear mixSiS
mixSiS=input(mixSiS_index,:);
clear mixAlSiS
mixAlSiS=input(mixAlSiS_index,:);
clear mixClS
mixClS=input(mixClS_index,:);
clear mixNaClSi
mixNaClSi=input(mixNaClSi_index,:);
clear mixNaClAlSi
mixNaClAlSi=input(mixNaClAlSi_index,:);
clear mixCaSi
mixCaSi=input(mixCaSi_index,:);
clear mixCaAlSi
mixCaAlSi=input(mixCaAlSi_index,:);
clear Fe
Fe=input(Fe_index,:);
clear Ti
Ti=input(Ti_index,:);
clear silicateK7
silicateK7=input(silicateK7_index,:);
clear agedmd
agedmd=input(agedmd_index,:);
clear industry
industry=input(industry_index,:);
clear Mgdom
Mgdom=input(Mgdom_index,:);
clear Pdom
Pdom=input(Pdom_index,:);
clear Cadom
Cadom=input(Cadom_index,:);
clear Aldom
Aldom=input(Aldom_index,:);
clear AlSidom
AlSidom=input(AlSidom_index,:);
clear Cldom
Cldom=input(Cldom_index,:);
clear Cudom
Cudom=input(Cudom_index,:);
clear Sidom
Sidom=input(Sidom_index,:);
clear Kdom
Kdom=input(Kdom_index,:);

ncarbon=length(carbon(:,1));
nbio=length(bio(:,1));
nfreshss=length(freshss(:,1));
nagedss=length(agedss(:,1));
nsulphate=length(sulphate(:,1));
nCa=length(Ca_rich(:,1));
nmindust=length(mineral_dust(:,1));
ngypsum=length(gypsum(:,1));
nmix=length(mix(:,1));
nphos=length(phosphate(:,1));
nox=length(oxide(:,1));
nother=length(other(:,1));
nbiomass=length(Kdom(:,1));

sizelt05=find(input(:,3)<=0.5);
sizegt05=find(input(:,3)>0.5);
lt05=length(sizelt05);
gt05=length(sizegt05);
disp(['Number of particles less than 0.5micron: ',num2str(lt05)])
disp(['Number of particles greater than 0.5micron: ',num2str(gt05)])

%
%
%
%

clsd=[nmindust;nmix;nCa;nfreshss;nagedss;nsulphate;ngypsum;...
    ncarbon;nbio;nphos;nox;nbiomass;nother];

classes=['Silicates      ','Mixed Silicates','Ca Rich           ','Fresh Chlorides    ','Aged Chlorides     ',...
    'Sulphate          ','Gypsum            ','Carbonaceous      ','Biogenic          ',...
    'Phosphate         ','Metallic          ','Biomass Tracers   ','Other             '];


OUTPUT.All.Raw=input;
for i=1:length(OUTPUT.All.Raw(:,1)),
    OUTPUT.All.Normalised(i,:)=OUTPUT.All.Raw(i,7:final_entry)./sum(OUTPUT.All.Raw(i,7:final_entry));
end
OUTPUT.All.Mean=nanmean(OUTPUT.All.Normalised,1);
OUTPUT.All.Number=length(OUTPUT.All.Raw(:,1));
OUTPUT.Silicate.All_Indices=mineral_dust_index;
OUTPUT.Silicate.Raw=mineral_dust;
for i=1:length(OUTPUT.Silicate.All_Indices),
    OUTPUT.Silicate.Normalised(i,:)=OUTPUT.Silicate.Raw(i,7:final_entry)./sum(OUTPUT.Silicate.Raw(i,7:final_entry));
end
if length(OUTPUT.Silicate.All_Indices)>0, OUTPUT.Silicate.Mean=nanmean(OUTPUT.Silicate.Normalised,1); end
OUTPUT.Silicate.Number=length(OUTPUT.Silicate.Raw(:,1));
OUTPUT.MixedSil.All_Indices=mix_index;
OUTPUT.MixedSil.Raw=mix;
for i=1:length(OUTPUT.MixedSil.All_Indices)
    OUTPUT.MixedSil.Normalised(i,:)=OUTPUT.MixedSil.Raw(i,7:final_entry)./sum(OUTPUT.MixedSil.Raw(i,7:final_entry));
end
if length(OUTPUT.MixedSil.All_Indices)>0, OUTPUT.MixedSil.Mean=nanmean(OUTPUT.MixedSil.Normalised,1); end
OUTPUT.MixedSil.Number=length(OUTPUT.MixedSil.Raw(:,1));
OUTPUT.CaRich.All_Indices=Ca_rich_index;
OUTPUT.CaRich.Raw=Ca_rich;
for i=1:length(OUTPUT.CaRich.All_Indices)
    OUTPUT.CaRich.Normalised(i,:)=OUTPUT.CaRich.Raw(i,7:final_entry)./sum(OUTPUT.CaRich.Raw(i,7:final_entry));
end
if length(OUTPUT.CaRich.All_Indices)>0, OUTPUT.CaRich.Mean=nanmean(OUTPUT.CaRich.Normalised,1); end
OUTPUT.CaRich.Number=length(OUTPUT.CaRich.Raw(:,1));
OUTPUT.FreshCl.All_Indices=seasalt_index;
OUTPUT.FreshCl.Raw=freshss;
for i=1:length(OUTPUT.FreshCl.All_Indices)
    OUTPUT.FreshCl.Normalised(i,:)=OUTPUT.FreshCl.Raw(i,7:final_entry)./sum(OUTPUT.FreshCl.Raw(i,7:final_entry));
end
if length(OUTPUT.FreshCl.All_Indices)>0, OUTPUT.FreshCl.Mean=nanmean(OUTPUT.FreshCl.Normalised,1);end
OUTPUT.FreshCl.Number=length(OUTPUT.FreshCl.Raw(:,1));
OUTPUT.AgedCl.All_Indices=agseasalt_index;
OUTPUT.AgedCl.Raw=agedss;
for i=1:length(OUTPUT.AgedCl.All_Indices)
    OUTPUT.AgedCl.Normalised(i,:)=OUTPUT.AgedCl.Raw(i,7:final_entry)./sum(OUTPUT.AgedCl.Raw(i,7:final_entry));
end
if length(OUTPUT.AgedCl.All_Indices)>0, OUTPUT.AgedCl.Mean=nanmean(OUTPUT.AgedCl.Normalised,1); end
OUTPUT.AgedCl.Number=length(OUTPUT.AgedCl.Raw(:,1));
OUTPUT.Sulph.All_Indices=sulphate_index;
OUTPUT.Sulph.Raw=sulphate;
for i=1:length(OUTPUT.Sulph.All_Indices)
    OUTPUT.Sulph.Normalised(i,:)=OUTPUT.Sulph.Raw(i,7:final_entry)./sum(OUTPUT.Sulph.Raw(i,7:final_entry));
end
if length(OUTPUT.Sulph.All_Indices)>0, OUTPUT.Sulph.Mean=nanmean(OUTPUT.Sulph.Normalised,1); end
OUTPUT.Sulph.Number=length(OUTPUT.Sulph.Raw(:,1));
OUTPUT.Gypsum.All_Indices=gypsum_index;
OUTPUT.Gypsum.Raw=gypsum;
for i=1:length(OUTPUT.Gypsum.All_Indices)
    OUTPUT.Gypsum.Normalised(i,:)=OUTPUT.Gypsum.Raw(i,7:final_entry)./sum(OUTPUT.Gypsum.Raw(i,7:final_entry));
end
if length(OUTPUT.Gypsum.All_Indices)>0, OUTPUT.Gypsum.Mean=nanmean(OUTPUT.Gypsum.Normalised,1); end
OUTPUT.Gypsum.Number=length(OUTPUT.Gypsum.Raw(:,1));
OUTPUT.Carbon.All_Indices=carbon_index;
OUTPUT.Carbon.Raw=carbon;
for i=1:length(OUTPUT.Carbon.All_Indices)
    OUTPUT.Carbon.Normalised(i,:)=OUTPUT.Carbon.Raw(i,7:final_entry)./sum(OUTPUT.Carbon.Raw(i,7:final_entry));
end
if length(OUTPUT.Carbon.All_Indices)>0, OUTPUT.Carbon.Mean=nanmean(OUTPUT.Carbon.Normalised,1); end
OUTPUT.Carbon.Number=length(OUTPUT.Carbon.Raw(:,1));
OUTPUT.Bio.All_Indices=bio_index;
OUTPUT.Bio.Raw=bio;
for i=1:length(OUTPUT.Bio.All_Indices)
    OUTPUT.Bio.Normalised(i,:)=OUTPUT.Bio.Raw(i,7:final_entry)./sum(OUTPUT.Bio.Raw(i,7:final_entry));
end
if length(OUTPUT.Bio.All_Indices)>0, OUTPUT.Bio.Mean=nanmean(OUTPUT.Bio.Normalised,1); end
OUTPUT.Bio.Number=length(OUTPUT.Bio.Raw(:,1));
OUTPUT.Metal.All_Indices=oxide_index;
OUTPUT.Metal.Raw=oxide;
for i=1:length(OUTPUT.Metal.All_Indices)
   OUTPUT.Metal.Normalised(i,:)=OUTPUT.Metal.Raw(i,7:final_entry)./sum(OUTPUT.Metal.Raw(i,7:final_entry));
end
if length(OUTPUT.Metal.All_Indices)>0, OUTPUT.Metal.Mean=nanmean(OUTPUT.Metal.Normalised,1); end
OUTPUT.Metal.Number=length(OUTPUT.Metal.Raw(:,1));
OUTPUT.Phos.All_Indices=phosphate_index;
OUTPUT.Phos.Raw=phosphate;
for i=1:length(OUTPUT.Phos.All_Indices)
    OUTPUT.Phos.Normalised(i,:)=OUTPUT.Phos.Raw(i,7:final_entry)./sum(OUTPUT.Phos.Raw(i,7:final_entry));
end
if length(OUTPUT.Phos.All_Indices)>0, OUTPUT.Phos.Mean=nanmean(OUTPUT.Phos.Normalised,1); end
OUTPUT.Phos.Number=length(OUTPUT.Phos.Raw(:,1));
OUTPUT.Biomass.All_Indices=Kdom_index;
OUTPUT.Biomass.Raw=Kdom;
if OUTPUT.Biomass.All_Indices>0
    for i=1:length(OUTPUT.Biomass.All_Indices)
        OUTPUT.Biomass.Normalised(i,:)=OUTPUT.Biomass.Raw(i,7:final_entry)./sum(OUTPUT.Biomass.Raw(i,7:final_entry));
    end
    OUTPUT.Biomass.Mean=nanmean(OUTPUT.Biomass.Normalised,1);
else OUTPUT.Biomass.Mean=zeros(length(text),1)';
end
OUTPUT.Biomass.Number=length(OUTPUT.Biomass.Raw(:,1));
OUTPUT.Other.All_Indices=other_index;
OUTPUT.Other.Raw=other;
for i=1:length(OUTPUT.Other.All_Indices)
    OUTPUT.Other.Normalised(i,:)=OUTPUT.Other.Raw(i,7:final_entry)./sum(OUTPUT.Other.Raw(i,7:final_entry));
end
if length(OUTPUT.Other.All_Indices)>0, OUTPUT.Other.Mean=nanmean(OUTPUT.Other.Normalised,1); end
OUTPUT.Other.Number=length(OUTPUT.Other.Raw(:,1));

end
