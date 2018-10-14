%% Dissimilarity Measures Module (DMM) Demo
%    Copyright (C) 2012  Juan M. Banda, Rafal A. Angryk from Montana State University
%    Contact: juan@jmbanda.com
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% MORE INFO
%  Stand-alone script to ilustrate the usage of the Feature Extraction
%  Module. Remember to set your path to the correct place where the demo
%  subset of the dataset has been extracted
%  For more details
%  Juan M. Banda's dissertation:
%  "FRAMEWORK FOR CREATING LARGE-SCALE CONTENT-BASED IMAGE RETRIEVAL SYSTEM
%  (CBIR) FOR SOLAR DATA ANALYSIS"
%  http://www.jmbanda.com/Dissertation/
%
%  Notes on this DEMO:
%  http://www.jmbanda.com/Framework/Demo/
%
%  Forced Outputs: Saved variables for the KLD, JSD, CHI2 Measures
%  Configurable Outputs: 2D/3D MDS plots, Curve Fitting Plots, Components
%  Plots, Weka files for Hard-threshold and tanget threshold experiments

clearvars;
clases=8;                           %Number of classes
Csize=20;                            %Number of images in each folder
total_im=clases*Csize;
segments=8;                         %Number of Columns/rows  NbyN
class_n=20;                         %class size
sample_size=clases*class_n;         %TOTAL number of images
cells=segments*segments;            %Grid Cell Size
parameters=10;                      %Number of Image Parameters
distances=12;
index=segments*segments*parameters; %Number of elements in the Feature Vector
dataSet='Solar';                    %Dataset Being Manipulated
tang_thres=135;                     %tangent thresholding for components
component_threshold=10;
dncNM={'euclidean' 'seuclidean' 'mahalanobis' 'cityblock' 'cosine' 'correlation' 'spearman' 'chebychev' 'hausdorff' 'KLD' 'CHI2' 'JSD'};
labels = {'\color{red} 1','\color{red} 2','\color{red} 3','\color{red} 4','\color{red} 5','\color{red} 6','\color{red} 7','\color{red} 8','\color{red} 9','\color{red} 10','\color{red} 11','\color{red} 12','\color{red} 13','\color{red} 14','\color{red} 15','\color{red} 16','\color{red} 17','\color{red} 18','\color{red} 19','\color{red} 20','\color{green} 21','\color{green} 22','\color{green} 23','\color{green} 24','\color{green} 25','\color{green} 26','\color{green} 27','\color{green} 28','\color{green} 29','\color{green} 30','\color{green} 31','\color{green} 32','\color{green} 33','\color{green} 34','\color{green} 35','\color{green} 36','\color{green} 37','\color{green} 38','\color{green} 39','\color{green} 40','\color{blue} 41','\color{blue} 42','\color{blue} 43','\color{blue} 44','\color{blue} 45','\color{blue} 46','\color{blue} 47','\color{blue} 48','\color{blue} 49','\color{blue} 50','\color{blue} 51','\color{blue} 52','\color{blue} 53','\color{blue} 54','\color{blue} 55','\color{blue} 56','\color{blue} 57','\color{blue} 58','\color{blue} 59','\color{blue} 60','\color{yellow} 61','\color{yellow} 62','\color{yellow} 63','\color{yellow} 64','\color{yellow} 65','\color{yellow} 66','\color{yellow} 67','\color{yellow} 68','\color{yellow} 69','\color{yellow} 70','\color{yellow} 71','\color{yellow} 72','\color{yellow} 73','\color{yellow} 74','\color{yellow} 75','\color{yellow} 76','\color{yellow} 77','\color{yellow} 78','\color{yellow} 79','\color{yellow} 80','\color{magenta} 81','\color{magenta} 82','\color{magenta} 83','\color{magenta} 84','\color{magenta} 85','\color{magenta} 86','\color{magenta} 87','\color{magenta} 88','\color{magenta} 89','\color{magenta} 90','\color{magenta} 91','\color{magenta} 92','\color{magenta} 93','\color{magenta} 94','\color{magenta} 95','\color{magenta} 96','\color{magenta} 97','\color{magenta} 98','\color{magenta} 99','\color{magenta} 100','\color{gray} 101','\color{gray} 102','\color{gray} 103','\color{gray} 104','\color{gray} 105','\color{gray} 106','\color{gray} 107','\color{gray} 108','\color{gray} 109','\color{gray} 110','\color{gray} 111','\color{gray} 112','\color{gray} 113','\color{gray} 114','\color{gray} 115','\color{gray} 116','\color{gray} 117','\color{gray} 118','\color{gray} 119','\color{gray} 120','\color{orange} 121','\color{orange} 122','\color{orange} 123','\color{orange} 124','\color{orange} 125','\color{orange} 126','\color{orange} 127','\color{orange} 128','\color{orange} 129','\color{orange} 130','\color{orange} 131','\color{orange} 132','\color{orange} 133','\color{orange} 134','\color{orange} 135','\color{orange} 136','\color{orange} 137','\color{orange} 138','\color{orange} 139','\color{orange} 140',' 141',' 142',' 143',' 144',' 145',' 146',' 147',' 148',' 149',' 150',' 151',' 152',' 153',' 154',' 155',' 156',' 157',' 158',' 159',' 160'};
exten='tif';
pathh='C:\Demo\';                   %Root folder of datasets
%Save to path
pathhTOS='C:\Demo\';                   
clasesNames={'Active Region' 'Coronal Jet' 'Emerging Flux' 'Filament' 'Filament Activation' 'Filament Eruption' 'Flare' 'Oscillation'};
exlen=length(char(exten))+1;
plt=0;       %0 for plots , 1 for no plots
%If there is need to re-run a section  (DEVELOPER MODE)  Default all should be 1
weka_write=1;                       %0 for no weka writting, 1 for weka writing
%Image Parameter List 
imgParam={'Entropy' 'Mean' 'StandardDev' 'FracDim' 'Moment3' 'Moment4' 'Uniformity' 'RelSm' 'TamDir' 'TamCont'};
%%%%%%%%%%%%%%%%%%%%%% END OF GLOBAL VARIABLES 
%% Get Data from the FE files to manipulate
Histogram_Array=zeros(parameters,Csize*clases,segments*segments);  %Pre-allocate
imageCounter=1;
clear labelsData;
for fv=1:clases
folder=char([strcat(pathh,dataSet,'\',clasesNames(fv),'\*',exten)]);
folder2=char([strcat(pathh,dataSet,'\',clasesNames(fv),'\')]);
list=dir(folder);
    for ii=1:class_n
        filenam=strcat(folder2,list(ii,1).name);
        s2=length(char(filenam));
        filenam=strcat(filenam(1:s2-exlen),'-',int2str(segments),'x',int2str(segments),'.txt');
        [sample2(:,1),sample2(:,2),sample2(:,3),sample2(:,4),sample2(:,5),sample2(:,6),sample2(:,7),sample2(:,8),sample2(:,9),sample2(:,10),sample2(:,11)] = textread(filenam,'%s %s %s %s %s %s %s %s %s %s %s','delimiter', '\t','bufsize', 8191);
        for sgems=1:cells
            for parmsts=1:10
                Histogram_Array(parmsts,imageCounter,sgems)=str2double(sample2(sgems,parmsts));
            end
            labelsData(imageCounter)=sample2(sgems,11);
        end
        imageCounter=imageCounter+1;
    end %Looop files in class
end %loop classes
%% Manipulation of Feature Vectors
%Parameters loop
for curr_param=1:parameters  %Change to parameters
    skip=0;
    clear m;
    tmpH=zeros(total_im,cells);
    for tmx=1:total_im
        for tmy=1:cells
            tmpH(tmx,tmy)=Histogram_Array(curr_param,tmx,tmy);
        end
    end
    % Histogram sum to 1 Normalization (for crashing KLD and such)
    for nmH=1:total_im
        intBy=sum(tmpH(nmH,:));
        for tmy=1:cells
            tmpH(nmH,tmy)=tmpH(nmH,tmy)/intBy;
            if tmpH(nmH,tmy)==0
                tmpH(nmH,tmy)=0.000000000000000000000001;
            end
        end
    end
    
%Loop through Measures per parameter
for dnc=1:distances 
    skip=0;
    % For details on these measures please refer to Juan M. Banda's
    % Ph.D dissertation available at: http://www.jmbanda.com
    if dnc <= 8
        if dnc==1       %Calculate Euclidean Dissimilarity Measure (DM)
            m=pdist(tmpH(:,:),'euclidean');   
        end
        if dnc==2       %Calculate Standarized Euclidean DM
            m=pdist(tmpH(:,:),'seuclidean');  
        end
        if dnc==3       %Calculate Mahalanobis DM
            m=pdist(tmpH(:,:),'mahalanobis'); 
        end
        if dnc==4       %Calculate City Block DM
            m=pdist(tmpH(:,:),'cityblock');  
        end
        if dnc==5       %Calculate Cosine Distance DM
            m=pdist(tmpH(:,:),'cosine');  
        end
        if dnc==6       %Calculate Correlation DM
            m=pdist(tmpH(:,:),'correlation');
        end      
        if dnc==7       %Calculate Spearman DM
            m=pdist(tmpH(:,:),'spearman'); 
        end
        if dnc==8       %Calculate Chebychev DM
            m=pdist(tmpH(:,:),'chebychev');
        end        
    else
        if dnc==9     %Calculate Hausdorff DM
            Hausdorff_DM=zeros(total_im,total_im);
            %D2 = ones(size(tmpH(1,:)));  
            for xx=1:total_im
                for yy=xx:total_im
                    if xx==yy
                        Hausdorff_DM(yy,xx)=0;
                    else

                        for binsA=1:cells
                            for binsB=1:cells
                                dstM(binsB)=abs(tmpH(xx,binsA) - tmpH(yy,binsB));
                            end
                            relm1(binsA)=min(dstM);
                        end
                        HD1=max(relm1);
                        for binsA=1:cells
                            for binsB=1:cells
                                dstM(binsB)=abs(tmpH(yy,binsA) - tmpH(xx,binsB));
                            end
                            relm2(binsA)=min(dstM);
                        end
                        HD2=max(relm2);
                        HDD=max(HD1,HD2);                        
                        Hausdorff_DM(xx,yy)=HDD;
                        Hausdorff_DM(yy,xx)=Hausdorff_DM(xx,yy);
                    end
                end
            end
            tmc6=strcat(dncNM(dnc),imgParam(curr_param));
            tmc7=strcat(pathhTOS,dataSet,'-',tmc6,'-',int2str(segments),'x',int2str(segments),'.mat');
            save(char(tmc7),'Hausdorff_DM');            
            clear m;
            m=Hausdorff_DM;
        end
            if dnc==10      %Calculate KLD DM
              KLDM=zeros(total_im,total_im);
              for xx=1:total_im   
                Hist1=tmpH(xx,:);
                for yy=1:total_im
                    KLD=0;
                    Hist2=tmpH(yy,:);
                    if xx==yy
                       KLDM(xx,yy)=0;                
                    else
                        for m=1:cells
                            temp=Hist1(m)*log2((Hist1(m)/Hist2(m)));
                            if isnan(temp)
                                temp=0;
                            end
                            if isinf(temp)
                                temp=0;
                            end            
                            KLD=KLD+temp;
                        end
                        KLDM(xx,yy)=KLD;        
                    end
                end
              end
              %Save Variable
              tmc6=strcat(dncNM(dnc),imgParam(curr_param));
              tmc7=strcat(pathhTOS,dataSet,'-',tmc6,'-',int2str(segments),'x',int2str(segments),'.mat');
              save(char(tmc7),'KLDM');            
              clear m;
              m=KLDM;
            end
            if dnc==11       %Calculate JSD DM
              JDM=zeros(total_im,total_im);
              for xx=1:total_im   
                Hist1=tmpH(xx,:);
                for yy=xx:total_im
                    KLD=0;
                    Hist2=tmpH(yy,:);
                    if xx==yy
                       JDM(xx,yy)=0;                
                    else
                        JD=0;
                        for m=1:cells
                            temp=Hist1(m)*log2( ( (2*Hist1(m)) / (Hist1(m)+Hist2(m)) )) + Hist2(m)*log2( ( (2*Hist2(m)) / (Hist2(m)+Hist1(m)) ));
                            if isnan(temp)
                                temp=0;
                            end            
                            JD=JD+temp;
                        end
                        JDM(xx,yy)=JD;        
                        JDM(yy,xx)=JD;        
                    end
                end
              end
              tmc6=strcat(dncNM(dnc),imgParam(curr_param));
              tmc7=strcat(pathhTOS,dataSet,'-',tmc6,'-',int2str(segments),'x',int2str(segments),'.mat');
              save(char(tmc7),'JDM');                
              clear m;
              m=JDM;
            end            
            if dnc==12      %Calculate CH2 DM
              CH2_M=zeros(total_im,total_im);
              for xx=1:total_im   
                Hist1=tmpH(xx,:);
                for yy=xx:total_im
                    KLD=0;
                    Hist2=tmpH(yy,:);
                    if xx==yy
                       CH2_M(xx,yy)=0;                
                    else
                        CH2=0;
                        for m=1:cells
                            temp= (Hist1(m)-Hist2(m)) / (Hist1(m)+Hist2(m)) ;
                            if isnan(temp)
                                temp=0;
                            end               
                            CH2=CH2+temp;
                        end
                        CH2=abs(CH2);
                        CH2_M(xx,yy)=CH2;        
                        CH2_M(yy,xx)=CH2;        
                    end
                end
              end
              tmc6=strcat(dncNM(dnc),imgParam(curr_param));
              tmc7=strcat(pathhTOS,dataSet,'-',tmc6,'-',int2str(segments),'x',int2str(segments),'.mat');
              save(char(tmc7),'CH2_M');                
              clear m;
              m=CH2_M;
            end       
            
    end
        % MDS Calculation IF it fails output a warning
        clear YMDS;
        clear eigvals2;            
        try
            [YMDS,eigvals2] = cmdscale(m);
        catch
            [d1 d2]=size(m);
        if d1 ~= d2
            mN=squareform(m);  
        else
            mN=m;
        end
        for iX=1:total_im
            for iY=iX:total_im
                if isnan(mN(iX,iY))
                    mN(iX,iY)=0;
                    mN(iY,iX)=0;
                end
            end
        end
        m=mN;
        try
            [YMDS,eigvals2] = cmdscale(m);
        catch
            if dnc== 10   %Avoid sending a message with KLD measure until it is symmetrized and fails
                skip=1;
            else            %Not KLD? then output which matrix is bad
                oytd=strcat('bad matrix-',dataSet,' Distance - ',dncNM(dnc))
                skip=1;
            end
        end
        %plt=1;
        %skip=1;
        end
    %PLOTS       
    if plt==1 || skip==1
        %skip=0;
        %plt=0;
    else    
        %%2D MDS PLOT
        figure(1);
        set(1, 'Visible', 'off');
        plot(YMDS(:,1),YMDS(:,2),'LineStyle','none');
        axis(max(max(abs(YMDS))) * [-1.1,1.1,-1.1,1.1]); axis('square');
        text(YMDS(:,1),YMDS(:,2),labels,'HorizontalAlignment','left');
        hx = graph2d.constantline(0, 'LineStyle','-','Color',[.7 .7 .7]);
        hx = changedependvar(hx,'x');
        hy = graph2d.constantline(0,'LineStyle','-','Color',[.7 .7 .7]);
        hy = changedependvar(hy,'y');
        tmp=strcat(dataSet,' Distance - ',dncNM(dnc));
        tmp2=strcat(tmp,' MDS map for ');
        name=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
        tmpTLT1=char([strcat(dataSet,' Distance - ',dncNM(dnc))]);
        tmpTLT2=char([strcat(tmpTLT1,' MDS map for ')]);
        tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
        title(char(tmpTLT));        
        saveas(1,name,'jpg');
        close(1);
        %3D MSD Plots
        figure(6);
        set(6, 'Visible', 'off');
        plot3(YMDS(:,1),YMDS(:,2),YMDS(:,3),'LineStyle','none');
        text(YMDS(:,1),YMDS(:,2),YMDS(:,3),labels,'HorizontalAlignment','left');
        grid on;
        tmp=strcat(dataSet,' Distance - ',dncNM(dnc));
        tmp2=strcat(tmp,' 3D MDS map for ');
        name=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
        tmpTLT1=char([strcat(dataSet,' Distance - ',dncNM(dnc))]);
        tmpTLT2=char([strcat(tmpTLT1,' 3D MDS map for ')]);
        tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
        title(char(tmpTLT));        
        saveas(6,name,'jpg');
        close(6);      
        %% PLOT Components
        figure(2);
        set(2, 'Visible', 'off');
        [ymds1 ymds2]=size(YMDS);
        tempCMP=zeros(ymds2);
        tempCMP2=zeros(ymds2);
        for yy=1:ymds2
            tempCMP(yy)=max(abs(YMDS(:,yy)));
            tempCMP2(yy)=sum(abs(YMDS(:,yy)));
        end
        plot(tempCMP);
        tmp=strcat(dataSet,' Distance - ',dncNM(dnc));
        tmp2=strcat(tmp,' Components Plot ');
        nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
        tmpTLT1=char([strcat(dataSet,' Distance - ',dncNM(dnc))]);
        tmpTLT2=char([strcat(tmpTLT1,' Components Plot-')]);
        tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
        title(char(tmpTLT));    
        saveas(2,nameR,'jpg');
        close(2);
        %SUM of components PLOT
        figure(3);
        set(3, 'Visible', 'off');
        plot(tempCMP2);
        tmp=strcat(dataSet,' Distance - ',dncNM(dnc));
        tmp2=strcat(tmp,' Sum Of Components Plot-');
        nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
        tmpTLT1=char([strcat(dataSet,' Distance - ',dncNM(dnc))]);
        tmpTLT2=char([strcat(tmpTLT1,' Sum Of Components Plot-')]);
        tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);  
        title(char(tmpTLT)); 
        saveas(3,nameR,'jpg');
        close(3);
        [sr1 sr2]=size(m);
        if sr1 < sr2
            m2=squareform(m);
        else
            m2=m;
        end
        figure(4);
        set(4, 'Visible', 'off');
        imagesc(m2);
        colorbar;
        tmp=strcat(dataSet,' Distance - ',dncNM(dnc));
        tmp2=strcat(tmp,' Scaled Image Plot ');
        nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
        tmpTLT1=char([strcat(dataSet,' Distance - ',dncNM(dnc))]);
        tmpTLT2=char([strcat(tmpTLT1,' Scaled Image Plot ')]);
        tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
        title(char(tmpTLT));       
        saveas(4,nameR,'jpg');
        close(4);    
        figure(5);
        set(5, 'Visible', 'off');
        J = mat2gray(m2);
        imshow(J);
        axis on;
        colorbar;
        tmp=strcat(dataSet,' Distance - ',dncNM(dnc));
        tmp2=strcat(tmp,' Grayscaled Image Plot ');
        nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
        tmpTLT1=char([strcat(dataSet,' Distance - ',dncNM(dnc))]);
        tmpTLT2=char([strcat(tmpTLT1,' Grayscaled Image Plot ')]);
        tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
        title(char(tmpTLT));   
        saveas(5,nameR,'jpg');
        close(5);           
    end %Plots yes or no
    
    % Curve Fitting for component Thresholding
    FLNnameCP=strcat(dataSet,'-Distance-',dncNM(dnc),'-Feature-',imgParam(curr_param),'-',int2str(segments),'x',int2str(segments));    %Name for components plot
   if skip==1
       skip=0;
   else    
    clear x;
    syms x;
    [ymds1 ymds2]=size(YMDS);
    z=1:ymds2;
    tempCMP2=zeros(1,ymds2);
    for yy=1:ymds2
        tempCMP2(yy)=sum(abs(YMDS(:,yy)));
    end    
    xy=1:ymds2;
    f = ezfit(xy, tempCMP2, 'exp');
    if plt==1 || skip==1
        %skip=0;
        %plt=0;
    else  
        figure(6);
        set(6, 'Visible', 'off');
        plot(xy,tempCMP2,'x')
        hold on
        showfit(f,'fitcolor','blue');
        hold off
        tmpTLT=char([strcat('Exponential Curve fit for Sum Components',FLNnameCP)]);
        title(char(tmpTLT));   
        nameR=char([strcat(pathhTOS,tmpTLT)]);
        saveas(6,nameR,'jpg');
        close(6);
    end
        a=f.m(1);
        b=f.m(2);
        xp2=diff(a*exp((b*x)));
        %%Inline
        xpin2=inline(char(xp2));
        Slopes=zeros(1,ymds2);
        for xrt=1:ymds2
            Slopes(xrt)=xpin2(xrt);
        end
    if plt==1 || skip==1
        skip=0;
        %plt=0;
    else          
        figure(7);
        set(7, 'Visible', 'off');
        hold on;
        plot(Slopes)
        hold off;
        tmpTLT=char([strcat('Slopes of Exponential Curve fit for Sum Components',FLNnameCP)]);
        title(char(tmpTLT));   
        nameR=char([strcat(pathhTOS,tmpTLT)]);
        saveas(7,nameR,'jpg');
        close(7);            
    end
        %Determine number of components to get
        components_to_get=0;
        SlopesAT=zeros(1,ymds2);
        for sm=1:ymds2
            SlopesAT(sm)=radtodeg(atan(Slopes(sm)));
            if SlopesAT(sm) < 0 
                SlopesAT(sm)=180 + SlopesAT(sm);
            end    
            if SlopesAT(sm) < tang_thres
            else
                components_to_get=sm;
                break;
            end
        end        
        %If the number of components with the threshold is very big (or
        %all) we need select the number of components minus one
        if components_to_get==0
            components_to_get=cells-1;  %Worst Case scenario grab cells minus one
        end
        %Write this file to WEKA
        if weka_write==1
            clear tempData;
            clear featureNames;   
            SlopesNumberComp=components_to_get;
            for images=1:total_im
                for components=1:SlopesNumberComp
                    tempData(images,components)=num2cell(YMDS(images,components));
                end
                tempData(images,SlopesNumberComp+1)=labelsData(images);
            end
            clear tmc1;
            variable=1;
            for components=1:SlopesNumberComp
                tmp=num2str(variable);
                featureNames(components) = {tmp};
                variable=variable+1;
            end
            featureNames(SlopesNumberComp+1)={'label'};
            classindex=SlopesNumberComp+1; %Columns
            tmc1=strcat(dataSet,'-TangSlopeComponents-',num2str(SlopesNumberComp),'-Distance-',FLNnameCP);
            try 
                tempCon = matlab2weka(tmc1,featureNames,tempData,classindex);
            catch
                msg=strcat('Tangent Threshold WEKA Object Creation Failed for: ',FLNameCP)
            end
            tmc5=strcat(pathhTOS,tmc1,'.arff');
            try
                saveARFF(tmc5,tempCon);     
            catch
                msg=strcat('Tangent Threshold WEKA File saving failed for: ',FLNameCP)
                tmc5
            end            
            % Save the first N-Components (hard threshold)
            % NOTE: If you tweak the first side of the for loop you can
            % actually generate from 1 to N files with components, useful for
            % brute-force comparisons
            for componts=component_threshold:component_threshold
             SlopesNumberComp=componts;
             clear tempData;
             clear featureNames;   
             for images=1:total_im
                for components=1:SlopesNumberComp
                    tempData(images,components)=num2cell(YMDS(images,components));
                end
                tempData(images,SlopesNumberComp+1)=labelsData(images);
             end
             clear tmc1;
             variable=1;
             for components=1:SlopesNumberComp
                tmp=num2str(variable);
                featureNames(components) = {tmp};
                variable=variable+1;
             end
             featureNames(SlopesNumberComp+1)={'label'};
             classindex=SlopesNumberComp+1; %Columns
             tmc1=strcat(dataSet,'-SlopeComponents-',num2str(SlopesNumberComp),'-Distance-',FLNnameCP);
             try 
                 tempCon = matlab2weka(tmc1,featureNames,tempData,classindex);
             catch
                 msg=strcat('Hard Threshold WEKA Object Creation Failed for: ',FLNameCP)
             end
             tmc5=strcat(pathhTOS,tmc1,'.arff');
             try
                 saveARFF(tmc5,tempCon);     
             catch
                 msg=strcat('Hard Threshold WEKA File saving failed for: ',FLNameCP)
                tmc5
             end        
            end         %end hard components threshold
        end %Weka Write
   end  %Skip if matrix fails 
   if dnc==10   %KLD Treat it differently
       
        mU=triu(KLDM);   % Upper Triangular A-B
        mL=tril(KLDM);   % Lower Triangular B-A

        % Take the first part A-B
        m=zeros(total_im,total_im);
        for iX=1:total_im
            for iY=iX:total_im
                m(iX,iY)=real(mU(iX,iY));
                m(iY,iX)=real(mU(iX,iY));
            end
        end
        clear YMDS;
        clear eigvals2;
        [YMDS,eigvals2] = cmdscale(m);
        if plt==0        
            figure(1);
            set(1, 'Visible', 'off');
            plot(YMDS(:,1),YMDS(:,2),'LineStyle','none');
            axis(max(max(abs(YMDS))) * [-1.1,1.1,-1.1,1.1]); axis('square');
            text(YMDS(:,1),YMDS(:,2),labels,'HorizontalAlignment','left');
            hx = graph2d.constantline(0, 'LineStyle','-','Color',[.7 .7 .7]);
            hx = changedependvar(hx,'x');
            hy = graph2d.constantline(0,'LineStyle','-','Color',[.7 .7 .7]);
            hy = changedependvar(hy,'y');
            tmp=strcat(dataSet,' Distance - ','KLD A-B');
            tmp2=strcat(tmp,' MDS map for ');
            name=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD A-B')]);
            tmpTLT2=char([strcat(tmpTLT1,' MDS map for ')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
            title(char(tmpTLT));        
            saveas(1,name,'jpg');
            close(1);
            
            figure(8);
            set(8, 'Visible', 'off');
            plot3(YMDS(:,1),YMDS(:,2),YMDS(:,3),'LineStyle','none');
            text(YMDS(:,1),YMDS(:,2),YMDS(:,3),labels,'HorizontalAlignment','left');
            grid on;
            tmp=strcat(dataSet,' Distance - ','KLD A-B');
            tmp2=strcat(tmp,' 3D MDS map for ');
            name=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD A-B')]);
            tmpTLT2=char([strcat(tmpTLT1,' 3D MDS map for ')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
            title(char(tmpTLT));        
            saveas(8,name,'jpg');
            close(8);          
        end
        if weka_write==1 
            m2=m;
        %% PLOT Components
        if plt==0
            figure(2);
            set(2, 'Visible', 'off');
            tempCMP=zeros(cells-1);
            tempCMP2=zeros(cells-1);
            for yy=1:63
                tempCMP(yy)=max(abs(YMDS(:,yy)));
                tempCMP2(yy)=sum(abs(YMDS(:,yy)));
            end
            plot(tempCMP);
            tmp=strcat(dataSet,' Distance - ','KLD A-B');
            tmp2=strcat(tmp,' Components Plot ');
            nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD A-B')]);
            tmpTLT2=char([strcat(tmpTLT1,' Components Plot-')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
            title(char(tmpTLT));    
            saveas(2,nameR,'jpg');
            close(2);
            
            figure(3);
            set(3, 'Visible', 'off');
            plot(tempCMP2);
            tmp=strcat(dataSet,' Distance - ','KLD A-B');
            tmp2=strcat(tmp,' Sum Of Components Plot-');
            nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD A-B')]);
            tmpTLT2=char([strcat(tmpTLT1,' Sum Of Components Plot-')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);  
            title(char(tmpTLT)); 
            saveas(3,nameR,'jpg');
            close(3);    
   
            figure(4);
            set(4, 'Visible', 'off');
            imagesc(m2);
            colorbar;
            tmp=strcat(dataSet,' Distance - ','KLD A-B');
            tmp2=strcat(tmp,' Scaled Image Plot ');
            nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD A-B')]);
            tmpTLT2=char([strcat(tmpTLT1,' Scaled Image Plot ')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
            title(char(tmpTLT));       
            saveas(4,nameR,'jpg');
            close(4);    
            figure(5);
            set(5, 'Visible', 'off');
            J = mat2gray(m2);
            imshow(J);
            axis on;
            colorbar;
            tmp=strcat(dataSet,' Distance - ','KLD A-B');
            tmp2=strcat(tmp,' Grayscaled Image Plot ');
            nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD A-B')]);
            tmpTLT2=char([strcat(tmpTLT1,' Grayscaled Image Plot ')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
            title(char(tmpTLT));   
            saveas(5,nameR,'jpg');
            close(5);            
        end    
        %WEKA Writing sections
        FLNname=strcat(dataSet,'Distance-KLD A-B','-Feature-',imgParam(curr_param),'-',int2str(segments),'x',int2str(segments));
        clear x;
        syms x;
        z=1:cells-1;
        tempCMP2=zeros(1,cells-1);
        for yy=1:cells-1
            tempCMP2(yy)=sum(abs(YMDS(:,yy)));
        end    
        xy=1:cells-1;
        f = ezfit(xy, tempCMP2, 'exp');
        if plt==0        
            figure(6);
            set(6, 'Visible', 'off');
            plot(xy,tempCMP2,'x')
            hold on
            showfit(f,'fitcolor','blue');
            hold off
            tmpTLT=char([strcat(dataSet,'-Exponential Curve fit for Sum Components',FLNname)]);
            title(char(tmpTLT));   
            nameR=char([strcat(pathhTOS,tmpTLT)]);
            saveas(6,nameR,'jpg');
            close(6);
        end
        a=f.m(1);
        b=f.m(2);
        xp2=diff(a*exp((b*x)));
        %%Inline
        xpin2=inline(char(xp2));
        Slopes=zeros(1,cells-1);
        for xrt=1:cells-1
            Slopes(xrt)=xpin2(xrt);
        end
        if plt==0                
            figure(7);
            set(7, 'Visible', 'off');
            hold on;
            plot(Slopes)
            hold off;
            tmpTLT=char([strcat(dataSet,'-Slopes of Exponential Curve fit for Sum Components',FLNname)]);
            title(char(tmpTLT));   
            nameR=char([strcat(pathhTOS,tmpTLT)]);
            saveas(7,nameR,'jpg');
            close(7);        
        end
        components_to_get=0;
        SlopesAT=zeros(1,cells-1);
        for sm=1:cells-1
            SlopesAT(sm)=radtodeg(atan(Slopes(sm)));
            if SlopesAT(sm) < 0 
                SlopesAT(sm)=180 + SlopesAT(sm);
            end    
            if SlopesAT(sm) < tang_thres
            else
                components_to_get=sm;
                break;
            end
        end        
        if components_to_get==0
            components_to_get=cells-1;
        end
        
        clear tempData;
        clear featureNames;   
        SlopesNumberComp=components_to_get;
        for images=1:total_im
            for components=1:SlopesNumberComp
                tempData(images,components)=num2cell(YMDS(images,components));
            end
            tempData(images,SlopesNumberComp+1)=labelsData(images);
        end
        clear tmc1;
        variable=1;
        for components=1:SlopesNumberComp
            tmp=num2str(variable);
            featureNames(components) = {tmp};
            variable=variable+1;
        end
        featureNames(SlopesNumberComp+1)={'label'};
        classindex=SlopesNumberComp+1; %Columns
        tmc1=strcat(dataSet,'-TangSlopeComponents-',num2str(SlopesNumberComp),'-Distance-',FLNname);
        try 
            tempCon = matlab2weka(tmc1,featureNames,tempData,classindex);
        catch
            msg=char(strcat('Tangent Threshold WEKA Object Creation Failed for: ',FLName))
        end
        tmc5=strcat(pathhTOS,tmc1,'.arff');
        try
            saveARFF(tmc5,tempCon);     
        catch
            msg=char(strcat('Tangent Threshold WEKA Object Creation Failed for: ',FLName))
        end    
     
        %% N to N components only
       for componts=component_threshold:component_threshold
            SlopesNumberComp=componts;
            %%SLOPES 1
            clear tempData;
            clear featureNames;   
            %tempData=zeros(1600,Slopes2NumberComp);
            for images=1:total_im
                for components=1:SlopesNumberComp
                    tempData(images,components)=num2cell(YMDS(images,components));
                end
                tempData(images,SlopesNumberComp+1)=labelsData(images);
            end
        clear tmc1;
        variable=1;
        for components=1:SlopesNumberComp
            tmp=num2str(variable);
            featureNames(components) = {tmp};
            variable=variable+1;
        end
        featureNames(SlopesNumberComp+1)={'label'};
        classindex=SlopesNumberComp+1; %Columns
        tmc1=strcat(dataSet,'-ALL-SlopeComponents-',num2str(SlopesNumberComp),'-Distance-',FLNname);
        try 
            tempCon = matlab2weka(tmc1,featureNames,tempData,classindex);
        catch
            msg=char(strcat('Hard Threshold WEKA Object Creation Failed for: ',FLName))
        end
        tmc5=strcat(pathhTOS,tmc1,'.arff');
        try
            saveARFF(tmc5,tempCon);     
        catch
            msg=char(strcat('Hard Threshold WEKA Saving Failed for: ',FLNname))
        end        
       end  %END the N-Components
     end %Weka Option
    %B-A
    m=zeros(total_im,total_im);
    for iY=1:total_im
        for iX=1:iY
            m(iX,iY)=real(mL(iY,iX));
            m(iY,iX)=real(mL(iY,iX));
        end
    end
        clear YMDS;
        clear eigvals2;
        [YMDS,eigvals2] = cmdscale(m);
        if plt==0        
            figure(1);
            set(1, 'Visible', 'off');
            plot(YMDS(:,1),YMDS(:,2),'LineStyle','none');
            axis(max(max(abs(YMDS))) * [-1.1,1.1,-1.1,1.1]); axis('square');
            text(YMDS(:,1),YMDS(:,2),labels,'HorizontalAlignment','left');
            hx = graph2d.constantline(0, 'LineStyle','-','Color',[.7 .7 .7]);
            hx = changedependvar(hx,'x');
            hy = graph2d.constantline(0,'LineStyle','-','Color',[.7 .7 .7]);
            hy = changedependvar(hy,'y');
            tmp=strcat(dataSet,' Distance - ','KLD B-A');
            tmp2=strcat(tmp,' MDS map for ');
            name=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD B-A')]);
            tmpTLT2=char([strcat(tmpTLT1,' MDS map for ')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
            title(char(tmpTLT));        
            saveas(1,name,'jpg');
            close(1);
            
            figure(8);
            set(8, 'Visible', 'off');
            plot3(YMDS(:,1),YMDS(:,2),YMDS(:,3),'LineStyle','none');
            text(YMDS(:,1),YMDS(:,2),YMDS(:,3),labels,'HorizontalAlignment','left');
            grid on;
            tmp=strcat(dataSet,' Distance - ','KLD B-A');
            tmp2=strcat(tmp,' 3D MDS map for ');
            name=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD B-A')]);
            tmpTLT2=char([strcat(tmpTLT1,' 3D MDS map for ')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
            title(char(tmpTLT));        
            saveas(8,name,'jpg');
            close(8);          
        end
        if weka_write==1 
            m2=m;
        %% PLOT Components
        if plt==0
            figure(2);
            set(2, 'Visible', 'off');
            tempCMP=zeros(cells-1);
            tempCMP2=zeros(cells-1);
            for yy=1:63
                tempCMP(yy)=max(abs(YMDS(:,yy)));
                tempCMP2(yy)=sum(abs(YMDS(:,yy)));
            end
            plot(tempCMP);
            tmp=strcat(dataSet,' Distance - ','KLD B-A');
            tmp2=strcat(tmp,' Components Plot ');
            nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD B-A')]);
            tmpTLT2=char([strcat(tmpTLT1,' Components Plot-')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
            title(char(tmpTLT));    
            saveas(2,nameR,'jpg');
            close(2);
            
            figure(3);
            set(3, 'Visible', 'off');
            plot(tempCMP2);
            tmp=strcat(dataSet,' Distance - ','KLD B-A');
            tmp2=strcat(tmp,' Sum Of Components Plot-');
            nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD B-A')]);
            tmpTLT2=char([strcat(tmpTLT1,' Sum Of Components Plot-')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);  
            title(char(tmpTLT)); 
            saveas(3,nameR,'jpg');
            close(3);    
   
            figure(4);
            set(4, 'Visible', 'off');
            imagesc(m2);
            colorbar;
            tmp=strcat(dataSet,' Distance - ','KLD B-A');
            tmp2=strcat(tmp,' Scaled Image Plot ');
            nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD B-A')]);
            tmpTLT2=char([strcat(tmpTLT1,' Scaled Image Plot ')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
            title(char(tmpTLT));       
            saveas(4,nameR,'jpg');
            close(4);    
            figure(5);
            set(5, 'Visible', 'off');
            J = mat2gray(m2);
            imshow(J);
            axis on;
            colorbar;
            tmp=strcat(dataSet,' Distance - ','KLD B-A');
            tmp2=strcat(tmp,' Grayscaled Image Plot ');
            nameR=char([strcat(pathhTOS,tmp2,imgParam(curr_param),'-',int2str(segments),'x',int2str(segments))]);
            tmpTLT1=char([strcat(dataSet,' Distance - ','KLD B-A')]);
            tmpTLT2=char([strcat(tmpTLT1,' Grayscaled Image Plot ')]);
            tmpTLT=char([strcat(tmpTLT2,imgParam(curr_param))]);
            title(char(tmpTLT));   
            saveas(5,nameR,'jpg');
            close(5);            
        end    
        %WEKA Writing sections
        FLNname=strcat(dataSet,'Distance-KLD B-A','-Feature-',imgParam(curr_param),'-',int2str(segments),'x',int2str(segments));
        clear x;
        syms x;
        z=1:cells-1;
        tempCMP2=zeros(1,cells-1);
        for yy=1:cells-1
            tempCMP2(yy)=sum(abs(YMDS(:,yy)));
        end    
        xy=1:cells-1;
        f = ezfit(xy, tempCMP2, 'exp');
        if plt==0        
            figure(6);
            set(6, 'Visible', 'off');
            plot(xy,tempCMP2,'x')
            hold on
            showfit(f,'fitcolor','blue');
            hold off
            tmpTLT=char([strcat(dataSet,'-Exponential Curve fit for Sum Components',FLNname)]);
            title(char(tmpTLT));   
            nameR=char([strcat(pathhTOS,tmpTLT)]);
            saveas(6,nameR,'jpg');
            close(6);
        end
        a=f.m(1);
        b=f.m(2);
        xp2=diff(a*exp((b*x)));
        %%Inline
        xpin2=inline(char(xp2));
        Slopes=zeros(1,cells-1);
        for xrt=1:cells-1
            Slopes(xrt)=xpin2(xrt);
        end
        if plt==0                
            figure(7);
            set(7, 'Visible', 'off');
            hold on;
            plot(Slopes)
            hold off;
            tmpTLT=char([strcat(dataSet,'-Slopes of Exponential Curve fit for Sum Components',FLNname)]);
            title(char(tmpTLT));   
            nameR=char([strcat(pathhTOS,tmpTLT)]);
            saveas(7,nameR,'jpg');
            close(7);        
        end
        components_to_get=0;
        SlopesAT=zeros(1,cells-1);
        for sm=1:cells-1
            SlopesAT(sm)=radtodeg(atan(Slopes(sm)));
            if SlopesAT(sm) < 0 
                SlopesAT(sm)=180 + SlopesAT(sm);
            end    
            if SlopesAT(sm) < tang_thres
            else
                components_to_get=sm;
                break;
            end
        end        
        if components_to_get==0
            components_to_get=cells-1;
        end
        
        clear tempData;
        clear featureNames;   
        SlopesNumberComp=components_to_get;
        for images=1:total_im
            for components=1:SlopesNumberComp
                tempData(images,components)=num2cell(YMDS(images,components));
            end
            tempData(images,SlopesNumberComp+1)=labelsData(images);
        end
        clear tmc1;
        variable=1;
        for components=1:SlopesNumberComp
            tmp=num2str(variable);
            featureNames(components) = {tmp};
            variable=variable+1;
        end
        featureNames(SlopesNumberComp+1)={'label'};
        classindex=SlopesNumberComp+1; %Columns
        tmc1=strcat(dataSet,'-TangSlopeComponents-',num2str(SlopesNumberComp),'-Distance-',FLNname);
        try 
            tempCon = matlab2weka(tmc1,featureNames,tempData,classindex);
        catch
            msg=char(strcat('Tangent Threshold WEKA Object Creation Failed for: ',FLName))
        end
        tmc5=strcat(pathhTOS,tmc1,'.arff');
        try
            saveARFF(tmc5,tempCon);     
        catch
            msg=char(strcat('Tangent Threshold WEKA Object Creation Failed for: ',FLName))
        end    
     
        %% N to N components only
       for componts=component_threshold:component_threshold
            SlopesNumberComp=componts;
            %%SLOPES 1
            clear tempData;
            clear featureNames;   
            %tempData=zeros(1600,Slopes2NumberComp);
            for images=1:total_im
                for components=1:SlopesNumberComp
                    tempData(images,components)=num2cell(YMDS(images,components));
                end
                tempData(images,SlopesNumberComp+1)=labelsData(images);
            end
        clear tmc1;
        variable=1;
        for components=1:SlopesNumberComp
            tmp=num2str(variable);
            featureNames(components) = {tmp};
            variable=variable+1;
        end
        featureNames(SlopesNumberComp+1)={'label'};
        classindex=SlopesNumberComp+1; %Columns
        tmc1=strcat(dataSet,'-ALL-SlopeComponents-',num2str(SlopesNumberComp),'-Distance-',FLNname);
        try 
            tempCon = matlab2weka(tmc1,featureNames,tempData,classindex);
        catch
            msg=char(strcat('Hard Threshold WEKA Object Creation Failed for: ',FLName))
        end
        tmc5=strcat(pathhTOS,tmc1,'.arff');
        try
            saveARFF(tmc5,tempCon);     
        catch
            msg=char(strcat('Hard Threshold WEKA Saving Failed for: ',FLNname))
        end        
       end  %END the N-Components
     end %Weka Option    
       
   end  %KLD Distance Option
end  % distances Loop
end  %Parameters loop

msg='Dissimilarity Measure Module Demo has been completed. Check your output folder for results'