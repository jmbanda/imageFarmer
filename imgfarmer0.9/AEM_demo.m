%% Attribute Evaluation Module (AEM) - DEMO
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
%  Stand-alone script to ilustrate the usage of the Attribute Evaluation
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

clear;
rehash;    
%%%%% Main Variables
parameters=10;  %Image Parameters
imgParam={'Entropy' 'Mean' 'StandardDev' 'FracDim' 'Moment3' 'Moment4' 'Uniformity' 'RelSm' 'TamDir' 'TamCont'};
clases=8;
Csize=20;   %Number of files in each class
segments=8;
class_n=20;   
index=segments*segments*parameters;
dataSet='Solar';
clasesNames={'Active Region' 'Coronal Jet' 'Emerging Flux' 'Filament' 'Filament Activation' 'Filament Eruption' 'Flare' 'Oscillation'};
exten='tif';
exlen=length(char(exten))+1;
pathh='C:\Demo\';
cells=segments*segments;            %Grid Cell Size
%Save to path
pathhTOS='C:\Demo\';     
%% Actual Code
%Manipulation of feature vector (using a function that
%reads the txt files of the Feature extraction to save up on loading a
%unique feature vector variable
Histogram_Array=zeros(parameters,Csize*clases,segments*segments);  %Pre-allocate
imageCounter=1;
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
        end
        imageCounter=imageCounter+1;
    end %Looop files in class
end %loop classes
intraClass=zeros(parameters,1);  %Pre-allocate
interClass=zeros(parameters,1);  %Pre-allocate

% Calculate intra and inter class correlations by looping through parameters
clear corrMat;
corrMat=zeros(clases,parameters,parameters);
for cls =1:clases
    clear tst;
    if cls==1      %SHIFT THE CLASS MARKER
        rng1=1;
        rng2=Csize;  %class size
    else
        rng1=rng1+Csize;
        rng2=rng2+Csize;
    end
    for parms=1:parameters    %Calculate the intraclass correlation of the first class for all 10 parameters
            tst(parms,:,:)=squareform(pdist(squeeze(Histogram_Array(parms,rng1:rng2,:)),'correlation'));  %Calculate intra class correlation Active Region            
    end
    %Now we need to generate the intra-class correlation matrix
    for matx=1:parameters
        for maty=1:parameters
            corrMat(cls,matx,maty)=abs(mean(mean(squeeze(tst(matx,:,:)-tst(maty,:,:)))));   %populate the matrix
            corrMat(cls,matx,maty)=abs(1-corrMat(cls,matx,maty));
        end
    end
    
    %INTRA-CLASS PLOTS
    name=char([strcat(clasesNames(cls),' Intra-class Correlation-DS-',dataSet)]);
    tempName=char([clasesNames(cls) ' Intra-class Correlation-DS-' dataSet]);
    figure('NumberTitle','off','Name','Intra-class Correlation','Colormap',[1 1 1;0.9841 0.9841 0.9841;0.9683 0.9683 0.9683;0.9524 0.9524 0.9524;0.9365 0.9365 0.9365;0.9206 0.9206 0.9206;0.9048 0.9048 0.9048;0.8889 0.8889 0.8889;0.873 0.873 0.873;0.8571 0.8571 0.8571;0.8413 0.8413 0.8413;0.8254 0.8254 0.8254;0.8095 0.8095 0.8095;0.7937 0.7937 0.7937;0.7778 0.7778 0.7778;0.7619 0.7619 0.7619;0.746 0.746 0.746;0.7302 0.7302 0.7302;0.7143 0.7143 0.7143;0.6984 0.6984 0.6984;0.6825 0.6825 0.6825;0.6667 0.6667 0.6667;0.6508 0.6508 0.6508;0.6349 0.6349 0.6349;0.619 0.619 0.619;0.6032 0.6032 0.6032;0.5873 0.5873 0.5873;0.5714 0.5714 0.5714;0.5556 0.5556 0.5556;0.5397 0.5397 0.5397;0.5238 0.5238 0.5238;0.5079 0.5079 0.5079;0.4921 0.4921 0.4921;0.4762 0.4762 0.4762;0.4603 0.4603 0.4603;0.4444 0.4444 0.4444;0.4286 0.4286 0.4286;0.4127 0.4127 0.4127;0.3968 0.3968 0.3968;0.381 0.381 0.381;0.3651 0.3651 0.3651;0.3492 0.3492 0.3492;0.3333 0.3333 0.3333;0.3175 0.3175 0.3175;0.3016 0.3016 0.3016;0.2857 0.2857 0.2857;0.2698 0.2698 0.2698;0.254 0.254 0.254;0.2381 0.2381 0.2381;0.2222 0.2222 0.2222;0.2063 0.2063 0.2063;0.1905 0.1905 0.1905;0.1746 0.1746 0.1746;0.1587 0.1587 0.1587;0.1429 0.1429 0.1429;0.127 0.127 0.127;0.1111 0.1111 0.1111;0.09524 0.09524 0.09524;0.07937 0.07937 0.07937;0.06349 0.06349 0.06349;0.04762 0.04762 0.04762;0.03175 0.03175 0.03175;0.01587 0.01587 0.01587;0 0 0]);
    set(gcf, 'Visible', 'off'); 
    labelsY=1:parameters;
    axes1=axes('YTickLabel',imgParam,'YDir','reverse','XTick',labelsY,'XAxisLocation','top','Layer','top');
    xlim([0.5 10.5]);
    ylim([0.5 10.5]);
    hold('all');
    image(squeeze(corrMat(cls,:,:)),'Parent',axes1,'CDataMapping','scaled');
    colorbar('peer',axes1);    
    title([clasesNames(cls) 'Intra-class Correlation-DS-' dataSet]);
    saveas(gcf,strcat(pathhTOS,name),'jpg');
    close(gcf);    
    msg=strcat(name,' plot has been saved')
    
    tmp2=squeeze(corrMat(cls,:,:));
    %Calculate MDS for the MDS 2 component plots
    [Y2,eigvals2] = cmdscale(tmp2);
    for i=1:parameters
        labels(i) = {int2str(i)};
    end
    set(gca, 'Visible', 'off');
    plot(Y2(:,1),Y2(:,2),'LineStyle','none');
    axis(max(max(abs(Y2))) * [-1.1,1.1,-1.1,1.1]); 
    axis('square');
    text(Y2(:,1),Y2(:,2),labels,'HorizontalAlignment','left');
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    name=char(strcat(clasesNames(cls),' MDS Intra-class Correlation-DS-',dataSet));
    title([clasesNames(cls) ' MDS Intra-class Correlation-DS-' dataSet]);
    saveas(gca,strcat(pathhTOS,name),'jpg');
    msg=strcat(name,' plot has been saved')
    %close(2);    

end   %Classes for Intraclass

%For the interclass correlation we create blocks of the classes
mx=1;
blocks=zeros(parameters,clases,Csize,Csize);
for imgPrm=1:parameters
    cls=squareform(pdist(squeeze(Histogram_Array(imgPrm,:,:)),'correlation'));
    for clss=1:clases
        if clss==1      %SHIFT THE CLASS MARKER
            rng1=1;
            rng2=Csize;  %class size
        else
            rng1=rng1+Csize;
            rng2=rng2+Csize;
        end
        blocks(imgPrm,clss,:,:)=cls(rng1:rng2,rng1:rng2);
    end   %this gives me all blocks of classes for the parameters
end

%Create the interclass corelation creating the blocks for each class
corrMat=zeros(clases,parameters);
for clss=1:clases    %clases  
    for imgPrm=1:parameters  %
        averga=0;
        for clss2=1:clases  %class versus rest
            averga=averga+abs((mean(mean(squeeze(blocks(imgPrm,clss,:,:))))-mean(mean(squeeze(blocks(imgPrm,clss2,:,:))))));
        end
        corrMat(clss,imgPrm)=averga/clases;
    end
end %loop classes

%Match the classes between each other
corrMat2=zeros(clases,parameters,parameters);
for cls=1:clases
    for matx=1:parameters
        for maty=1:parameters
            corrMat2(cls,matx,maty)=abs(mean(mean(squeeze(corrMat(cls,matx)-corrMat(cls,maty)))));   %populate the matrix
            corrMat2(cls,matx,maty)=abs(1-corrMat2(cls,matx,maty));
        end
    end
    
    %Create Figure plots
    name=char([strcat(clasesNames(cls),' Inter-class Correlation-DS-',dataSet)]);
    tempName=char([clasesNames(cls) ' Inter-class Correlation-DS-' dataSet]);
    figure('NumberTitle','off','Name','Inter-class Correlation','Colormap',[1 1 1;0.9841 0.9841 0.9841;0.9683 0.9683 0.9683;0.9524 0.9524 0.9524;0.9365 0.9365 0.9365;0.9206 0.9206 0.9206;0.9048 0.9048 0.9048;0.8889 0.8889 0.8889;0.873 0.873 0.873;0.8571 0.8571 0.8571;0.8413 0.8413 0.8413;0.8254 0.8254 0.8254;0.8095 0.8095 0.8095;0.7937 0.7937 0.7937;0.7778 0.7778 0.7778;0.7619 0.7619 0.7619;0.746 0.746 0.746;0.7302 0.7302 0.7302;0.7143 0.7143 0.7143;0.6984 0.6984 0.6984;0.6825 0.6825 0.6825;0.6667 0.6667 0.6667;0.6508 0.6508 0.6508;0.6349 0.6349 0.6349;0.619 0.619 0.619;0.6032 0.6032 0.6032;0.5873 0.5873 0.5873;0.5714 0.5714 0.5714;0.5556 0.5556 0.5556;0.5397 0.5397 0.5397;0.5238 0.5238 0.5238;0.5079 0.5079 0.5079;0.4921 0.4921 0.4921;0.4762 0.4762 0.4762;0.4603 0.4603 0.4603;0.4444 0.4444 0.4444;0.4286 0.4286 0.4286;0.4127 0.4127 0.4127;0.3968 0.3968 0.3968;0.381 0.381 0.381;0.3651 0.3651 0.3651;0.3492 0.3492 0.3492;0.3333 0.3333 0.3333;0.3175 0.3175 0.3175;0.3016 0.3016 0.3016;0.2857 0.2857 0.2857;0.2698 0.2698 0.2698;0.254 0.254 0.254;0.2381 0.2381 0.2381;0.2222 0.2222 0.2222;0.2063 0.2063 0.2063;0.1905 0.1905 0.1905;0.1746 0.1746 0.1746;0.1587 0.1587 0.1587;0.1429 0.1429 0.1429;0.127 0.127 0.127;0.1111 0.1111 0.1111;0.09524 0.09524 0.09524;0.07937 0.07937 0.07937;0.06349 0.06349 0.06349;0.04762 0.04762 0.04762;0.03175 0.03175 0.03175;0.01587 0.01587 0.01587;0 0 0]);
    set(gcf, 'Visible', 'off'); 
    labelsY=1:parameters;
    axes1=axes('YTickLabel',imgParam,'YDir','reverse','XTick',labelsY,'XAxisLocation','top','Layer','top');
    xlim([0.5 10.5]);
    ylim([0.5 10.5]);
    hold('all');
    image(squeeze(corrMat2(cls,:,:)),'Parent',axes1,'CDataMapping','scaled');
    colorbar('peer',axes1);
    title([clasesNames(cls) 'Inter-class Correlation-DS-' dataSet]);
    saveas(gcf,strcat(pathhTOS,name),'jpg');
    close(gcf);     
    msg=strcat(name,' plot has been saved')
    
    %Calculate MDS for the MDS 2 component plots
    tmp2=squeeze(corrMat2(cls,:,:));
    [Y2,eigvals2] = cmdscale(tmp2);
    for i=1:parameters
        labels(i) = {int2str(i)};
    end
    set(gca, 'Visible', 'off');
    plot(Y2(:,1),Y2(:,2),'LineStyle','none');
    axis(max(max(abs(Y2))) * [-1.1,1.1,-1.1,1.1]); 
    axis('square');
    text(Y2(:,1),Y2(:,2),labels,'HorizontalAlignment','left');
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    name=char(strcat(clasesNames(cls),' MDS Inter-class Correlation-DS-',dataSet));
    title([clasesNames(cls) ' MDS Inter-class Correlation-DS-' dataSet]);
    saveas(gca,strcat(pathhTOS,name),'jpg');
    msg=strcat(name,' plot has been saved')
end %loop clases
msg='Attribute Evaluation Demo has been completed. Check your output folder for results'