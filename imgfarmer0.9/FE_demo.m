%% Feature Extraction Module (FEM) - DEMO
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
% MORE INFO:
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
%  http://www.imagefarmer.org

clear;
rehash;    
%%%%% Main Variables
clases=8;                           %Number of classes
Csize=200;                            %Number of images in each folder
segments=8;                        %Number of Columns/rows  NbyN
class_n=200;                         %class size
sample_size=clases*class_n;         %TOTAL number of images
cells=segments*segments;            %Grid Cell Size
parameters=10;                      %Number of Image Parameters
index=segments*segments*parameters; %Number of elements in the Feature Vector
dataSet='';                    %Dataset Being Manipulated
exten='tif';
%dataSet='ClefMed05';
%exten='png';
pathh='f:\SDO Research\Datasets\TRACE - Versions\Original\';                   %Root folder of datasets
pathhTOS='f:\SDO Research\Datasets\TRACE - Versions\Original\';     %Path to SAVE files to (except FE files)
%solar Class names (folder names)
clasesNames={'Active Region' 'Coronal Jet' 'Emerging Flux' 'Filament' 'Filament Activation' 'Filament Eruption' 'Flare' 'Oscillation'};

%If there is need to re-run a section  (DEVELOPER MODE)  Default all should be 1
prm_vis=0;                          %0 for no Visualization, 1 for visualization
feat_ext=1;                         %0 for no feature extraction, 1 for feature extraction
weka_write=0;                       %0 for no weka writting, 1 for weka writing

%Image Parameter List 
imgParam={'Entropy' 'Mean' 'StandardDev' 'FracDim' 'Moment3' 'Moment4' 'Uniformity' 'RelSm' 'TamDir' 'TamCont'};
i=0;
%Feature Vector Structure
if feat_ext==1
clear feature_vector;
clear Histogram_Array;
feature_vector = struct('number', {}, 'name', {}, 'entropy', {}, 'mean' , {}, 'standarddev' ,{},'HS_3moment'   ,{},'HS_4moment'   ,{},'HS_uniformity'   ,{},'HS_rs',{}, 'fractal_dimension', {}, 'segment', {}, 'Tamura_directionality', {} ,'Tamura_contrast' , {},'label'   ,{});
for fv=1:clases
folder=char([strcat(pathh,dataSet,'\',clasesNames(fv),'\*',exten)]);
folder2=char([strcat(pathh,dataSet,'\',clasesNames(fv),'\')]);
list=dir(folder);
for ii=1:class_n;
    i=i+1;
	I=imread([folder2 list(ii,1).name]);
    imageInfo = imfinfo([folder2 list(ii,1).name]);
    name=list(ii,1).name;
    s2=length(char(name));
    tmc7=strcat(folder2,list(ii,1).name(1:s2-4),'-',int2str(segments),'x',int2str(segments),'.txt');    
    fid = fopen(tmc7,'w');    
    % Create variable size grids
    im_y=imageInfo.Height;
    im_x=imageInfo.Width;
    chunkX = round(im_x/segments);
    chunkY = round(im_y/segments);
    
    x=1;
    y=1;   
    x2=chunkX;
    y2=chunkY;
    for ir=1:segments; %rows 
     for iCol=1:segments;  %columns         
            feature_vector(fv,ii,(ir*segments)-(segments-iCol)).number=i;
            feature_vector(fv,ii,(ir*segments)-(segments-iCol)).name=list(ii,1).name;
            TEMP=I(y:y2,x:x2);  
            feature_vector(fv,ii,(ir*segments)-(segments-iCol)).segment=[iCol, ir];
            feature_vector(fv,ii,(ir*segments)-(segments-iCol)).entropy=entropy(TEMP);
            fprintf(fid, '%f\t', feature_vector(fv,ii,(ir*segments)-(segments-iCol)).entropy);
            feature_vector(fv,ii,(ir*segments)-(segments-iCol)).mean=mean2(TEMP);
            fprintf(fid, '%f\t', feature_vector(fv,ii,(ir*segments)-(segments-iCol)).mean);
            feature_vector(fv,ii,(ir*segments)-(segments-iCol)).standarddev=std2(TEMP);    
            fprintf(fid, '%f\t', feature_vector(fv,ii,(ir*segments)-(segments-iCol)).standarddev);              
            feature_vector(fv,ii,(ir*segments)-(segments-iCol)).label=clasesNames(fv);
            feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension=new_bc(TEMP,3,0);
            if isnan(feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension) 
              feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension=1.646;
            else
                if isfinite(feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension)
                    if feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension < 0 
                        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension=0;
                    else
                        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension=feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension;
                    end
                else
                    if feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension < 0 
                        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension=0;
                    else
                        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension=1.646;
                    end
                end
            end
            fprintf(fid, '%f\t', feature_vector(fv,ii,(ir*segments)-(segments-iCol)).fractal_dimension);
        p=imhist(TEMP);
        p=p./numel(TEMP) ;
        L=length(p);
        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).HS_3moment=skewness(p);
        fprintf(fid, '%f\t',feature_vector(fv,ii,(ir*segments)-(segments-iCol)).HS_3moment);
        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).HS_4moment=kurtosis(p);   
        fprintf(fid, '%f\t',feature_vector(fv,ii,(ir*segments)-(segments-iCol)).HS_4moment);
        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).HS_uniformity=sum(p.^ 2);           
        fprintf(fid, '%f\t',feature_vector(fv,ii,(ir*segments)-(segments-iCol)).HS_uniformity); 
        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).HS_rs=1-(1/ (1+(std(p)^2))); 
        fprintf(fid, '%f\t', feature_vector(fv,ii,(ir*segments)-(segments-iCol)).HS_rs);          
        try   %Error handling
            temp_tamp=TamuraTextures(TEMP);
        %Individual
        catch
            %If it fails.... make zeros
            temp_tamp=[0 0];
        end
        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).Tamura_directionality = temp_tamp(1);                        
        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).Tamura_contrast = temp_tamp(2);
        fprintf(fid, '%f\t',temp_tamp(1));                        
        fprintf(fid, '%f\t',temp_tamp(2));        
        feature_vector(fv,ii,(ir*segments)-(segments-iCol)).label = clasesNames(fv);
        fprintf(fid, '%s\t',char(feature_vector(fv,ii,(ir*segments)-(segments-iCol)).label)); 

        x=x2;
        x2=x2+chunkX;
        if x2>im_x
            x2=im_x;
        end     
        fprintf(fid, '\n');
     end % columns
        y=y2;
        y2=y2+chunkY;
        %Check for out of bounds in rows
        if y2>im_y
            y2=im_y;
        end
        %Initialize column
        x=1;
        x2=chunkX;
    end % rows
    imgDone=i
    fclose(fid);
end %class elements number
end %class in data set
msg='Image Feature Extraction Done...'
end %Feature Extraction

if weka_write==1
%%Save WEKA file
imageCounter=1;
imageCounter2=1;
TEST=zeros(sample_size,index);
for classNumber=1:clases
    for imageNumber=1:class_n
        smt=0;
        for segment=1:cells
            Histogram_Array(1,imageCounter,segment)=feature_vector(classNumber,imageNumber,segment).entropy;
            Histogram_Array(2,imageCounter,segment)=feature_vector(classNumber,imageNumber,segment).mean;
            Histogram_Array(3,imageCounter,segment)=feature_vector(classNumber,imageNumber,segment).standarddev;
            Histogram_Array(4,imageCounter,segment)=feature_vector(classNumber,imageNumber,segment).fractal_dimension;
            Histogram_Array(5,imageCounter,segment)=feature_vector(classNumber,imageNumber,segment).HS_3moment;
            Histogram_Array(6,imageCounter,segment)=feature_vector(classNumber,imageNumber,segment).HS_4moment;
            Histogram_Array(7,imageCounter,segment)=feature_vector(classNumber,imageNumber,segment).HS_uniformity;
            Histogram_Array(8,imageCounter,segment)=feature_vector(classNumber,imageNumber,segment).HS_rs;
            Histogram_Array(9,imageCounter,segment)=feature_vector(classNumber,imageNumber,segment).Tamura_directionality;
            Histogram_Array(10,imageCounter,segment)=feature_vector(classNumber,imageNumber,segment).Tamura_contrast;
            if segment==cells
            else
                smt=segment*parameters; 
            end
        end
        imageCounter2=imageCounter2+1;
        imageCounter=imageCounter+1;
    end
end

    %tmpH  is the traning Data
    imag=1;
    cnt=1;
    imsr=1;
    tm=1;
    tmpH=zeros(sample_size,index);
    for imgk=1:sample_size
        for tmy=1:cells
            for cont=1:parameters
                tmpH(imsr,(((tmy-1)*parameters)+cont))=Histogram_Array(cont,imag,tmy);
            end
        end      
        if cnt==(Csize+1)
            tm=tm+1;
            cnt=1;
        end
        if imag <=(Csize*tm)
            labelsData(imsr)=clasesNames(tm);
            cnt=cnt+1;
        end
        imsr=imsr+1;
        imag=imag+1;
    end
    
    clear tempData;
    clear featureNames;   
    NumberComp=index;
    %size(mappedX,1)
    for images=1:sample_size
        for components=1:NumberComp
            tempData(images,components)=num2cell(tmpH(images,components));
            tempData(images,NumberComp+1)=labelsData(images);
        end
    end
    clear tmc1;
    variable=1;
    for components=1:NumberComp
        tmp=num2str(variable);
        featureNames(components) = {tmp};
        variable=variable+1;
    end
    featureNames(NumberComp+1)={'label'};
    classindex=NumberComp+1; %Columns
    
    tmc1=strcat('ALL-',dataSet,'-Feature-ALL-',char(num2str(segments)),'by',char(num2str(segments)));
    try 
        tempCon = matlab2weka(tmc1,featureNames,tempData,classindex);
    catch
        msg='RE-Save Weka File'
    end
    tmc5=strcat(pathhTOS,tmc1,'.arff');

    try
        saveARFF(tmc5,tempCon);     
    catch
        msg='RE-Save Weka File'
        tmc5
    end
msg='WEKA File Saved...'
end  %Write Weka File

%PARAMETER VISUALIZATION
if prm_vis==1
 for IMAGE_C=1:sample_size
  for param=1:parameters
    pMatrix=zeros(segments,segments);
    cellCnt=1;
    for xC=1:segments
        for yC=1:segments
            pMatrix(xC,yC)=Histogram_Array(param,IMAGE_C,cellCnt);
            cellCnt=cellCnt+1;
        end
    end
    figure(1);
    set(1, 'Visible', 'off');
    imagesc(pMatrix);
    titul=char([dataSet,'-ImageParam-',char(num2str(param)),'-imageno-',char(num2str(IMAGE_C))]);
    title(titul)
    experiment=char('ScaledImagePlotImageParam');
    name=char([strcat(pathhTOS,experiment'',num2str(param),'-',char(num2str(segments)),'by',char(num2str(segments)),'-',dataSet,'-imageno-',char(num2str(IMAGE_C)))]);
    saveas(1,strcat(name),'jpg');
    close(1);  
  end
 end
end   %PRM_VIS
msg='Feature Extraction Demo has been completed. Check your output folder for results'