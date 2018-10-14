%% WEB UI - SQL statement generation for imageFARMER
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

%Extract Data from Histogram_Array and make it into SQL statements for the
%demo
% ---- Normalize Data First ----------- %
clear
sample_size = 200;
cells = 8 * 8;
parameters = 10;
create_tables=1;
create_values=1;
create_distance_matrix =1;

load('E:\SDO Research\Matlab CODE\Framework Demo Code\med-histogram_array_10x1600x64.mat');
stat= 'File Loaded!'
clear maxM;
clear minM;
maxM =zeros(10,1);
minM =zeros(10,1);
%Normalize  between 0 and 1
for parm=1:10
    maxM(parm)=max(max(Histogram_Array(parm,:,:)));  % Max of Image Param 1
    minM(parm)=min(min(Histogram_Array(parm,:,:)));  % Max of Image Param 1
end
norm=1;   % Change to skip this!!!

%Normal Alternative
for nmH=1:sample_size
    for tmy=1:cells
        for clse=1:parameters
            if (maxM(clse)+abs(minM(clse))) == 0
                Histogram_Array(clse,nmH,tmy)=0; %Handle Pesky Zeroes
            else
                Histogram_Array(clse,nmH,tmy)=(Histogram_Array(clse,nmH,tmy) + abs(minM(clse))) / (maxM(clse)+abs(minM(clse)));
            end
        end
    end  
end
param = 10;
rows=8;
cols=8;
stat = 'Normalized!'
dataSETNAME='MED05';



if create_tables==1     
    fil2=strcat('Tables-',dataSETNAME,'-',int2str(cells/8) ,'by',int2str(cells/8),'-1600-Param-ALL.txt');
    fid = fopen(fil2, 'w');    
    filnm= 'CREATE TABLE IF NOT EXISTS `params` ( `paramID` int(11) NOT NULL AUTO_INCREMENT, `paramName` varchar(100) NOT NULL, PRIMARY KEY (`paramID`)) ENGINE=MyISAM  DEFAULT CHARSET=latin1 AUTO_INCREMENT=11 ;';
    fprintf(fid, '%s', filnm);
    fprintf(fid, '\n');    
    filnm1=' INSERT INTO `params` (`paramID`, `paramName`) VALUES (1, ''Entropy''), (2, ''Mean''), (3, ''STDEV''), (4, ''FD''), (5, ''skewness''), (6, ''kurtosis''), (7, ''uniformity''), (8, ''rs''), (9, ''TDirectionality''), (10, ''TContrast'');';
    fprintf(fid, '%s', filnm1);
    fprintf(fid, '\n');    
    %Go through all the cells and create a huge string
    fill = 'CREATE TABLE IF NOT EXISTS `vals` (`valsID` int(11) NOT NULL AUTO_INCREMENT, `imageID` int(11) NOT NULL, `paramID` int(11) NOT NULL,';
    %Loop adding the rows and cols
    for rowSW= 0: rows-1
        for colSW = 0: cols-1
            fill = strcat(fill,'`R',int2str(rowSW),'-C',int2str(colSW),'` float NOT NULL, ');
        end
    end
  fill = strcat(fill, 'PRIMARY KEY (`valsID`)) ENGINE=MyISAM DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;');
  fprintf(fid, '%s', fill);
  fprintf(fid, '\n');
  
  fill = 'CREATE TABLE IF NOT EXISTS `dist` (`paramID` int(11) NOT NULL, `imageID` int(11) NOT NULL,';
    %Loop adding the rows and cols
    for rowSW= 1: 1600
            fill = strcat(fill,'`I',int2str(rowSW),'` float NOT NULL');
            if rowSW ~= 1600
                    fill=strcat(fill,',');
            end
    end
  fill = strcat(fill, ' ) ENGINE=MyISAM DEFAULT CHARSET=latin1 ;');  
  fprintf(fid, '%s', fill);
  fprintf(fid, '\n');  
  
  %Image O Table
  fill ='CREATE TABLE IF NOT EXISTS `images_o` ( `imageID` int(11) NOT NULL, `imagePath` varchar(150) NOT NULL, UNIQUE KEY `imageID` (`imageID`)) ENGINE=MyISAM DEFAULT CHARSET=latin1;';
  fprintf(fid, '%s', fill);
  fprintf(fid, '\n');    
  %fclose(fid);    
  %Image H Table
  fill ='CREATE TABLE IF NOT EXISTS `images_h` ( `imageID` int(11) NOT NULL, `paramID` int(11) NOT NULL, `imagePath` varchar(150) NOT NULL) ENGINE=MyISAM DEFAULT CHARSET=latin1;';
  fprintf(fid, '%s', fill);
  fprintf(fid, '\n');    
  %fclose(fid);      
  %Image H Table
  fill ='CREATE TABLE IF NOT EXISTS `images_v` ( `imageID` int(11) NOT NULL, `paramID` int(11) NOT NULL, `imagePath` varchar(150) NOT NULL) ENGINE=MyISAM DEFAULT CHARSET=latin1;';
  fprintf(fid, '%s', fill);
  fprintf(fid, '\n');    
  fclose(fid);       
  statu='Tables Created'
end 


if create_values ==1
%Start looping through each parameter and get SQL output
recrd=1;
for parm=1:10
    fil2=strcat('Values-',dataSETNAME,'-',int2str(cells/8) ,'by',int2str(cells/8),'-1600-Param-',int2str(parm),'.txt');
    fid = fopen(fil2, 'w');        
    for imgR= 1: 1600
        fill = strcat('INSERT INTO vals VALUES (',int2str(recrd),',',int2str(imgR),',',int2str(parm),',');
        for rowSW= 1: rows*cols
                fill = strcat(fill,num2str(Histogram_Array(parm,imgR,rowSW)));
                if rowSW ~= 64
                    fill=strcat(fill,',');
                end
        end
        fill = strcat(fill ,');');
        fprintf(fid, '%s', fill);
        fprintf(fid, '\n');          
        recrd=recrd+1;
    end
    fclose(fid);  
end

end
statu='Values Created'

if create_distance_matrix ==1
    
    for parm=1:10
        clear M;
        tmp=squeeze(Histogram_Array(parm,:,:));
        M=pdist(tmp);
        M=squareform(M);    
        fil2=strcat('Distances-',dataSETNAME,'-',int2str(cells/8) ,'by',int2str(cells/8),'-1600-Param-',int2str(parm),'.txt');
        fid = fopen(fil2, 'w');        
        for imgR= 1: 1600
            fill = strcat('INSERT INTO dist VALUES (',int2str(parm),',',int2str(imgR),',');
            for imgR2=1:1600
                fill = strcat(fill,num2str(M(imgR,imgR2)));
                if imgR2 ~= 1600
                    fill=strcat(fill,',');
                end
            end
            fill = strcat(fill ,');');
            fprintf(fid, '%s', fill);
            fprintf(fid, '\n');          
        end
        fclose(fid);  
        status = strcat('Param-',int2str(parm),'-done');
    end
end
statu='Distances Created'



