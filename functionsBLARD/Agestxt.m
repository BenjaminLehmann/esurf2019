function Agestxt( filename,A,ParamCell,textPR )
% Function that writes the CREp program ages results

% Delete commas
if strcmp(ParamCell{2},'LSD (Lifton et al., 2014)')==1;
    ParamCell{2}='LSD (Lifton et al 2014)';
end
if strcmp(ParamCell{3},'ERA40 (Uppala et al., 2005)')==1;
    ParamCell{3}='ERA40 (Uppala et al 2005)';
end
if strcmp(ParamCell{4},'Muscheler et al., 2005')==1;
    ParamCell{4}='Muscheler et al (2005)';
end

% First lines
Time=clock;
Heading=sprintf('CREp program results\n\n %d-%d-%d \n%dh%dmin\n',Time(1),Time(2),Time(3),Time(4),Time(5));
Param=sprintf('Parameters\nNuclide, %s\nScaling scheme, %s\nAtmosphere, %s\nMagnetic data, %s\nProduction rate, %s',ParamCell{1},ParamCell{2},ParamCell{3},ParamCell{4},textPR);

% Ages data
[Nl,~]=size(A);
% Names line
Names=sprintf('Samples, Scaling Factor, Age (ka), 1s (ka), 1s withour PR uncert., Status');

% First dataline
Agesdat=sprintf('%s, %s, %s, %s, %s, %s',A{1,1},A{1,2},A{1,3},A{1,4},A{1,5},A{1,6});

% Following lines
for i=2:Nl;
    Rowdat=sprintf('%s, %s, %s, %s, %s, %s',A{i,1},A{i,2},A{i,3},A{i,4},A{i,5},A{i,6});
    Agesdat=sprintf('%s \n%s',Agesdat,Rowdat);
end

Agesdat=sprintf('%s\n%s\n\n%s\n%s',Heading,Param,Names,Agesdat);

filename=strcat(filename,'.csv');
fileID=fopen(filename,'w');
fprintf(fileID,'%s',Agesdat);
fclose(fileID);

end

