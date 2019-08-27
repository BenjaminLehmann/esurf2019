function PDFtxt(filename,A,ParamCell,textPR)
% Function that writes the PDF csv document for the crep program

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
[Nl,Nc]=size(A);
% Names line
Names=A{1,1};
for i=2:Nc;
    Names=strcat(Names,', ',A{1,i});
end

% First dataline
Agesdat=sprintf('%4.3f',A{2,1});
for i=2:Nc;
    Agesdat=strcat(Agesdat,', ',sprintf('%4.3e',A{2,i}));
end

% Following lines
for i=3:Nl;
    Rowdat=sprintf('%4.3f',A{i,1});
    for j=2:Nc;
        Rowdat=strcat(Rowdat,', ',sprintf('%4.3e',A{i,j}));
    end
    Agesdat=sprintf('%s \n %s',Agesdat,Rowdat);
end

Agesdat=sprintf('%s\n%s\n\n%s\n%s',Heading,Param,Names,Agesdat);

filename=strcat(filename,'.csv');
fileID=fopen(filename,'w');
fprintf(fileID,'%s',Agesdat);
fclose(fileID);

end

