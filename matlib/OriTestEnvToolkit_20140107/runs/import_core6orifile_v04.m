%%automatize analysis
clear;
clc;

iter_start = 0;
iter_end = 0;

for ii=iter_start:iter_end
 
%% init workspace
%clear;
%clc; %clear your screen and the console
 
 %% enable profiling of this script
optProfile = 0; % 0 - off, 1 - on;

if (optProfile == 1) 
    profile on;
end;

%% Import CORE6.IMM input files and plot arbitrary sections

% define symmetry
CS = symmetry('m-3m');
SS = symmetry('-1');
%SS = symmetry('mmm');

plotx2north
resol = 2.5*degree;
psi = kernel('de la Vallee Poussin','halfwidth',resol);

% define plotting phi2 sections
phi2_1 = 0*degree;
phi2_2 = 45*degree;
phi2_3 = 65*degree;

isolevel_min = 1;
isolevel_max = 20;
isolevels10 = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10];
isolevels20 = [1 2 3 4 5 6 7 8 9 10 12 14 16 18 20];
% ideal orientation resolution
optVolRes = 10*degree;

% plot also orthorhombic ODF?
optPlotOrtho = 0; % 0 - off, 1 - on;

%% path to file which should be imported
pname = ['\\cifs\cluster\Home\mk958655\OriTestEnvToolkit\build\'];
fname = ['randomOriShoemake.IMM.ori'];

%fname = [pname 'oribucket.' num2str(ii) '.304.IMM.ori'];
    
%% END OF USER INTERACTION

%% our beloved "Standardlagen"
%qref = euler2quat(59.0*degree,36.7*degree,63.4*degree,'Bunge');
Cube = orientation('Cube',CS,SS);
CubeTwin = orientation('Euler',45*degree,75*degree,45*degree,CS,SS);
CubeND22 = orientation('CubeND22',CS,SS);
CubeND22b = orientation('CubeND22b',CS,SS);
CubeRD22 = orientation('CubeRD22',CS,SS);
CubeRD22b = orientation('CubeRD22b',CS,SS);
CubeTD22 = orientation('CubeTD22',CS,SS);
CubeTD22b = orientation('CubeTD22b',CS,SS);
% CubeND45 = C orientation in bcc
CubeND45 = orientation('CubeND45',CS,SS);
% In bcc Copper = D1,2 orientations, also close to E, and -E (shear)
Copper = orientation('Copper',CS,SS);
Copper2 = orientation('Copper2',CS,SS);
%CopperGIA = orientation('CopperGIA',CS,SS);
%CopperGIA2 = orientation('CopperGIA2',CS,SS);
S = orientation('S',CS,SS);
S2 = orientation('S2',CS,SS);
S3 = orientation('S3',CS,SS);
S4 = orientation('S4',CS,SS);
R = orientation('R',CS,SS);
R2 = orientation('R2',CS,SS);
R3 = orientation('R3',CS,SS);
R4 = orientation('R4',CS,SS);
%SHirsch = orientation('SHirsch',CS,SS);
%SHirsch2 = orientation('SHirsch2',CS,SS);
%SHirsch3 = orientation('SHirsch3',CS,SS);
%SHirsch4 = orientation('SHirsch4',CS,SS);
% according to 1985_hirsch and 1983_itolücke
% R_Al {123}<63-4>
% S_Al {236}<32-2>  --> SHirsch
Brass = orientation('Brass',CS,SS);
Brass2 = orientation('Brass2',CS,SS);
% In bcc also called F orientation
Goss = orientation('Goss',CS,SS);
PLage = orientation('PLage',CS,SS);
PLage2 = orientation('PLage2',CS,SS);
QLage = orientation('QLage',CS,SS);
QLage2 = orientation('QLage2',CS,SS);
QLage3 = orientation('QLage3',CS,SS);
QLage4 = orientation('QLage4',CS,SS);
invGoss = orientation('invGoss',CS,SS);
% Dillamore orientation more suitable for Al than Copper orientation
Dill = orientation('Dill',CS,SS);
Dill2 = orientation('Dill2',CS,SS);
CTwin = orientation('Euler',15.4*degree,47.2*degree,67.9*degree,CS,SS);


%% Import the text data, 
% format: phi1, PHI, phi2, normalized weight, ignored value
% all angles in degree

% import data
newData1 = importdata(fname);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

noris = length(data);
%noris = 2; %10;

%MANUAL ADJUSTMENT OF THE DATA

%quat = zeros(1,noris);
%weights = zeros(1,noris);

for j=1:noris
    if (mod(j,100) == 0) %plot to estimate where you are...
       j
    end  
    quat(j) = euler2quat(data(j,1),data(j,2),data(j,3),'Bunge');
    weights(j) = data(j,4);  
end;

%% generate SO3 grid from your quaternion set
s3g = SO3Grid(quat,CS,SS);

%% define, e.g. calculate ODF directly
%odf = ODF(s3g,weights,psi,CS,SS);

odf = calcODF(s3g,'weight',weights,'halfwidth',5*degree); %corrects for ghosts automatically

%odf_mmm = set(odf_det16,'SS',symmetry('mmm')); %force ODF to have %orthorhombic sample symmetry

%% inject other unimodal, uniform or fiber odfs here if you like...
%odf = uniformODF(CS,SS); %uniform

%% plot tricline ohi2 slices
%figure;
%plotodf(odf,'sections',1,'gray','contourf',[1 2 5 10 25 50 100 200 300 500],...
%        'projection','plain','fontsize',8,'PHI2',phi2_1,'Bunge','grid','minmax','off'); 
%annotate(qref,'fontsize',10,'gray','marker','s','MarkerFaceColor','red'); 
%colorbar;
%savefigure([fname '.phi2.0deg.tricline.png']);
%close;

%figure;
%plotodf(odf,'sections',1,'gray','contourf',[1 2 5 10 25 50 100 200 300 500],...
%        'projection','plain','fontsize',8,'PHI2',phi2_2,'Bunge','grid','minmax','off'); 
%annotate(qref,'fontsize',10,'gray','marker','s','MarkerFaceColor','red'); %'label','S1'
%colorbar;
%savefigure([fname '.phi2.45deg.tricline.png']);
%close;

%figure;
%plotodf(odf,'sections',1,'gray','contourf',[1 2 5 10 25 50 100 200 300 500],...
 %       'projection','plain','fontsize',8,'PHI2',phi2_3,'Bunge','grid','minmax','off'); 
%annotate(qref,'fontsize',10,'gray','marker','s','MarkerFaceColor','red');
%colorbar;
%savefigure([fname '.phi2.65deg.tricline.png']);
%close;

%if (optPlotOrtho == 1) 
    %% plot orthorhombic phi2 slices 
%    odf_mmm = set(odf_det16,'SS',symmetry('mmm'));
    
    %odf_mmm = odf;

    figure;
    plotodf(odf,'phi2',[phi2_1]*degree,'contourf',isolevels20,'gray','projection','plain',...
      'fontsize',8,'PHI2',phi2_1,'Bunge','minmax','off','colorrange',[isolevel_min isolevel_max]);
    colorbar;
    savefigure([fname '.phi2.0deg.tricline.png']);
    close;
    
    
    figure;
    plotodf(odf,'phi2',[phi2_2]*degree,'contourf',isolevels20,'gray','projection','plain',...
      'fontsize',8,'PHI2',phi2_2,'Bunge','minmax','off','colorrange',[isolevel_min isolevel_max]);
    colorbar;
    savefigure([fname '.phi2.45deg.tricline.png']);
    close;    
    
    figure;
    plotodf(odf,'phi2',[phi2_3]*degree,'contourf',isolevels20,'gray','projection','plain',...
      'fontsize',8,'PHI2',phi2_3,'Bunge','minmax','off','colorrange',[isolevel_min isolevel_max]);
    colorbar;
    savefigure([fname '.phi2.65deg.tricline.png']);
    close;    
    
    
  %  figure;
  %  plotodf(odf_mmm,'sections',1,'gray','contourf',[1 2 5 10 25 50 100 200 300 500],...
  %          'projection','plain','fontsize',8,'PHI2',phi2_1,'Bunge','grid','minmax','off'); 
    %annotate(qref,'fontsize',10,'marker','s','MarkerFaceColor','red');
  %  colorbar;
  %  savefigure([fname '.phi2.0deg.ortho.png']);
  %  close;

  %  figure;
  %  plotodf(odf_mmm,'sections',1,'gray','contourf',[1 2 5 10 25 50 100 200 300 500],...
  %        'projection','plain','fontsize',8,'PHI2',phi2_2,'Bunge','grid','minmax','off'); 
    %annotate(qref,'fontsize',10,'gray','marker','s','MarkerFaceColor','red');
  %  colorbar;
  %  savefigure([fname '.phi2.45deg.ortho.png']);
  %  close;

  %  figure;
  %  plotodf(odf_mmm,'sections',1,'gray','contourf',[1 2 5 10 25 50 100 200 300 500],...
  %         'projection','plain','fontsize',8,'PHI2',phi2_3,'Bunge','grid','minmax','off'); 
    %annotate(qref,'fontsize',10,'gray','marker','s','MarkerFaceColor','red');
  %  colorbar;
  %  savefigure([fname '.phi2.65deg.ortho.png']);
  %  close;
%end;

%% calculate maximum intensity and orientation at this maximum value

[MAXIMUMINTENSITY, MAXORI] = max(odf);
MAXIMUMINTENSITY
MAXORI

%% get the texture index
TextureIndex = textureindex(odf);  
TextureIndex

m111 = vector3d(1,1,1,CS);
figure
plotpdf(odf,m111);
colorbar;
savefigure ([fname '.ipdf.tricline.png']);
close;

%compare to random ipdf
uniformodf = uniformODF(CS,SS);
figure 
plotpdf(uniformodf,m111);
colorbar;
savefigure ([fname '.random.ipdf.tricline.png']);
close;

mackenzie_rnd = calcAngleDistribution(uniformodf);
degg = zeros(300,1);
for k=1:300
    degg(k) = k/300 * 62.8;
end;
mackenzie_shoemake = calcAngleDistribution(odf);

figure;
plot(degg,mackenzie_rnd,'color',[1 0.5 0]);
hold on
plot(degg,mackenzie_shoemake,'color',[0 0 0]);
savefigure([fname 'orangeMacKenzieUniRandom_blackOriLibRandomShoemake.png']);

close;

end;

%% calculate volume fractions you are interested
VolCubeComp = 0.0 + volume(odf,Cube,optVolRes);
    VolCubeComp = VolCubeComp + volume(odf,CubeND22,optVolRes);
    VolCubeComp = VolCubeComp + volume(odf,CubeND22b,optVolRes);
    VolCubeComp = VolCubeComp + volume(odf,CubeRD22,optVolRes);
    VolCubeComp = VolCubeComp + volume(odf,CubeRD22b,optVolRes);
    VolCubeComp = VolCubeComp + volume(odf,CubeTD22,optVolRes);
    VolCubeComp = VolCubeComp + volume(odf,CubeTD22b,optVolRes);
    VolCubeComp = VolCubeComp + volume(odf,CubeND45,optVolRes);

VolPQComp = 0.0 + volume(odf,PLage,optVolRes);
    VolPQComp = VolPQComp + volume(odf,PLage2,optVolRes);
    VolPQComp = VolPQComp + volume(odf,QLage,optVolRes);
    VolPQComp = VolPQComp + volume(odf,QLage2,optVolRes);
    VolPQComp = VolPQComp + volume(odf,QLage3,optVolRes);
    VolPQComp = VolPQComp + volume(odf,QLage4,optVolRes);

VolCopper = volume(odf,Copper,optVolRes) + volume(odf,Copper2,optVolRes);
VolDillamore = volume(odf,Dill,optVolRes) + volume(odf,Dill2,optVolRes);

VolGoss = volume(odf,Goss,optVolRes);


VolBrass = volume(odf,Brass,optVolRes) + volume(odf,Brass2,optVolRes);

VolS = 0.0 + volume(odf,S,optVolRes);
    VolS = VolS + volume(odf,S2,optVolRes);
    VolS = VolS + volume(odf,S3,optVolRes);
    VolS = VolS + volume(odf,S4,optVolRes);


%% Write Volume fraction comparison

fid = fopen ( [fname '.Protocol.txt'], 'w+' );
fprintf(fid, 'Sample Name: %s\n', fname);
fprintf(fid, 'Crystal Symmetry: %s\n', Laue(CS));
fprintf(fid, 'Sample Symmetry: %s\n', Laue(SS));
fprintf(fid, 'Kernel: %s, halfwidth/°%.2f\n', get(psi,'name'), get(psi,'hw'));
% Volume Fractions
%if (length(SS) == 1)
fprintf(fid, '---------------initial ODF vs. discretized from EBSD---------------\n');
fprintf(fid, 'orthorhombic volume fractions( %.2f°) \n\n', optVolRes/degree);
% fprintf(fid, 'Texture Type: %s\n', Texturetype);
fprintf(fid, 'VolBrass     \t%.3f\n', VolBrass);
fprintf(fid, 'VolCopper    \t%.3f\n', VolCopper);
fprintf(fid, 'VolCube      \t%.3f\n', VolCubeComp);
fprintf(fid, 'VolDillamore \t%.3f\n', VolDillamore);
fprintf(fid, 'VolGoss      \t%.3f\n', VolGoss);
fprintf(fid, 'VolPQ         \t%.3f\n', VolPQComp);
fprintf(fid, 'VolS         \t%.3f\n', VolS);
fprintf(fid, 'TextureIndex \t%.3f \n', TextureIndex);
%end;
%if (length(SS) == 4)
%    fprintf(fid, '---------------initial ODF vs. discretized from EBSD---------------\n');
%fprintf(fid, 'orthorhombic Volume Fractions volumesphere( %.2f°) \n\n', optVolRes/degree);
% fprintf(fid, 'Texture Type: %s\n', Texturetype);
%fprintf(fid, 'Cube     \t%.3f\t\t%.3f\n', InitVolCube,DiscrVolCube);
%fprintf(fid, 'CubeTwin \t%.3f\t\t%.3f\n', InitVolCubeTwin,DiscrVolCubeTwin);
%fprintf(fid, 'CubeND22 \t%.3f\t\t%.3f\n', InitVolCubeND22,DiscrVolCubeND22);
%fprintf(fid, 'CubeRD22 \t%.3f\t\t%.3f\n', InitVolCubeRD22,DiscrVolCubeRD22);
%fprintf(fid, 'CubeTD22 \t%.3f\t\t%.3f\n', InitVolCubeTD22,DiscrVolCubeTD22);
%fprintf(fid, 'CubeND45 \t%.3f\t\t%.3f\n', InitVolCubeND45,DiscrVolCubeND45);
%%fprintf(fid, 'Goss    \t%.3f\t\t%.3f\n', InitVolGoss,DiscrVolGoss);
%%fprintf(fid, 'invGoss \t%.3f\t\t%.3f\n', InitVolinvGoss,DiscrVolinvGoss);
%fprintf(fid, 'Copper   \t%.3f\t\t%.3f\n', InitVolCopper, DiscrVolCopper);
%fprintf(fid, 'CTwin    \t%.3f\t\t%.3f\n', InitVolCTwin, DiscrVolCTwin);
%fprintf(fid, 'CGIA     \t%.3f\t\t%.3f\n', InitVolCopperGIA,DiscrVolCopperGIA);
%fprintf(fid, 'S       \t%.3f\t\t%.3f\n', InitVolS5,DiscrVolS5);
%fprintf(fid, 'SHirsch   \t%.3f\t\t%.3f\n', InitVolSHirsch, DiscrVolSHirsch);
%fprintf(fid, 'R       \t%.3f\t\t%.3f\n', InitVolR5, DiscrVolR5);
%fprintf(fid, 'Brass   \t%.3f\t\t%.3f\n', InitVolBrass, DiscrVolBrass);
%fprintf(fid, 'PLage   \t%.3f\t\t%.3f\n', InitVolPLage, DiscrVolPLage);
%fprintf(fid, 'QLage   \t%.3f\t\t%.3f\n', InitVolQLage, DiscrVolQLage);
%fprintf(fid, 'Dillamore \t%.3f\t\t%.3f\n', InitVolDill, DiscrVolDill);
%fprintf(fid, '\nTI\t\t\t%.3f\t\t%.3f\n', InitialTI, DiscrTI);
%end;
fclose (fid);    

    
%end of automatization
%end

%% profiling ended, save profiling information
if (optProfile == 1)
    profile viewer
    p = profile(fname);
    %profsave(p,'profile_results')
    profile off
end;
