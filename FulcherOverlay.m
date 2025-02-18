%Generates subfigures 5a, b, c, g, h, i

diag = struct();
diag.Rp = [1.9, 1.9];
diag.Zp = [-2, 0];
diag.Lp = [0,0];
diag.Dp = [0.05, 0.05];
diag.ap = [-pi/2, -pi/2];

fpath = '/common/projects/physics/SOLPS/dmoulton/mastu_sxd_46860_450ms/puff=13.0e21_pump=0.001_divchiconstfix_bcmom2_visper1_rxf0p1_parmvsa2_redoutpfrtrans_pufflfs_pin1.0';
%fpath = '/common/projects/physics/SOLPS/dmoulton/mastu_sxd_46860_450ms/puff=0.3e21_pump=0.001_nodrifts_bcmom2_parmvsa2_redoutpfrtrans';
output_SXD = SOLPS_ExtractParams_ForOpticalModel_balance([fpath, '/balance.nc'],[fpath,'/fort.44'],[fpath,'/fort.46'],[fpath,'/fort.33'],[fpath,'/fort.34']);
output_SXD = SOLPS_Balmer_Model(output_SXD);
output_SXD = SOLPS_SynDiag_Baratron2(output_SXD,diag,1);

fpath = '/common/projects/physics/SOLPS/dmoulton/mastu_ed_47079_540ms/puff=13.0e21_pump=0.001_divchiconstfix_bcmom2_visper1_rxf0p1_parmvsa2_redoutpfrtrans_pufflfs_pin1.0';
output_ED = SOLPS_ExtractParams_ForOpticalModel_balance([fpath, '/balance.nc'],[fpath,'/fort.44'],[fpath,'/fort.46'],[fpath,'/fort.33'],[fpath,'/fort.34']);
output_ED = SOLPS_Balmer_Model(output_ED);
output_ED = SOLPS_SynDiag_Baratron2(output_ED,diag,1);

fpath = '/common/projects/physics/SOLPS/dmoulton/mastu_cd_46866_450ms/puff=13.0e21_pump=0.001_divchiconstfix_bcmom2_visper1_rxf0p1_parmvsa2_redoutpfrtrans_pufflfs_pin1.0';
output_CD = SOLPS_ExtractParams_ForOpticalModel_balance([fpath, '/balance.nc'],[fpath,'/fort.44'],[fpath,'/fort.46'],[fpath,'/fort.33'],[fpath,'/fort.34']);
output_CD = SOLPS_Balmer_Model(output_CD);
output_CD = SOLPS_SynDiag_Baratron2(output_CD,diag,1);

rwall = [1.6,1.2503,1.3483,1.47,1.47,1.45,1.45,1.3214,1.1904,0.89296,0.86938,0.83981,0.82229,0.81974,0.81974,0.82734,0.8548,0.89017,0.91974,0.94066,1.555,1.85,2,2,2,2,1.3188,1.7689,1.7301,1.35,1.09,1.09,0.90576,0.90717,0.53948,0.53594,0.5074,0.5074,0.4788,0.4788,0.333,0.333,0.334,0.261,0.261,0.261,0.261,0.261,0.261];
zwall = [1,1,0.86,0.86,0.81,0.81,0.82,0.82,1.007,1.304,1.3312,1.3826,1.4451,1.4812,1.4936,1.5318,1.5696,1.5891,1.5936,1.5936,1.567,1.08,1.08,1.7,2.035,2.169,2.169,1.7189,1.68,2.06,2.06,2.06,1.8786,1.8772,1.5095,1.5017,1.4738,1.4738,1.4458,1.4458,1.303,1.1,1.1,0.502,0.348,0.348,0.146,0.146,0];

output_CD.Baratron.Pressure.P.*1.602e-19
output_ED.Baratron.Pressure.P.*1.602e-19
output_SXD.Baratron.Pressure.P.*1.602e-19

save('SOLPS_shape_scan','output_CD','output_ED','output_SXD')

%make plot

figure
colormap hot
pcolor(output_SXD.params.cr, output_SXD.params.cz, output_SXD.AMJUEL.Emiss.FulcherEmiss_au)
shading interp
hold on
TeNorm = 5;
plot(output_SXD.params.cr(:,output_SXD.params.sep),output_SXD.params.cz(:,output_SXD.params.sep),'LineWidth',4,'Color',[0,0,1,0.7])
contour(output_SXD.params.cr, output_SXD.params.cz, output_SXD.params.te,[TeNorm,TeNorm],'b--','LineWidth',3)
plot(rwall,-zwall,'k')
%plot([0.95,0.95],[-2.1,-1.4],'Color',[1,0,1,0.5],'LineWidth',3)
axis([0.6,1.8,-2.1,-1.5])
caxis([0 5e19])
pbaspect([(1.8-0.6) (2.1-1.5) 1])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
export_fig('-dpng','-r300','SXD_Fulcher.png')

figure
colormap hot
pcolor(output_ED.params.cr, output_ED.params.cz, output_ED.AMJUEL.Emiss.FulcherEmiss_au)
shading interp
hold on
TeNorm = 5;
plot(output_ED.params.cr(:,output_ED.params.sep),output_ED.params.cz(:,output_ED.params.sep),'LineWidth',4,'Color',[0,1,0,0.5])
contour(output_ED.params.cr, output_ED.params.cz, output_ED.params.te,[TeNorm,TeNorm],'g--','LineWidth',3,'EdgeAlpha',0.7)
%contour(output_SXD.params.cr, output_SXD.params.cz, output_SXD.params.te,[TeNorm,TeNorm],'b--','LineWidth',3,'EdgeAlpha',0.7)
plot(rwall,-zwall,'k')
%plot([0.95,0.95],[-2.1,-1.4],'Color',[1,0,1,0.5],'LineWidth',3)
axis([0.6,1.8,-2.1,-1.5])
caxis([0 5e19])
pbaspect([(1.8-0.6) (2.1-1.5) 1])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
export_fig('-dpng','-r300','ED_Fulcher.png')

figure
colormap hot
pcolor(output_CD.params.cr, output_CD.params.cz, output_CD.AMJUEL.Emiss.FulcherEmiss_au)
shading interp
hold on
TeNorm = 5;
plot(output_CD.params.cr(:,output_CD.params.sep),output_CD.params.cz(:,output_CD.params.sep),'LineWidth',4,'Color',[1,0,0,0.6])
contour(output_CD.params.cr, output_CD.params.cz, output_CD.params.te,[TeNorm,TeNorm],'r--','LineWidth',3,'EdgeAlpha',0.7)
plot(rwall,-zwall,'k')
%plot([0.95,0.95],[-2.1,-1.4],'Color',[1,0,1,0.5],'LineWidth',3)
axis([0.6,1.8,-2.1,-1.5])
caxis([0 5e19])
pbaspect([(1.8-0.6) (2.1-1.5) 1])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
export_fig('-dpng','-r300','CD_Fulcher.png')

% plot Fulcher = SXD

load('MAST_U vsl.mat')
crG = output_SXD.params.crG;
czG = output_SXD.params.czG;
A = output_SXD.ADAS.Hydro.Reaction.IonMap;

% Assuming crG, czG, and A are of size r x z x 4
[r, z, ~] = size(crG); % Get dimensions of the grid
% Initialize a figure
figure;
colormap hot
hold on;
% Loop over each element in the grid
for i = 1:r-1
for j = 1:z-1
% Get the coordinates of the 4 vertices of the polygon
x_coords = squeeze(crG(i,j,:)); % x-coordinates of the 4 vertices
y_coords = squeeze(czG(i,j,:)); % y-coordinates of the 4 vertices
B = boundary(x_coords,y_coords);
x_coords=x_coords(B);
y_coords=y_coords(B);
% Get the value for A in the current polygon (this could be an average or a specific one)
A_value = A(i,j); % Assuming you want to plot the value corresponding to the center of the element
% Create a polygon using the coordinates
patch(x_coords, y_coords, A_value, 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 1);
end
end
% Add color bar to show the values of A
colorbar;
caxis([0,2e23])
axis equal;
axis([0.3,1.75,-2.1,-0.45])
%plot([vsl(1,:);vsl(3,:);],[vsl(2,:);vsl(4,:);],'k','LineWidth',2)
%plot([0.3,1.75],[-1.6,-1.6],'r','LineWidth',2)
%plot([0.3,1.75],[-1.07,-1.07],'m','LineWidth',2)
hold off;
set(gcf, 'Renderer', 'Painters')
print2eps IonSource_SXD_SOLPS

figure
plot([vsl(1,:);vsl(3,:);],[vsl(2,:);vsl(4,:);],'k','LineWidth',2)
hold on
plot([0.3,1.75],[-1.6,-1.6],'r','LineWidth',2)
plot([0.3,1.75],[-1.07,-1.07],'m','LineWidth',2)
axis equal
axis([0.3,1.75,-2.1,-0.45])
print2eps SXD_SOLPS_vessel

A = output_SXD.AMJUEL.Emiss.FulcherEmiss_au;

% Assuming crG, czG, and A are of size r x z x 4
[r, z, ~] = size(crG); % Get dimensions of the grid
% Initialize a figure
figure;
colormap hot
hold on;
% Loop over each element in the grid
for i = 1:r-1
for j = 1:z-1
% Get the coordinates of the 4 vertices of the polygon
x_coords = squeeze(crG(i,j,:)); % x-coordinates of the 4 vertices
y_coords = squeeze(czG(i,j,:)); % y-coordinates of the 4 vertices
B = boundary(x_coords,y_coords);
x_coords=x_coords(B);
y_coords=y_coords(B);
% Get the value for A in the current polygon (this could be an average or a specific one)
A_value = A(i,j); % Assuming you want to plot the value corresponding to the center of the element
% Create a polygon using the coordinates
patch(x_coords, y_coords, A_value, 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 1);
end
end
% Add color bar to show the values of A
colorbar;
caxis([0,1e20])
axis equal;
axis([0.3,1.75,-2.1,-0.45])
plot([vsl(1,:);vsl(3,:);],[vsl(2,:);vsl(4,:);],'k','LineWidth',2)
hold off;
set(gcf, 'Renderer', 'Painters')
print2eps Fulcher_SXD_SOLPS_XPI

% plot Fulcher = SXD

load('MAST_U vsl.mat')
crG = output_CD.params.crG;
czG = output_CD.params.czG;
A = output_CD.ADAS.Hydro.Reaction.IonMap;

% Assuming crG, czG, and A are of size r x z x 4
[r, z, ~] = size(crG); % Get dimensions of the grid
% Initialize a figure
figure;
colormap hot
hold on;
% Loop over each element in the grid
for i = 1:r-1
for j = 1:z-1
% Get the coordinates of the 4 vertices of the polygon
x_coords = squeeze(crG(i,j,:)); % x-coordinates of the 4 vertices
y_coords = squeeze(czG(i,j,:)); % y-coordinates of the 4 vertices
B = boundary(x_coords,y_coords);
x_coords=x_coords(B);
y_coords=y_coords(B);
% Get the value for A in the current polygon (this could be an average or a specific one)
A_value = A(i,j); % Assuming you want to plot the value corresponding to the center of the element
% Create a polygon using the coordinates
patch(x_coords, y_coords, A_value, 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 1);
end
end
% Add color bar to show the values of A
colorbar;
caxis([0,2e23])
axis equal;
axis([0.3,1.75,-2.1,-0.45])
%plot([vsl(1,:);vsl(3,:);],[vsl(2,:);vsl(4,:);],'k','LineWidth',2)
%plot([0.3,1.75],[-1.6,-1.6],'r','LineWidth',2)
%plot([0.3,1.75],[-1.07,-1.07],'m','LineWidth',2)
hold off;
set(gcf, 'Renderer', 'Painters')
print2eps IonSource_CD_SOLPS

A = output_CD.AMJUEL.Emiss.FulcherEmiss_au;

% Assuming crG, czG, and A are of size r x z x 4
[r, z, ~] = size(crG); % Get dimensions of the grid
% Initialize a figure
figure;
colormap hot
hold on;
% Loop over each element in the grid
for i = 1:r-1
for j = 1:z-1
% Get the coordinates of the 4 vertices of the polygon
x_coords = squeeze(crG(i,j,:)); % x-coordinates of the 4 vertices
y_coords = squeeze(czG(i,j,:)); % y-coordinates of the 4 vertices
B = boundary(x_coords,y_coords);
x_coords=x_coords(B);
y_coords=y_coords(B);
% Get the value for A in the current polygon (this could be an average or a specific one)
A_value = A(i,j); % Assuming you want to plot the value corresponding to the center of the element
% Create a polygon using the coordinates
patch(x_coords, y_coords, A_value, 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 1);
end
end
% Add color bar to show the values of A
colorbar;
caxis([0,1e20])
axis equal;
axis([0.3,1.75,-2.1,-0.45])
%plot([vsl(1,:);vsl(3,:);],[vsl(2,:);vsl(4,:);],'k','LineWidth',2)
hold off;
set(gcf, 'Renderer', 'Painters')
print2eps Fulcher_CD_SOLPS_XPI

% plot Fulcher = ED

load('MAST_U vsl.mat')
crG = output_ED.params.crG;
czG = output_ED.params.czG;
A = output_ED.ADAS.Hydro.Reaction.IonMap;

% Assuming crG, czG, and A are of size r x z x 4
[r, z, ~] = size(crG); % Get dimensions of the grid
% Initialize a figure
figure;
colormap hot
hold on;
% Loop over each element in the grid
for i = 1:r-1
for j = 1:z-1
% Get the coordinates of the 4 vertices of the polygon
x_coords = squeeze(crG(i,j,:)); % x-coordinates of the 4 vertices
y_coords = squeeze(czG(i,j,:)); % y-coordinates of the 4 vertices
B = boundary(x_coords,y_coords);
x_coords=x_coords(B);
y_coords=y_coords(B);
% Get the value for A in the current polygon (this could be an average or a specific one)
A_value = A(i,j); % Assuming you want to plot the value corresponding to the center of the element
% Create a polygon using the coordinates
patch(x_coords, y_coords, A_value, 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 1);
end
end
% Add color bar to show the values of A
colorbar;
caxis([0,2e23])
axis equal;
axis([0.3,1.75,-2.1,-0.45])
%plot([vsl(1,:);vsl(3,:);],[vsl(2,:);vsl(4,:);],'k','LineWidth',2)
%plot([0.3,1.75],[-1.6,-1.6],'r','LineWidth',2)
%plot([0.3,1.75],[-1.07,-1.07],'m','LineWidth',2)
hold off;
set(gcf, 'Renderer', 'Painters')
print2eps IonSource_ED_SOLPS

A = output_ED.AMJUEL.Emiss.FulcherEmiss_au;

% Assuming crG, czG, and A are of size r x z x 4
[r, z, ~] = size(crG); % Get dimensions of the grid
% Initialize a figure
figure;
colormap hot
hold on;
% Loop over each element in the grid
for i = 1:r-1
for j = 1:z-1
% Get the coordinates of the 4 vertices of the polygon
x_coords = squeeze(crG(i,j,:)); % x-coordinates of the 4 vertices
y_coords = squeeze(czG(i,j,:)); % y-coordinates of the 4 vertices
B = boundary(x_coords,y_coords);
x_coords=x_coords(B);
y_coords=y_coords(B);
% Get the value for A in the current polygon (this could be an average or a specific one)
A_value = A(i,j); % Assuming you want to plot the value corresponding to the center of the element
% Create a polygon using the coordinates
patch(x_coords, y_coords, A_value, 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 1);
end
end
% Add color bar to show the values of A
colorbar;
caxis([0,1e20])
axis equal;
axis([0.3,1.75,-2.1,-0.45])
%plot([vsl(1,:);vsl(3,:);],[vsl(2,:);vsl(4,:);],'k','LineWidth',2)
hold off;
set(gcf, 'Renderer', 'Painters')
print2eps Fulcher_ED_SOLPS_XPI

%ionisation fractions

Ion_CD = sum(sum(output_CD.ADAS.Hydro.Reaction.IonMap(output_CD.params.omp:end,:).*output_CD.params.vol(output_CD.params.omp:end,:)));
Ion_ED = sum(sum(output_ED.ADAS.Hydro.Reaction.IonMap(output_ED.params.omp:end,:).*output_ED.params.vol(output_ED.params.omp:end,:)));
Ion_SXD = sum(sum(output_SXD.ADAS.Hydro.Reaction.IonMap(output_SXD.params.omp:end,:).*output_SXD.params.vol(output_SXD.params.omp:end,:)));
Ion_CD_U = sum(sum(output_CD.ADAS.Hydro.Reaction.IonMap(output_CD.params.omp:end,:).*output_CD.params.vol(output_CD.params.omp:end,:).*(output_CD.params.cz(output_CD.params.omp:end,:)>-1.3)));
Ion_ED_U = sum(sum(output_ED.ADAS.Hydro.Reaction.IonMap(output_ED.params.omp:end,:).*output_ED.params.vol(output_ED.params.omp:end,:).*(output_ED.params.cz(output_ED.params.omp:end,:)>-1.3)));
Ion_SXD_U = sum(sum(output_SXD.ADAS.Hydro.Reaction.IonMap(output_SXD.params.omp:end,:).*output_SXD.params.vol(output_SXD.params.omp:end,:).*(output_SXD.params.cz(output_SXD.params.omp:end,:)>-1.3)));

%calculate core ionisation
Ion_CD_core = sum(sum(output_CD.ADAS.Hydro.Reaction.IonMap(output_CD.params.omp:output_CD.params.xcut(end)-1,1:output_CD.params.sep).*output_CD.params.vol(output_CD.params.omp:output_CD.params.xcut(end)-1,1:output_CD.params.sep)));
Ion_ED_core = sum(sum(output_ED.ADAS.Hydro.Reaction.IonMap(output_ED.params.omp:output_ED.params.xcut(end)-1,1:output_ED.params.sep).*output_ED.params.vol(output_ED.params.omp:output_ED.params.xcut(end)-1,1:output_ED.params.sep)));
Ion_SXD_core = sum(sum(output_SXD.ADAS.Hydro.Reaction.IonMap(output_SXD.params.omp:output_SXD.params.xcut(end)-1,1:output_SXD.params.sep).*output_SXD.params.vol(output_SXD.params.omp:output_SXD.params.xcut(end)-1,1:output_SXD.params.sep)));

(Ion_SXD_U-Ion_SXD_core+output_SXD.params.upstream_ionflux)/(Ion_SXD-Ion_SXD_core+output_SXD.params.upstream_ionflux)
(Ion_ED_U-Ion_ED_core+output_ED.params.upstream_ionflux)/(Ion_ED-Ion_ED_core+output_ED.params.upstream_ionflux)
(Ion_CD_U-Ion_CD_core+output_CD.params.upstream_ionflux)/(Ion_CD-Ion_CD_core+output_CD.params.upstream_ionflux)

Ion_CD_U = sum(sum(output_CD.ADAS.Hydro.Reaction.IonMap(output_CD.params.omp:output_CD.params.xcut(end)-1,output_CD.params.sep:end).*output_CD.params.vol(output_CD.params.omp:output_CD.params.xcut(end)-1,output_CD.params.sep:end).*(output_CD.params.cz(output_CD.params.omp:output_CD.params.xcut(end)-1,output_CD.params.sep:end)>-1.3)));
Ion_ED_U = sum(sum(output_ED.ADAS.Hydro.Reaction.IonMap(output_ED.params.omp:output_ED.params.xcut(end)-1,output_ED.params.sep:end).*output_ED.params.vol(output_ED.params.omp:output_ED.params.xcut(end)-1,output_ED.params.sep:end).*(output_ED.params.cz(output_ED.params.omp:output_ED.params.xcut(end)-1,output_ED.params.sep:end)>-1.3)));
Ion_SXD_U = sum(sum(output_SXD.ADAS.Hydro.Reaction.IonMap(output_SXD.params.omp:output_SXD.params.xcut(end)-1,output_SXD.params.sep:end).*output_SXD.params.vol(output_SXD.params.omp:output_SXD.params.xcut(end)-1,output_SXD.params.sep:end).*(output_SXD.params.cz(output_SXD.params.omp:output_SXD.params.xcut(end)-1,output_SXD.params.sep:end)>-1.3)));

(Ion_SXD_U+output_SXD.params.upstream_ionflux)/(Ion_SXD-Ion_SXD_core+output_SXD.params.upstream_ionflux)
(Ion_ED_U+output_ED.params.upstream_ionflux)/(Ion_ED-Ion_ED_core+output_ED.params.upstream_ionflux)
(Ion_CD_U+output_CD.params.upstream_ionflux)/(Ion_CD-Ion_CD_core+output_CD.params.upstream_ionflux)

%obtain ionisation fraction covered by DMS

L=2;
angadj=0;

load([matlab_home,'/diagnostic_files/MastU_York_V1_MU01.mat'])
%get cells captured by DMS chords
SpecView_V1=[DCD.RD(1), DCD.ZD(1); DCD.RD(1) - L*cos(DCD.PVD(1)+angadj), DCD.ZD(1) + L*sin(DCD.PVD(1)+angadj); DCD.RD(10) - L*cos(DCD.PVD(10)), DCD.ZD(10) + L*sin(DCD.PVD(10)); DCD.RD(10), DCD.ZD(10)];
ISpecView_V1_SXD = logical(inpolygon(output_SXD.params.cr,output_SXD.params.cz,SpecView_V1(:,1),SpecView_V1(:,2)));
ISpecView_V1_ED = logical(inpolygon(output_ED.params.cr,output_ED.params.cz,SpecView_V1(:,1),SpecView_V1(:,2)));
ISpecView_V1_CD = logical(inpolygon(output_CD.params.cr,output_CD.params.cz,SpecView_V1(:,1),SpecView_V1(:,2)));

load([matlab_home,'/diagnostic_files/MastU_York_V2_MU01.mat'])
%get cells captured by DMS chords
SpecView_V2=[DCD.RD(1), DCD.ZD(1); DCD.RD(1) - L*cos(DCD.PVD(1)+angadj), DCD.ZD(1) + L*sin(DCD.PVD(1)+angadj); DCD.RD(10) - L*cos(DCD.PVD(10)), DCD.ZD(10) + L*sin(DCD.PVD(10)); DCD.RD(10), DCD.ZD(10)];
ISpecView_V2_SXD = logical(inpolygon(output_SXD.params.cr,output_SXD.params.cz,SpecView_V2(:,1),SpecView_V2(:,2)));
ISpecView_V2_ED = logical(inpolygon(output_ED.params.cr,output_ED.params.cz,SpecView_V2(:,1),SpecView_V2(:,2)));
ISpecView_V2_CD = logical(inpolygon(output_CD.params.cr,output_CD.params.cz,SpecView_V2(:,1),SpecView_V2(:,2)));

ISpecView_SXD = logical(logical(ISpecView_V1_SXD) + logical(ISpecView_V2_SXD));
ISpecView_ED = logical(logical(ISpecView_V1_ED)+logical(ISpecView_V2_SXD));
ISpecView_CD = logical(logical(ISpecView_V1_CD)+logical(ISpecView_V2_CD));

Ion_CD_U = sum(sum(output_CD.ADAS.Hydro.Reaction.IonMap(ISpecView_CD).*output_CD.params.vol(ISpecView_CD)));
Ion_ED_U = sum(sum(output_ED.ADAS.Hydro.Reaction.IonMap(ISpecView_ED).*output_ED.params.vol(ISpecView_ED)));
Ion_SXD_U = sum(sum(output_SXD.ADAS.Hydro.Reaction.IonMap(ISpecView_SXD).*output_SXD.params.vol(ISpecView_SXD)));

(Ion_SXD_U)/(Ion_SXD-Ion_SXD_core+output_SXD.params.upstream_ionflux)
(Ion_ED_U)/(Ion_ED-Ion_ED_core+output_ED.params.upstream_ionflux)
(Ion_CD_U)/(Ion_CD-Ion_CD_core+output_CD.params.upstream_ionflux)

disp('...')

