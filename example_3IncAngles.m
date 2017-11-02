% Script to create the plots TB H-pol for three incidence angles e.g. 20, 30, 40°
% Region min to max coordinates [47.5 7.5 50 10]: Case study Neckar Catchment.
% This script has been used to produce Figure 2 for JORS paper.
% Copiright 2017, Pablo Saavedra G. (see LICENCE.TXT)
%

% Geographical limits:
lat_lon_box = [47.5 7.5 50 10];
latlim = lat_lon_box([1,3]);
lonlim = lat_lon_box([2,4]);
% selecting DBL file to process and plot:
file_example = ['/tmp/smos/' ...
                  'SM_OPER_MIR_SCLF1C_20150702T042618_20150702T051937_620_001_1.DBL'];
if exist(file_example,'file'),
  %TSI = Process_SMOSxL1C(file_example,lat_lon_box);  % no MAT-file output!
else
  % Open file browser GUI:
  %TSI = Process_SMOSxL1C;
end

idx_pol =1; % 1 for H, 2 for V
inc_ang = [35 40 45]; % indices for the Incidence Angles
N_ang = length(inc_ang);
figure;
set(gcf,'PaperPositionMode','auto','Position',[117 500 1326 395]);
for idx=1:N_ang,
  ax(idx) = subplot(1,N_ang,idx);
  scatter(TSF.GridPoint_Longitude,TSF.GridPoint_Latitude,150,...
	  TSF.TB_Fixed_IncAngle(:,inc_ang(idx)+1,idx_pol),'o','filled');
  tl=title(sprintf('TB_H [K] @ \\theta_{inc} = %2d^\\circ',...
			   inc_ang(idx)),...
	   'FontName','utopia','FontSize',15,'render','latex');
  tmp=get(ax(idx),'XTickLabel');
  set(ax(idx),'XTickLabel',[tmp repmat(' °E',size(tmp,1),1)]);
  tmp=get(ax(idx),'YTickLabel');
  set(ax(idx),'YTickLabel',[tmp repmat(' °N',size(tmp,1),1)]);
  xlabel('LONGITUD');
  if idx==1, ylabel('LATITUDE'); end
end
set(ax,'FontName','utopia','FontSize',19,'LineWidth',1,'Color','none',...
    'Box','on','GridLineStyle','--','XLim',lonlim,'YLim',latlim,...
    'TickDir','out','XMinorTick','on','YMinorTick','on','Layer','top',...
    'LineWidth',1.5,'TickLength',[0.02 0.05],'CLim',[205 285]);

colormap(jet);
hb = colorbar('FontSize',15,'FontName','utopia');;
return;

ax = worldmap([-90 90],[-90 100]); %([0 80],[-35 45]);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.5 0.5 0.5]);
pp=scatterm(TSI.GridPoint_Latitude([1:15:end]), ...
            TSI.GridPoint_Longitude([1:15:end]),5,TBxy([1:15:end]),'o','filled');
caxis([220 290]);
title('L1C SMOS TB X-pol [K]','FontSize',20);
colorbar('EastOutSide');
setm(ax,'FontSize',19);

figure(2);
set(gcf,'PaperPositionMode','auto');
bx = worldmap('germany');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(bx, land, 'FaceColor', [0.5 0.5 0.5]);
pp=scatterm(TSI.GridPoint_Latitude([1:5:end]), ...
            TSI.GridPoint_Longitude([1:5:end]),30,TBxy([1:5:end]),'o','filled');
caxis([220 290]);
title('SMOS TB X-pol [K]','FontSize',20);
colorbar('EastOutSide');
setm(bx,'FontSize',19);


% end of script
