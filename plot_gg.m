function plot_gg(points,anchors,set_loc,set_unloc,ConnectivityM_GG,BoxScale,graphlabel)
points = points*BoxScale;
npoints = size(points,1);
p_anchors = points(anchors,:);
p_loc = points(set_loc,:);
p_unloc = points(set_unloc,:);
idx1 = 1; idx2 = 2;
markersize = 10;
r1 = 1; r2 = 2;
h = figure();
set(h,'name',graphlabel,'Numbertitle','off')
%% plot anchor triangle
hold on;
plot([p_anchors(1,1),p_anchors(2,1)],[p_anchors(1,2),p_anchors(2,2)],'k','linewidth',2);
hold on;
plot([p_anchors(1,1),p_anchors(3,1)],[p_anchors(1,2),p_anchors(3,2)],'k','linewidth',2);
hold on;
plot([p_anchors(3,1),p_anchors(2,1)],[p_anchors(3,2),p_anchors(2,2)],'k','linewidth',2);
%% plot edges of generated graph
hold on;
for j = 4:npoints
    if ~ismember(j,set_loc)
        %continue
    end
    nei = find(ConnectivityM_GG(j,:)==1);
    if ~isempty(nei) 
        len = length(nei);
        starts = repmat(points(j,:),len,1);
        ends = points(nei,:);
        hold on;
        plot([starts(:,r1) ends(:,r1)]',[starts(:,r2) ends(:,r2)]','k','linewidth',2);
    end
end
%% plot anchors
hold on;
p1 = plot(p_anchors(:,idx1),p_anchors(:,idx2),'d','markersize',markersize-2,...
    'color',[214,39,40]/255,'linewidth',2,'markerfacecolor',[214,39,40]/255);
%% plot localizable nodes 
hold on;
p2 = plot(p_loc(:,idx1),p_loc(:,idx2),'o','markersize',markersize,'markerfacecolor',[255,255,255]/255,...
    'color',[44,160,44]/255,'linewidth',2);
hold on;

%% plot nodes outside anchor triangle 
hold on;
p3 = plot(p_unloc(:,idx1),p_unloc(:,idx2),'s','markersize',markersize,'markerfacecolor',[255,255,255]/255,...
    'color',[31,119,180]/255,'linewidth',1.8);%[44,160,44]/255
%% plot node labels
hold on;
for i = 1:npoints
    text(points(i,r1)+3.5,points(i,r2)-1, num2str(i),'FontSize',12);
end
%% set legends
% l = legend([p1,p2,p3],'$\mathcal{A}$','$\mathcal{S}^{*}$',...
%     '$\mathcal{S}\setminus\mathcal{S}^{*}$',...
%     'NumColumnsMode','manual','NumColumns',1,'Location','SouthEast','FontSize',12);
% set(l,'Interpreter','latex');
xlim([-55 55]);
ylim([-55 55]);
x_Matrix = -50:25:50;
set(gca,'xtick', x_Matrix)
set(gca,'ytick', x_Matrix)

grid on
box on
axis('square');
%centerX = 500;centerY = 500;
%width = 300;
%height = 300;
%set(gcf,'position',[centerX,centerY,width,height])
end



