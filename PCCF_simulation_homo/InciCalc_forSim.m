function [SumPos, CPos, MPos, Nt_ori, Nt_mod, Ncluser, Nmono] = InciCalc(Pos, PixSize, InciTh, FigDraw)
% Convert threshold from nm to pix
Th = InciTh/PixSize;
% Generate tree diagram
Z = linkage(Pos(:,1:2),'single');

% Draw tree diagram
if FigDraw == 1
    figure(1)
    dendrogram(Z,0,'ColorThreshold',Th);
    ylabel('Distance (pix)');
    % Delete extra space
    ax = gca;
    ax.XTickLabel = []; 
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end

% Clacified to clusters 
T = cluster(Z,'Cutoff',Th, 'Criterion','distance');
% Find # of cluster
Nspots = max(T);
% Average positions of clusters
L=[]; L2=[]; Ncluser=1; Nmono=1; CPos=[]; MPos=[]; CPosOri{1}=[];
for cldx=1:Nspots
    L = T == cldx;
    L2 = find(L);
    % Pick original positions
    for m=1:length(L2)
        P{cldx}(m,:) = Pos(L2(m),:);
    end
    if length(L2) >= 2
        CPosOri{Ncluser} = P{cldx};
        CPos(Ncluser,1) = mean(P{cldx}(:,1));
        CPos(Ncluser,2) = mean(P{cldx}(:,2));
        % Put # of molecule in the cluster
        CPos(Ncluser,3) = length(P{cldx}(:,1));
        Ncluser = Ncluser+1;
    else
        MPos(Nmono,1) = P{cldx}(:,1);
        MPos(Nmono,2) = P{cldx}(:,2);
        MPos(Nmono,3) = 1;
        Nmono = Nmono+1;
    end
end
% Re-calc #
Ncluser = Ncluser-1;
Nmono = Nmono-1;
% Sumarize spots
SumPos = vertcat(CPos, MPos);
% Stat
Nt_ori = length(Pos(:,1));
Nt_mod = length(SumPos(:,1));

% Draw scatter plots
if FigDraw == 1
    % Generate color table
    CloTableOri = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE'};
    %{'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE'};
    CloTable = CloTableOri;
    for k=2:ceil(Ncluser/6)
        CloTable = horzcat(CloTable, CloTableOri);
    end
    % Cluster View
    figure(2); hold on;
    set(gcf,'Position',[0,0,500,500]);
    % Draw cluster
    for k=1:Ncluser
        scatter(CPosOri{k}(:,1),CPosOri{k}(:,2),'.', 'MarkerEdgeColor', CloTable{k});
    end
    % draw monomer
    if ~isempty(MPos)
        scatter(MPos(:,1), MPos(:,2),5,'k+');
    end
    pbaspect([1 1 1]);
    xlabel('(pix)');
    ylabel('(pix)');
    % Delete extra space
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

    % Sum View    
    figure(3);
    set(gcf,'Position',[0,0,500,500]);
    scatter(SumPos(:,1), SumPos(:,2),'k.');
    pbaspect([1 1 1]);
    xlabel('(pix)');
    ylabel('(pix)');
    % Delete extra space
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    % Original view   
    figure(4);
    set(gcf,'Position',[0,0,500,500]);
    scatter(Pos(:,1), Pos(:,2),'k.');
    pbaspect([1 1 1]);
    xlabel('(pix)');
    ylabel('(pix)');
    % Delete extra space
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    % Save figures
    saveas(figure(1), 'Tree.fig')
    saveas(figure(1), 'Tree.png')
    saveas(figure(2), 'ClusterId.fig')
    saveas(figure(2), 'ClusterId.png')
    saveas(figure(3), 'ModSpots.fig')
    saveas(figure(3), 'ModSpots.png')
    saveas(figure(4), 'OriSpots.fig')
    saveas(figure(4), 'OriSpots.png')   
end


end
