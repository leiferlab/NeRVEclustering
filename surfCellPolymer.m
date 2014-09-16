subplot(1,2,1);

    if size(cellF.Ff,1)==size(cellcoord.x,1)-1
        cellF.Ff=[cellF.Ff(end,:);cellF.Ff];
    end
    if size(cellF.Ffr,1)==size(cellcoord.x,1)-1
        cellF.Ffr=[cellF.Ffr(end,:);cellF.Ffr];
    end
    
    
surf(cellcoord.x,cellcoord.y,cellcoord.z,cellF.Ff,'facecolor','interp','facealpha',1,'edgecolor','none');
axis equal
hold on

 if ~isempty(polyFit);
            for pi=1:length(polyFit);
                plot3(polyFit(pi).xs,polyFit(pi).ys,...
                    polyFit(pi).zs,'black','LineWidth',3);
            end

 end

 subplot(1,2,2);
 
 h=pcolor([cellF.Ff;cellF.Ff;cellF.Ff]);
set(h,'edgecolor','none','facecolor','interp')
hold on
offset=2;
      for iPolymer=1:length(polyFit);          

        plot(polyFit(iPolymer).rc(:,1),...
            polyFit(iPolymer).rc(:,2)+offset,'black','linewidth',3);

         plot(polyFit(iPolymer).rc(:,1),...
            polyFit(iPolymer).rc(:,2)-size(cellF.Ff,1)+offset,'black','linewidth',3);
       
                 plot(polyFit(iPolymer).rc(:,1),...
            polyFit(iPolymer).rc(:,2)-2*size(cellF.Ff,1)+offset,'black','linewidth',3);
      end
        
        hold off
        set(gca,'ydir','normal')
        set(gca,'ytick',linspace(1,size(cellF.Ff,1),3));
        set(gca,'yticklabel',({'0', 'pi', '2pi'}));
        set(gca,'xticklabel',[]);
        ylim([1,size(cellF.Ff,1)]);
