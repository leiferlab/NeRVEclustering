function plotTrackData(im,trackData,t)
%%%%%% plotTrackData takes a single frame at time t which is part of an
%%%%%% image sequence used for tracking and scatters the detected peaks and
%%%%%% labels on them

imagesc(im);
hold on
colormap gray
tracks=trackData(trackData(:,end-1)==t,:);
scatter(tracks(:,1),tracks(:,2),'rx');
text(tracks(:,1),tracks(:,2),cellstr(num2str(tracks(:,end))),'VerticalAlignment'...
    ,'bottom', 'HorizontalAlignment','right','color','g');
text(1,1,['t=' num2str(t)],'VerticalAlignment'...
    ,'bottom', 'HorizontalAlignment','right','color','black','fontsize', 40);


axis equal
hold off


