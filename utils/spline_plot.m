function spline_plot(name,points,corners,coeff,method,color)
cornerpoints = points(:,corners); % coodinate
propsP = {'LineStyle','none','Marker','o','MarkerEdge','black','MarkerFaceColor','red','MarkerSize',6};
f = figure;
ax = axes(f);
axis(ax,'equal')
grid(ax,'off');
set(ax,'XTick',[]);
set(ax,'YTick',[]);
box(ax,'off');
axis(ax,'off');
daspect(ax,[1,1,1]);
corner_handle=...
line(ax,'XData',cornerpoints(1,:),'YData',cornerpoints(2,:),propsP{:});

if method(1) == 'p'
polygon_handle=line(ax,'XData',points(1,:),'YData',points(2,:),'PickableParts','none');
return
end
title(ax,name);
sh = spline_helper;
sh.set_points(corners);
sh.set_coeff(coeff);
sh.try_draw_coeff(ax,corners,coeff,0,method,color);

end