% plotnan -- packr([x y]) first

function h=plotnan(x,y,a,b,c,d,e);
data=packr([x y]);
xx=data(:,1);
yy=data(:,2);
plot(xx,yy,a,b,c,d,e);
