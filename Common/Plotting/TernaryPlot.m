function TernaryPlot(a,b,c,pointColors)

nP=1000;
% img=ones(nP,nP,3);
% [xMat,yMat]=meshgrid(linspace(0,1,nP),linspace(0,1,nP));
% 
% b=xMat-yMat/sqrt(3);
% c=2*yMat/sqrt(3);
% a=1-(b+c);
% disallowed=(a>1)|(a<0)|(b>1)|(b<0)|(c>1)|(c<0);
% a(disallowed)=1;b(disallowed)=1;c(disallowed)=1;
% a(~disallowed)=1-a(~disallowed);b(~disallowed)=1-b(~disallowed);c(~disallowed)=1-c(~disallowed);
% 
% colorExp=2;
% img(:,:,1)=a;img(:,:,2)=b;img(:,:,3)=c;
% %img=1-img;
% figure;
% image(img.^colorExp);

plot(nP*[0,1],nP*[0,0]+1,'-k','LineWidth',1);
hold on;
plot(nP*[0,1/2],nP*[0,sqrt(3)/2]+1,'-k','LineWidth',1);
plot(nP*[1,1/2],nP*[0,sqrt(3)/2]+1,'-k','LineWidth',1);

plot(nP*[1/4,3/4],nP*sqrt(3)/4*[1,1],'--k','LineWidth',1);
plot(nP*[1/2,3/4],nP*sqrt(3)/4*[0,1],'--k','LineWidth',1);
plot(nP*[1/2,1/4],nP*sqrt(3)/4*[0,1],'--k','LineWidth',1);
axis off equal xy;


x=(2*b+c)/2; x=x+0.005*randn(size(x));y=sqrt(3)*c/2;y=y+0.005*randn(size(y));

scatter(round(nP*x),round(nP*y),100,pointColors,'filled');