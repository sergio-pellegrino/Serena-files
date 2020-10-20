clear all
close all
clc

set(0,'defaultaxeslinewidth',2)
set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesfontweight','bold')
set(0,'defaulttextfontweight','bold')
set(0,'defaultaxesfontsize',14)
set(0,'defaulttextfontsize',14)

xxp = linspace(-1,1,100);
yyp = linspace(-1,1,100);
[xp,yp]=meshgrid(xxp,yyp);

x = (19.55 + 37.95.*xp + 10.35.*yp + 10.35.*xp.*yp).*(xp<-1/3)...
    + (16.1 + 6.9.*yp).*(xp>1/3)...
    + (13.4054 + 13.8.*xp + 6.9.*yp -17.1484.*xp.^2).*(xp>=(-1/3)).*(xp<=1/3);

y = (16.1 + 6.9.*yp).*(xp<-1/3)...
    + (19.55 - 37.95.*xp + 10.35.*yp - 10.35.*xp.*yp).*(xp>1/3)...
    + (13.4054 - 13.8.*xp + 6.9.*yp - 17.1484.*xp.^2).*(xp>=(-1/3)).*(xp<=1/3);

Nh = 0;
Nl = 2;
z = zeros(length(xxp),length(yyp));
for h = 0:Nh
    for l = 0:Nl
%         z = z + ((cos((2*h+1)*(pi/2).*xp)+sin((h)*(pi).*xp))...
%             .*(cos((2*l+1)*(pi/2).*yp)+sin((l)*(pi).*yp))).^2;

        z = z + ((cos((2*h+1)*(pi/2).*xp))...
            .*(cos((2*l+1)*(pi/2).*yp))).^2;
    end
end
% 
figure
surf(x,y,z)
xlabel('x')
ylabel('y')
zlabel('z')
%xlim([-20 25])
%ylim([-20 25])

figure
surf(xp,yp,z)
xlabel('x_p')
ylabel('y_p')
zlabel('z')

figure
z_zeros = zeros(length(x),length(y));
surf(x,y,z_zeros)
xlabel('x')
ylabel('y')
axis equal

figure
surf(xp,yp,z_zeros)
xlabel('x_p')
ylabel('y_p')
axis equal

figure
contour(x,y,z)
xlabel('x')
ylabel('y')
axis equal

figure
contour(xp,yp,z)
xlabel('x_p')
ylabel('y_p')
axis equal

%% cutout

figure
[con_mat, obj]=contour(x,y,z,[1.35 1.35]);
axis equal

z_indices = find(abs(z-0.8)<0.001);
xim = x(z_indices);
yim = y(z_indices);

figure
plot(xim,yim,'.')

figure
number = 1;
i = 1;
while i<=length(con_mat)
    lev = con_mat(1,i);
    cnt = con_mat(2,i);
    xtemp{number}=con_mat(1,i+(1:cnt));
    ytemp{number}=con_mat(2,i+(1:cnt));
    
    plot(xtemp{number},ytemp{number})
    axis([-15 25 -15 25])
    xlabel('x')
    ylabel('y')
    axis square
    hold on
    i = i+cnt+1;
    number=number+1;
end;
hold off

% % save text file
% for jj = 1:length(xtemp)
%     fid = fopen(join(['cutout-',num2str(jj),'-iter1']),'w+');
%     fprintf(fid, '%f\n',[xtemp{jj},ytemp{jj}]);
%     fclose(fid);
% end