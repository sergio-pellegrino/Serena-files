clear
clc
close all

iter = 1;

while exist('designVariables.txt') ~= 2
    pause(10)
end

while exist('designVariables.txt') == 2
    
    iter = iter + 1;
    
    clearvars -except iter
    close all
    
    cutPlane_data = importdata('designVariables.txt');
    %II = find(cutPlane_data{1}=='[') + 1;
    %JJ = find(cutPlane_data{1}==']') - 1;
    %cutPlane = str2num(cutPlane_data{1}(II:JJ));

    % cutPlane = 0.1;

    xxp = linspace(-1,1,100);
    yyp = linspace(-1,1,100);
    [xp,yp]=meshgrid(xxp,yyp);

    x = (19.55 + 37.95.*xp + 10.35.*yp + 10.35.*xp.*yp).*(xp<-1/3)...
        + (16.1 + 6.9.*yp).*(xp>1/3)...
        + (13.4054 + 13.8.*xp + 6.9.*yp -17.1484.*xp.^2).*(xp>=(-1/3)).*(xp<=1/3);

    y = (16.1 + 6.9.*yp).*(xp<-1/3)...
        + (19.55 - 37.95.*xp + 10.35.*yp - 10.35.*xp.*yp).*(xp>1/3)...
        + (13.4054 - 13.8.*xp + 6.9.*yp - 17.1484.*xp.^2).*(xp>=(-1/3)).*(xp<=1/3);

    Nh = 2;
    Nl = 0;
    z = zeros(length(xxp),length(yyp));
    for h = 0:Nh
        for l = 0:Nl
            z = z + ((cos((2*h+1)*(pi/2).*xp)+sin((h)*(pi).*xp))...
            .*(cos((2*l+1)*(pi/2).*yp)+sin((l)*(pi).*yp))).^2;
        end
    end
    % 
    % figure
    % surf(x,y,z)
    % figure
    [con_mat, obj]=contour(x,y,z,[cutPlane_data cutPlane_data]);
    % axis equal
    % 
    % z_indices = find(abs(z-0.8)<0.001);
    % xim = x(z_indices);
    % yim = y(z_indices);
    % 
    % figure
    % plot(xim,yim,'.')

    number = 1;
    i = 1;
    while i<=length(con_mat)
        lev = con_mat(1,i);
        cnt = con_mat(2,i);
        xtemp{number}=con_mat(1,i+(1:cnt));
        ytemp{number}=con_mat(2,i+(1:cnt));

        i = i+cnt+1;
        number=number+1;
    end;

    % figure
    % plot(xtemp{1},ytemp{1},'-b')
    % axis equal
    % hold on
    % plot(xtemp{2},ytemp{2},'-r')
    % plot(xtemp{3},ytemp{3},'-g')
    % plot(xtemp{4},ytemp{4},'-y')
    % plot(xtemp{5},ytemp{5},'-k')
    % hold off

    % save text file
    for jj = 1:length(xtemp)
        fid = fopen(join(['cutout-',num2str(jj),'-iter',num2str(iter)]),'w+');
        fprintf(fid, '%f\n',[xtemp{jj},ytemp{jj}]);
        fclose(fid);
    end

    % system('abaqus cae script=H:\Serena\Opt_Mapping1\geometry-2.py')
   
    delete('designVariables.txt')
    
    while exist('designVariables.txt') ~= 2
        pause(10)
    end
    
end