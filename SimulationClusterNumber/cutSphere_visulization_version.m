%R: the distance of cluster center to the space center (the space center is [0,0,0])
%r: the radius of the each cluster
%number_shape: the cluster numbers in one sub-nucleus(space center)
number_shape=22;
R=249;r=68;

%load the possible sites for the center of each cluster
%load matrix "Rcenter"
load('../../Rcenter_R249_r68.mat');

%{
Rcenter=zeros(1,3);
count0=0;
for x=Rrag
    for y=Rrag
        for z=Rrag
            if sqrt(x^2+y^2+z^2)<=R+0.5 && sqrt(x^2+y^2+z^2)>=R-0.5
                count0=count0+1;
                Rcenter(count0,1)=x;
                Rcenter(count0,2)=y;
                Rcenter(count0,3)=z;
            end
        end
    end
end
%}

%randomly choose the position of cluster center from the matrix "Rcenter"
Rselect=zeros(number_shape,3);
Rcenter_1=Rcenter;
for num=1:number_shape
    Rs=size(Rcenter_1);
    del=zeros(1);
    p=randperm(Rs(1));
    count1=0;
    Rselect(num,:)=Rcenter_1(p(1),:);
    for n=1:Rs(1)
        %be sure that the distance of the two cluster centers are over 1.5
        %fold of cluster radius
        if sqrt((Rcenter_1(n,1)-Rselect(num,1))^2+(Rcenter_1(n,2)-Rselect(num,2))^2+(Rcenter_1(n,3)-Rselect(num,3))^2)<=1.5*r
            count1=count1+1;
            del(count1)=n;
        end
    end
    Rcenter_1(del,:)=[];
end
%save('Rcenter.mat','Rcenter');
%save_Rcenter=fopen('Rcenter_save.txt','wt');
%{
Rs=size(Rcenter);
for n=1:Rs(1)
    fprintf(save_Rcenter,'%d\t%d\t%d\n',Rcenter(n,1),Rcenter(n,2),Rcenter(n,3));
end
%}                      

%shape=strcat('shape_',num2str(shape_num));
%mkdir(shape);
%dir=cd;
%dir_s=[dir,'/',shape];
%cd(dir_s);
%save_select=['Rselect-',shape,'.mat'];
%save(save_select,'Rselect');
%save_Rselect=fopen('Rselect.txt','wt');
%for n=1:number_shape
%    fprintf(save_Rselect,'%d\t%d\t%d\n',Rselect(n,1),Rselect(n,2),Rselect(n,3));
%end


s=size(Rselect);
radius=zeros(s(1),1);
radius(:)=r;
vert=Rselect;
a=vert(:,1);
b=vert(:,2);
c=vert(:,3);
%{
[x,y,z]=sphere(30);
for k=1:s(1)
    surf(radius(k)*x+a(k),radius(k)*y+b(k),radius(k)*z+c(k));
    %shading interp;
    hold on;
end
hold off;
%}

%randomly rotate the cluster system 200 times
number_ind=zeros(200,1);
for num=1:200
%randomly choose the rotate angle 
azel=[360*rand,360*rand];
alpha=360*rand;
newxyz=rotateEllipsoid(vert,azel,alpha);
nd=zeros(s(1),1);
nd(:)=rand*10;
newxyz(:,3)=newxyz(:,3)+nd;

%newxyz=vert;
%record: record all points of clusters around z=0 (get the cross section of
%z=0)
record=zeros(1,4);
count=0;
rag=linspace(-800,800,100);
zrag=linspace(-3,3,30);
for x=rag
    for y=rag
        for z=zrag
            for n=1:s(1)
                if (x-newxyz(n,1))^2+(y-newxyz(n,2))^2+(z-newxyz(n,3))^2<=radius(n)^2
                    count=count+1;
                    record(count,1)=x;record(count,2)=y;record(count,3)=z;
                    record(count,4)=n;
                end
            end
        end
    end
end 

ind=unique(record(:,4));
if length(ind)>1
    count5=hist(record(:,4),ind);
    ind_n=zeros(1);
    count6=0;
    for nn=1:length(ind)
        if count5(nn)>0;
            count6=count6+1;
            ind_n(count6)=ind(nn);
        end
    end

    
    ind_num=length(ind_n);
    ind_center=zeros(ind_num,3);
    for nn=1:ind_num
        ind_center(nn,1)=sum(record(record(:,4)==ind_n(nn),1))/length(record(record(:,4)==ind_n(nn),1));
        ind_center(nn,2)=sum(record(record(:,4)==ind_n(nn),2))/length(record(record(:,4)==ind_n(nn),2));
        ind_center(nn,3)=sum(record(record(:,4)==ind_n(nn),3))/length(record(record(:,4)==ind_n(nn),3));
    end
%
    ind_radius=zeros(ind_num,1);
    for nn=1:ind_num
        obj=record(record(:,4)==ind_n(nn),:);
        os=size(obj);
        max_o=0;
        for aa=1:os(1)
            rd=sqrt((obj(aa,1)-ind_center(nn,1))^2+(obj(aa,2)-ind_center(nn,2))^2+(obj(aa,3)-ind_center(nn,3))^2);
            if rd>max_o
                max_o=rd;
            end
        end
        ind_radius(nn)=max_o;
    end

    ind_label=ind_n;
    for kk=1:ind_num
        for zz=1:ind_num
            if ind_n(kk)~=ind_n(zz)
                d=sqrt((ind_center(kk,1)-ind_center(zz,1))^2+(ind_center(kk,2)-ind_center(zz,2))^2+(ind_center(kk,3)-ind_center(zz,3))^2);
                if d<max(ind_radius(kk),ind_radius(zz))
                    [md,idx]=max([ind_radius(kk),ind_radius(zz)]);
                    ind_radius(kk)=md;
                    ind_radius(zz)=md;
                    if idx==1
                        rx=ind_center(kk,1);
                        ry=ind_center(kk,2);
                        rz=ind_center(kk,3);
                    %ind_center(kk,:)=[rx,ry,rz];
                    %ind_center(zz,:)=[rx,ry,rz];
                    else
                        rx=ind_center(zz,1);
                        ry=ind_center(zz,2);
                        rz=ind_center(zz,3);
                    end
                    ind_center(ind_n==ind_n(kk),1)=rx;
                    ind_center(ind_n==ind_n(zz),1)=rx;
                    ind_center(ind_n==ind_n(kk),2)=ry;
                    ind_center(ind_n==ind_n(zz),2)=ry;
                    ind_center(ind_n==ind_n(kk),3)=rz;
                    ind_center(ind_n==ind_n(zz),3)=rz;
                    ind_radius(ind_n==ind_n(kk))=rd;
                    ind_radius(ind_n==ind_n(zz))=rd;
                    ind_n(zz)=ind_n(kk);
                end
            end
        end
    end
    ind_uniq=unique(ind_n);
    number_ind(num)=length(ind_uniq);
else
    number_ind(num)=0;
end


plot3(record(:,1),record(:,2),record(:,3),'.','color',[255/255,112/255,115/255]);
view(90,90);
axis([-900,900,-900,900,-900,900]);
hold on;
plot3(ind_center(:,1),ind_center(:,2),ind_center(:,3),'*','color',[255/255,112/255,115/255]);
title(num2str(length(ind_uniq)));
hold off;
saveas(gcf,strcat(num2str(num),'.tif'));

end
%cd(dir);


