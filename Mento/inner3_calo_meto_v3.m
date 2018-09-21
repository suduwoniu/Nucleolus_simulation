%This code is for exploring the optimal distribution of DFC processing factors

%One sub-nucleolus is set in a three dimentional space, whose volumn is
%n1*n2*n3. The pre-rRNA processing factory is set in the center of 3D space, which is a
%sphere with a diameter as the same measured in the real world.  And N
%processing factors are distributed in the 3D space in a sphere whose center is the center of 3D space
%but not in the sphere of pre-rRNA processing factory

%prameters set
n1=50;n2=50;n3=50;nprocess=1837;nrna=175;rrna=6.5;inteval=8;iter=40;pnum1=7;pnum2=25;ratio=11;
%n1,n2,n3:the length, width, height of the three dimential space (1 represent 10 nm)
%nprocess: the total number of the processing factors in the DFC region
%nrna: the total number of pol II which produce pre-rRNA
%pnum1: the radius of the sphere of pre-rRNA processing factory
%pnum2: the radius of the sphere of the total sub-nucleous
%ratio: the raito of nprocess:nrna 

% filedNum is the sites which processing factor could stand in.
filed=zeros(n1,n2,n3);
filedNum=zeros(prod(size(filed)),4);
count=1;
n11=(n1+1)/2;n22=(n2+1)/2;n33=(n3+1)/2;
for xx=1:n1
    for yy=1:n2
        for zz=1:n3
            %processing factors are not in the sphere of pre-rRNA
            %processing factory
            if ((xx-n11)^2+(yy-n21)^2+(zz-n31)^2)>=(pnum1^2)
                %if ((xx-n11)^2+(yy-n21)^2+(zz-n31)^2)<=(pnum2^2)
                filedNum(count,1)=count;
                filedNum(count,2)=xx;
                filedNum(count,3)=yy;
                filedNum(count,4)=zz;
                filedNum(count,5)=1;
                count=count+1;
               %end
            end
        end
    end
end

%random choose processing factors position in the filedNum
count=count-1;
p=randperm(count);
p1=p(1:nprocess);
randomFactor=zeros(nprocess,7);
count2=0;
for xx=1:count
    if (ismember(filedNum(xx,1),p1)==1)
        count2=count2+1;
        randomFactor(count2,1)=count2;
        randomFactor(count2,2)=filedNum(xx,2);
        randomFactor(count2,3)=filedNum(xx,3);
        randomFactor(count2,4)=filedNum(xx,4);
    end
end

%the possible positions for pre-rRNA producing place 
count0=0;
for xx=18:33
    for yy=18:33
        for zz=18:33
            if sqrt((xx-n11)^2+(yy-n22)^2+(zz-n33)^2)>=7 && sqrt((xx-n11)^2+(yy-n22)^2+(zz-n33)^2)<=8
                count0=count0+1;
                inner(count0,1:3)=[xx,yy,zz];
            end
        end
    end
end

%the pre-rRNAs are producing in a three-bunch sites
count1=0;
pol=zeros(1,5);
for k=1:count0
    [theta,phi,r]=cart2sph(inner(k,1)-25.5,inner(k,2)-25.5,inner(k,3)-25.5);
    if phi>-1/4*pi && phi<=1/4*pi 
        if (theta>=-1/8*pi && theta <=1/8*pi) || (theta>=13/24*pi && theta<=19/24*pi) || (theta>=-19/24*pi && theta<=-13/24*pi)
            count1=count1+1;
            pol(count1,1:3)=inner(k,1:3);
        end
    end
end

%randomly choose the pre-rRNA producing site
pinner=randperm(count1);
pinner1=pinner(1:nrna);
RRNA=pol(pinner1,:);
inner=RRNA;

%the pre-rRNA choosing which processing factor to bind, the pre-rRNA bind
%to the processing factor whose distance between them are shortest and the
%processing factor is no other pre-rRNA binding
nrna;
for n=1:nprocess
    dist0=99999999999;idx=0;
    for k=1:nrna
        if inner(k,4)<ratio
            dist1=sqrt((randomFactor(n,2)-inner(k,1))^2+(randomFactor(n,3)-inner(k,2))^2+(randomFactor(n,4)-inner(k,3))^2);
            if dist1<dist0
                dist0=dist1;
                idx=k;
            end
        end
    end
    randomFactor(n,5)=idx;
    randomFactor(n,6)=dist0;
    inner(idx,4)=inner(idx,4)+1;
end

%Calculation of each processing factor's binding ability to pre-rRNA
%randomFactor is the position of processing factors and the 7th column is
%the binding ability
%inner is the 3D position of pre-rRNA producing site
%n1,n2,n3 is the length, width and height of 3D space
%3 is range for calculation the density of one processing factors
[randomFactor_res]=CalEffect_all_new_v4(randomFactor,inner,n1,n2,n3,3);
%moveSteps(tem_track,inner,step,select_num,n1,n2,n3,ratio);
%[randomFactor_temp,inner_1]=moveSteps(randomFactor_res,inner,10,100,n1,n2,n3,ratio);
randomFactor_temp=randomFactor_res;
inner_temp=inner;

%factor_record: record the movement of each disturbtion
factor_record=zeros(nprocess,4,10000);
cn=1;
factor_record(:,1,cn)=randomFactor_res(:,2);factor_record(:,2,cn)=randomFactor_res(:,3);
factor_record(:,3,cn)=randomFactor_res(:,4);factor_record(:,4,cn)=randomFactor_res(:,7);
%give 10000 disturbtions for this system
for i=1:10000
    %function [rand_sort_res,inner]=moveSteps_max(tem_track,inner,step,select_num,max_point,n1,n2,n3,ratio,rrna);
    %choose the 500 processing factors to move to another place
    [randomFactor_temp_1,inner_temp_1]=moveSteps_max_v1(randomFactor_temp,inner_temp,4,500,500,n1,n2,n3,ratio,rrna);
   step=5;select_num=200;
   %[randomFactor_temp_1,inner]=moveSteps(randomFactor_temp_1,inner,step,select_num,n1,n2,n3,ratio,rrna);
    %[randomFactor_temp_1,inner_temp_1]=moveSteps(randomFactor_temp,inner_temp,10,200,n1,n2,n3,ratio,rrna);
    if sum(randomFactor_temp_1(:,7))>sum(randomFactor_temp(:,7))
        randomFactor_temp=randomFactor_temp_1;
        inner_temp=inner_temp_1;
        cn=cn+1;
        factor_record(:,1,cn)=randomFactor_temp(:,2);
        factor_record(:,2,cn)=randomFactor_temp(:,3);
        factor_record(:,3,cn)=randomFactor_temp(:,4);
        factor_record(:,4,cn)=randomFactor_temp(:,7);
    else
        rand_1=rand;
        if rand_1<0.3
         randomFactor_temp=randomFactor_temp_1;
        inner_temp=inner_temp_1;
        cn=cn+1;
        factor_record(:,1,cn)=randomFactor_temp(:,2);
        factor_record(:,2,cn)=randomFactor_temp(:,3);
        factor_record(:,3,cn)=randomFactor_temp(:,4);
        factor_record(:,4,cn)=randomFactor_temp(:,7);
        end
    end
end

%{
plot3(factor_record(:,1,cn),factor_record(:,2,cn),factor_record(:,3,cn),'.')
axis([0,50,0,50,0,50]);
hold on;
plot3(inner(:,1),inner(:,2),inner(:,3),'.');
hold off;
plot3(factor_record(:,1,cn),factor_record(:,2,cn),factor_record(:,3,cn),'.');
plot3(factor_record(:,1,cn),factor_record(:,2,cn),factor_record(:,3,cn),'.');
%}
save('random_v4_s3_4_500_500_res.mat','randomFactor_res');
save('radom_v4__s3_4_500_500.mat','factor_record');

%plot3(randomFactor_temp(:,2),randomFactor_temp(:,3),randomFactor(:,4),'.');
%{
type_w=fopen('model_res/effect_all_inball_outerSqure.txt','wt');
for n=1:nprocess
    fprintf(type_w,'%d\t%d\t%d\t%d\t%d\t%d\n',randomFactor(n,2),randomFactor(n,3),randomFactor(n,4),randomFactor(n,5),randomFactor(n,6),randomFactor_res(n,7)*(10^6));
end
%}

