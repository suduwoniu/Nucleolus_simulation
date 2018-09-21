load('random_v4_s3_4_500_500_res.mat');
s=size(factor_record);
aviobj=VideoWriter('change_color.avi');
aviobj.FrameRate=10;
open(aviobj);
qn=0;
for n=1:s(3)
    if sum(factor_record(:,1,n))>0
        tem1=factor_record(:,:,n);
        figure(1);
        %set(gcf,'color','none');
        plot3(tem1(:,1),tem1(:,2),tem1(:,3),'.','color',[79/255,197/255,210/255],'MarkerSize',20,'MarkerFaceColor','b');
        hold on;
        plot3(inner(:,1),inner(:,2),inner(:,3),'.','color','black','MarkerSize',8);
        set(gca,'color','white');
        hold on;
        set(gcf,'color','white');
        hold on;
        box on;
        set(gca,'LineWidth',3);
        set(gca,'Xcolor','black');
        set(gca,'Ycolor','black');
        set(gca,'Zcolor','black');
        %view(x,y);
        hold off;
        axis([0,50,0,50,0,50]);
        frame=getframe(gcf);
        im=frame2im(frame);
        writeVideo(aviobj,im);
        qn=qn+1;

    end
end
close(aviobj);

aviobj=VideoWriter('randomFactor_state_5_200_30_final.avi');
aviobj.FrameRate=10;
open(aviobj);
final=factor_record(:,:,qn);
for i=-37:2:323
    plot3(final(:,1),final(:,2),final(:,3),'.','color',[79/255,197/255,210/255],'MarkerSize',15);
    hold on;
    plot3(inner(:,1),inner(:,2),inner(:,3),'.','color','black','MarkerSize',8);
    view(i,30);
    set(gca,'color','white');
        hold on;
        set(gcf,'color','white');
        hold on;
        box on;
        set(gca,'LineWidth',3);
        set(gca,'Xcolor','black');
        set(gca,'Ycolor','black');
        set(gca,'Zcolor','black');
        %view(x,y);
        hold off;
        axis([0,50,0,50,0,50]);
        frame=getframe(gcf);
        im=frame2im(frame);
        writeVideo(aviobj,im);
    %pause(0.06);
end
close(aviobj);

%{
aviobj=VideoWriter('track_25.avi');
aviobj.FrameRate=10;
sa=size(time_track);
open(aviobj);
for n=1:53
    tem1=time_track(1:p_num,:,n);
    tem2=time_track(p_num+1:p_num+f_num,:,n);
    tem3=time_track(p_num+f_num+1:total_num,:,n);
    figure(1);
    scatter3(tem1(:,2),tem1(:,3),tem1(:,4),60,'b','.');
    hold on;
    scatter3(tem2(:,2),tem2(:,3),tem2(:,4),60,'m','.');
    grid off;
    hold off;
    frame=getframe(gcf);
    im=frame2im(frame);
    writeVideo(aviobj,im);
end
close(aviobj);
%}