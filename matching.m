%读入图像
close all;
clear all;
% run('Z:\matlab2016a\toolbox\vl_setup.m');
filee='';
file1='2';
file2='3';
Iaa=im2double(imread([file1 '.jpg']));
Ia=Iaa(:,:,1);
Ibb=im2double(imread([file2 '.jpg']));
Ib=Ibb(:,:,1); 
n=size(Ia,1);
m=size(Ia,2);
%调整两幅图像大小一致
k=2;
Ia=imresize(Iaa,[n/k m/k]);
Ib=imresize(Ibb,[n/k m/k]);
% 图像预处理
Ia1 = single(rgb2gray(Ia));
Ib1 = single(rgb2gray(Ib));
% 
% %==========================================================================
%提取SIFT特征点
[fa, da] = vl_sift(Ia1) ;
[fb, db] = vl_sift(Ib1) ;%Each column
thres=2;
%   of F is a feature frame and has the format [X;Y;S;TH], where X,Y
%   is the (fractional) center of the frame, S is the scale and TH is
%   the orientation (in radians).
%   [F,D] = VL_SIFT(I) computes the SIFT descriptors [1] as well. Each
%   column of D is the descriptor of the corresponding frame in F. A
%   descriptor is a 128-dimensional vector of class UINT8.
% %==========================================================================
% 挑选有用的SIFT特征点
%利用计算得到的SIFT特征da，db进行配准
[matches, scores] = vl_ubcmatch(da, db,thres) ; %matches是2xM的矩阵，matches（1，M）与matches（2，M)是对应点的序号
Mpoints_SIFT=[fa(1,matches(1,:))' fa(2,matches(1,:))' zeros(size(matches,2),1) da(:,matches(1,:))'];%每个关键点有三个信息：位置、尺度、方向
Fpoints_SIFT=[fb(1,matches(2,:))' fb(2,matches(2,:))' zeros(size(matches,2),1) db(:,matches(2,:))'];
Mpoints=[fa(1,matches(1,:))' fa(2,matches(1,:))' zeros(size(matches,2),1) ];
Fpoints=[fb(1,matches(2,:))' fb(2,matches(2,:))' zeros(size(matches,2),1) ];
k1=round(size(matches,2)*1.0);
k2=round(size(matches,2)*1.0);
Mpoints=Mpoints(1:k1,:); Fpoints=Fpoints(1:k2,:);
M_SIFT=double(Mpoints_SIFT(:,4:131));
F_SIFT=double(Fpoints_SIFT(:,4:131));
[row column]=size(Mpoints);

   
        %%%%%%%%GL-CATE
        opt.sao=double(M_SIFT);
        opt.sbo=double(F_SIFT);
        opt.outliers = 0.5;
        opt.viz = 0;
        opt.t = 0.9999;
        opt.sparse = 0;
        opt.nsc = 5;
        opt.normalize = 1;
        opt.fa=fa;
        opt.fb=fb;
        opt.da=da;
        opt.db=db;
        opt.thres=thres;
        opt.beta = 2;
        opt.lambda = 1;
        opt.eta=1;
        opt.tol = 1e-6;
        opt.K=5;
        [Transform, C, corr]=glcate_registerO(Fpoints(:,1:2),Mpoints(:,1:2), opt);
        vx=Transform.Y; 
        
% %         % % % % % % % % %Ma-CPD
% %         %Init full set of options %%%%%%%%%%
% %         mo='R';
% %         opt1.outliers = 0.5;
% %         opt1.viz =1;
% %         opt1.t = 0.9;
% %         opt1.sparse = 0;
% %         opt1.nsc = 5;
% %         opt1.corresp=1;
% %         opt1.SiftM=SiftM;
% %         opt1.normalize = 1;
% %         opt1.beta = 2;
% %         opt1.lambda = 3;
% %         opt1.tol = 1e-10;
% %         [Transform, C, corr] = prgls_register(Fpoints(:,1:2), Mpoints(:,1:2),opt1, mo);      
% %         vx=Transform.Y; 
%         Mpoints_warped=[Mpoints_warped zeros(size(Mpoints_warped,1),1)];
%         Fp_resize=[Fp_resize zeros(size(Fp_resize,1),1)];
        % %% %%%%%
        Pos1=[Mpoints(:,2) Mpoints(:,1) Mpoints(:,3)]; 
        % Pos2=[Fpoints(MF(:,1),2) Fpoints(MF(:,1),1) Fpoints(MF(:,1),3)];
        Pos2=[vx(:,2) vx(:,1) zeros(size((vx),1),1)]; 
        % Pos2=[Fpoints(:,2) Fpoints(:,1) Fpoints(:,3)]; 
        Pos1(:,3)=zeros(size(Pos1,1),1);
        Pos2(:,3)=zeros(size(Pos2,1),1);
        % Pos3(:,3)=zeros(size(Pos3,1),1);
        I1_warped=tps_warp(Ia,Pos1,Pos2,'bicubic');
        figure(10); 
%         set(0,'CurrentFigure',1) 
        set (gcf, 'color', [1 1 1]);
        hold off;
        imshow(I1_warped);title('The Transformed Image');
        hold on;
        plot(Pos2(:,2),Pos2(:,1),'.','Color',[1 1 0],'MarkerSize',10) %打点
        hold on;    

        sizea1=size(Iaa,1);    
        sizea2=size(Iaa,2);
        sizeb1=size(Ibb,1);
        sizeb2=size(Ibb,2);
        transVec=size(Ia1,2);
        figure(1); imshow(Ia);hold on;
        plot(Mpoints(:,1),Mpoints(:,2),'.','LineWidth',3,'color','r');hold on
        figure(2); imshow(Ib);hold on;
        plot(Fpoints(:,1),Fpoints(:,2),'.','LineWidth',3,'color','r');hold on
        I=[Ia Ib];
        figure(3); imshow(I);hold on;
        Fpoints1=Fpoints;
        Fpoints1(:,1)=Fpoints1(:,1)+transVec;
        plot(Mpoints(:,1),Mpoints(:,2),'.','Color','r','MarkerSize',10, 'LineWidth', 2);hold on
        plot(Fpoints1(:,1),Fpoints1(:,2),'.','Color','r','MarkerSize',10, 'LineWidth', 2);hold on
        
sizeM=size(Mpoints,1);
btn=-1; 
ith=1;
crt=zeros(sizeM,1);
yellow1=[];
yellow2=[];
while  ith<=sizeM
    senP=Mpoints(ith,:);
    refP=Fpoints(C(ith),:); 
    refPTran(:,1)=refP(:,1)+transVec;
    refPTran(:,2)=refP(:,2);   
    drawn=plot([senP(:,1) refPTran(:,1)],[senP(:,2) refPTran(:,2)],'Color','y','LineWidth',1);hold on
    yellow1=plot(senP(:,1),senP(:,2),'.','Color','y','MarkerSize',5);hold on
    yellow2=plot(refPTran(:,1),refPTran(:,2),'.','Color','y','MarkerSize',5);hold on
    disp(['Currently    ' num2str(ith) '   to     ' num2str(C(ith)) ',  There is ' num2str(sizeM-ith) ' Left' ]);
    [x1,y1,btn]=ginput(1); 
    if btn==2
       break; 
    elseif btn==3
        crt(ith,1)=0;
    elseif btn==1
        crt(ith,1)=ith;
    elseif btn==32
        ith=ith-2;
        delete(yellow1);
        delete(yellow2);
        if ith<=0
            ith=0;
        end
    end
% %     if  (ith==C(ith) ||  ith==C(ith)-1 ||ith==C(ith)+1 ) 
%     
% %     if  (ith==C(ith) || ith==C(ith)-2 || ith==C(ith)-1 ||ith==C(ith)+1 || ith==C(ith) +2) || isnumeric(find(corr==1)==1)==1
%     if  (ith==C(ith) || ith==C(ith)-1 ||ith==C(ith)+1 || isnumeric(find(corr==1)==1)==1)   
%         crt(ith,1)=ith;
%     end
    pause(0.5);
    ith=ith+1;  
    delete(drawn);
end
%===========================================================================
crt1=find(crt~=0);
results=crt1;
save([filee '_' file1 '-' file2 'Thres_' num2str(thres) 'k' num2str(k) '.mat'],'results');
% crt1=find(crt~=0);
% [precise, recall, corrRate] = evaluate(corr, crt1, k1);