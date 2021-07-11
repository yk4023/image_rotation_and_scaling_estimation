%--------------------------------------------------------------------------
%  Examples of usage of the method proposed in:
%
%  [1] C. Chen, J. Ni, and Z. Shen, "Effective estimation of image rotation
% angle using spectral method," Signal Processing Letters, IEEE, vol. 21,
% no. 7, pp. 890--894, July 2014.
%
% The 'Lena' image used in this example can be downloaded from the USC-SIPI
% image database on:
% http://sipi.usc.edu/database/download.php?vol=misc&img=4.2.04
% For copyright information, please go to:
% http://sipi.usc.edu/database/copyright.php
%
%--------------------------------------------------------------------------
% This code is provided only for research purposes.
%--------------------------------------------------------------------------
% Clear all variables and close all figures
clear all; 
close all;
Path = 'resamping/test';
File = dir(fullfile(Path,'*.png'));
name = {File.name}';
t = 0.085;
totle = 0;
table = zeros(length(name),9);
TTT = zeros(256);
% TTT = open('rotate.fig');
for i = 1:length(name)
    fprintf('%2.2f %%\n', i/length(name)*100);
    % 读取图像
    % RGB = imread('resamping/test/3_25_1.5f.png');
     fprintf('true img name: %s\n', name{i});
    img_name = [Path,'/',name{i}];
    img_r = imread(img_name);
    fac = 0;
    ro = 0;
    factor = zeros(1,4);
    angle = 0;
    d=0;
    % 获取绿色通道
    % img = RGB(:,:,2);
    % 角度设置
%     angle = 25;
%     factor = 1.5;
%     fprintf('true factor: %f\n', factor);
    % 插值核，kernel = 'bicubic';
%     fprintf('true angle: %d\n', angle);
    % 旋转图像
    % img_r = imrotate(img,angle);
    % img_r = imresize(img_r,factor,'bicubic');
    % compute Txx，计算Txx
    block_size = 256;
    Txx = CalculateTxx(img_r,block_size);
    Txx = fftshift(Txx);
    [h,w] = size(Txx);
    Txx(h/2-3:h/2+3,w/2-3:w/2+3)=0;
    Txx = fftshift(Txx);
    cha = getcha(Txx,0);
    table(i,3)=cha;
    TTT(Txx>60)=Txx(Txx>60);
%       if mod(i,10)==0 
%          figure('Name',name{i}),imshow(Txx,[]);
%      end
%      figure('Name',name{i}),imshow(Txx,[]);
%     imshow(Txx,[]);
%     hold on;
%     imshow(uint8(Txx))
%     imshow(uint8(fftshift(Txx)))
    % 判断是否经过重采样，经过什么样的重采样
    if cha>t
        [x,y] = getdp(Txx,20,h);
        if sum(x==1)
            y = y(x==1);
            y = y(1);
            d = y(1);
            %缩放系数可能有4种结果
             if y > block_size/2
                    y = block_size-y;
             end
             
            factor(1) = 1/(1-(y-1)/block_size);
            factor(2) = 1/((y-1)/block_size);
            factor(3) = 1/(2-(y-1)/block_size);
            factor(4) = 1/(1+(y-1)/block_size);
            fac = 1;
        elseif sum(y==1)
            x = x(y==1);
            x = x(1);
            d = x(1);
             if x > block_size/2
                    x = block_size-x;
             end
            factor(1) = 1/(1-(x-1)/block_size);
            factor(2) = 1/((x-1)/block_size);
            factor(3) = 1/(2-(x-1)/block_size);
            factor(4) = 1/(1+(x-1)/block_size);
            fac = 1;
        else
            if length(unique(x))~=length(x)
                [B, I] = unique(x, 'first');
                rep_Num = setdiff(1:numel(x), I);
                req_x = x(x==x(rep_Num(1)));
                req_y = y(x==x(rep_Num(1)));
%                 if length(req_x)>2
%                     req_x=req_x(end-1:end);
%                     req_y = sortrows(req_y',1);
%                     req_y=req_y(end-1:end)';
%                 end
                x_p = req_x(1);
                y_p = min(req_y);
                estimatedAngle= atan((x_p-1)/(block_size-y_p+2))*180/pi*2;
                angle = round(estimatedAngle);
            else
                estimatedAngle= EstimateAngle(Txx,0);
                angle = round(estimatedAngle(4));
            end
            
            
%             fprintf('angle estimated using the proposed method: %f\n', estimatedAngle(1));
            ro = 1;
        end
    end
    
    %检测是否经过第二次重采样
    if fac ==1
        Txx(1:end,1)=0;
        Txx(1:end,end)=0;
        Txx(1,1:end)=0;
        Txx(end,1:end)=0;
        factor(factor>2.5)=0;
        factor(factor<0.5)=0;
        cha=getcha(Txx,1);
        table(i,4)=cha;
        if cha>t
            ro = 1;
            [x,y] = getdp(Txx,30,h);
            xy = sortrows([x',y'],1);
            x = xy(1:end,1)';
            y = xy(1:end,2)';
            if length(unique(x))~=length(x)
                [B, I] = unique(x, 'first');
                rep_Num = setdiff(1:numel(x), I);
                req_x = x(x==x(rep_Num(1)));
                req_y = y(x==x(rep_Num(1)));
                if req_x(1)==d && length(rep_Num)~=1
                    req_x = x(x==x(rep_Num(2)));
                    req_y = y(x==x(rep_Num(2)));
                end
                if length(req_x)>2
                    req_x=req_x(end-1:end);
                    req_y = sortrows(req_y',1);
                    req_y=req_y(end-1:end)';
                end
                y_d = abs(req_y(2)-req_y(1));
                req = 2;
                while (y_d>d+2 || y_d<d-2)&&(d~=0)
                    req_x = x(x==x(rep_Num(req)));
                    req_y = y(x==x(rep_Num(req)));
                    if length(req_x)>2
                        req_x=req_x(end-1:end);
                        req_y = sortrows(req_y',1);
                        req_y=req_y(end-1:end)';
                    end
                    y_d = abs(req_y(2)-req_y(1));
                    req = req + 1;
                    if(req>length(rep_Num))
                        break;
                    end
                end
                req =1;
                x_d = req_x(1) + y_d;
                need_x = req_x;
                need_y = req_y;
                x_dd = x_d;
                y_dd = y_d;
                labelx=0;
                while (x_d<256)&&~sum(sum(Txx(x_dd-1:x_dd+1,max(need_y)-1:max(need_y)+1)>60))&&(d~=0)
                    if(req>length(rep_Num))
                        break;
                    end
                    need_x = x(x==x(rep_Num(req)));
                    if need_x(1)==labelx
                        for j = 1:length(need_x)
                        req = req+1;
                        end
                        continue;
                    end
                    need_y = y(x==x(rep_Num(req)));
                    if length(need_x)>2
                        need_x=need_x(end-1:end);
                        need_y = sortrows(need_y',1);
                        need_y=need_y(end-1:end)';
                    end
                    labelx=need_x(1);
                    y_dd = abs(need_y(2)-need_y(1));
                    x_dd = need_x(1)+ y_dd;
                    if x_dd>255
                        break;
                    end
                    req = req + 1;
                    
                    if sum(sum(Txx(x_dd-1:x_dd+1,max(need_y)-1:max(need_y)+1)>60))&&y_dd==d
                        req_x = need_x;
                        req_y = need_y;
                        y_d = y_dd;
                        x_d = need_x(1) + y_d;
                        break;
                    end
                end
                if (x_d<255)&&sum(sum(Txx(x_d-1:x_d+1,max(req_y)-1:max(req_y)+1)>60))
                    if max(req_y)>h/2
                        estimatedAngle= atan((block_size-max(req_y)+1)/(block_size-x_d+1))*180/pi*2;
                    else
                        estimatedAngle= atan((block_size-max(req_y)+2)/(block_size-x_d+2))*180/pi;
                    end
                elseif req_x(1)>y_d
                    if sum(sum(Txx(req_x(1)-y_d-1:req_x(1)-y_d+1,min(req_y)-1:min(req_y)+1)>60))||sum(sum(Txx(req_x(1)-y_d-1:req_x(1)-y_d+1,min(req_y)-1:min(req_y)+1)>60))
                        x_p = req_x(1)-y_d;
                        y_p = min(req_y);
                        estimatedAngle= atan((block_size-y_p+2)/(block_size-x_p+2))*180/pi;
                    end
                else
                   x_p = req_x(1);
                   y_p = min(req_y);
                   if y_p<h/2
                       estimatedAngle= atan((y_p-1)/(x_p-1))*180/pi*2;
                   else
                       estimatedAngle= atan((x_p-1)/(block_size-y_p))*180/pi*2;
                   end
                end
                if estimatedAngle>90
                    estimatedAngle=estimatedAngle-90;
                end
                angle = estimatedAngle;
            end

        end
    elseif ro ==1
        if fac == 0
            Txx(1:end,1)=0;
        Txx(1:end,end)=0;
        Txx(1,1:end)=0;
        Txx(end,1:end)=0;
        factor(factor>2.5)=0;
        factor(factor<0.5)=0;
        cha=getcha(Txx,1);
        table(i,4)=cha;
        if cha>t
            ro = 1;
            [x,y] = getdp(Txx,10,h);
            xy = sortrows([x',y'],1);
            x = xy(1:end,1)';
            y = xy(1:end,2)';
            if length(unique(x))~=length(x)
                [B, I] = unique(x, 'first');
                rep_Num = setdiff(1:numel(x), I);
                req_x = x(x==x(rep_Num(1)));
                req_y = y(x==x(rep_Num(1)));
                d = abs(req_y(2)-req_y(1));
            factor(1) = 1/(1-d/block_size);
            factor(2) = 1/(d/block_size);
            factor(3) = 1/(2-d/block_size);
            factor(4) = 1/(1+d/block_size);
            fac = 1;
            end

        end
        end
        for ii=1:2
            Txx(Txx==max(max(Txx))) = 0;
        end
        cha = getcha(Txx,0);
        table(i,4)=cha;
        if cha>t
            Txx = fftshift(Txx);
            Txx_r = imrotate(Txx,-angle/2);
            [ht,wt] = size(Txx_r);
            Txx_r = Txx_r(floor(ht/2-block_size/2)+1:floor(ht/2+block_size/2),floor(ht/2-block_size/2)+1:floor(ht/2+block_size/2));
            H_Txx = zeros(block_size);
            H_Txx(block_size/2-1:block_size/2+1,1:end) = 1;
            H_Txx(1:end, block_size/2-1:block_size/2+1) = 1;
            XY_Txx = Txx_r.*H_Txx;
            XY_Txx = fftshift(XY_Txx);
            [x,y]=find(XY_Txx==max(max(XY_Txx)));
            block_sizet = ht;
            if length(x) ~= 1
            x =x(1);
            y = y(1);
            end
            if x==1
                if y > block_size/2
                    y = block_size-y;
                end
                factor(1) = 1/(1-(y-1)/block_size);
                factor(2) = 1/((y-1)/block_size);
                factor(3) = 1/(2-(y-1)/block_size);
                factor(4) = 1/(1+(y-1)/block_size);
                fac = 1;
            elseif y==1
                 if x > block_size/2
                        x = block_size-x;
                 end
                factor(1) = 1/(1-(x-1)/block_size);
                factor(2) = 1/((x-1)/block_size);
                factor(3) = 1/(2-(x-1)/block_size);
                factor(4) = 1/(1+(x-1)/block_size);
                fac = 1;
            end
            factor(factor>2.5)=0;
            factor(factor<0.5)=0;
        end
    end
%     % 获取真实标签
   table(i,1)=fac;
   table(i,2)=ro;
   if length(angle)>1
       angle=angle(4);
   end
   table(i,5)=angle;
   table(i,6:9)=factor;
%     ture = strsplit(name{i},'_');
%     if length(ture)==1
%         ture_fac = 0;
%         ture_ro = 0;
%     elseif length(ture)==2
%         ture2 = str2num(ture{2}(1:end-4));
%         if ture2 > 3
%             ture_ro = 1;
%             ture_angle = ture2;
%         else
%             ture_fac = 1;
%             ture_factor = ture2;
%         end
%     else
%         ture_fac = 1;
%         ture_ro = 1;
%     end
%     
%     %计算准确率
%     if ture_fac==fac && ture_ro == ro
%         totle = totle + 1;
%     end
%     
     if fac == 1
%         fprintf('factor estimated: %f\n', factor');
        fprintf('factor estimated: %f %f %f %f\n', factor(1),factor(2),factor(3),factor(4));
    end
    if ro == 1
        fprintf('angle estimated: %f\n', angle);
    end
    % Padin et al.'s method
    % fprintf('angle estimated using Padin et al.''s method: %f\n', estimatedAngle(4));
end
