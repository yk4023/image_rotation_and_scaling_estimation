function [angle,factor]=estimate_factor_and_angle(Txx)
% -------------------------------------------------------------------------
%   Written Feb 10, 2021 by Kun Yu
%   Copyright 2021 by Kun Yu
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The end-user understands that the program was
% developed for research purposes and is advised not to rely exclusively on
% the program for any reason.
% -------------------------------------------------------------------------
% Contact: yukuncomputer@foxmail.com
% -------------------------------------------------------------------------
% Input:       Txx .... Cyclic correlation spectrum produced by
%                       function CalculateTxx
% Output:    angle .... the rotation angle of image
%           factor .... the scaling factor of image
% -------------------------------------------------------------------------
% Label initialization and get the size of Txx
    factor = 1;
    [block_size,~] = size(Txx);
    dev = 0.01;
% Obtain the coordinates of up to 40 points in the image whose peak value exceeds 80
    [x,y] = getdp(Txx,40,block_size,80);
    if x(1)<2 || y(1)<2
        if x(1)<2 && y(1)<block_size/2
            d = y(1)/block_size;
            factor = 1 / (1-d);
        elseif x(1)<2 && y(1)>block_size/2
            d = (1-y(1)/block_size);
            factor = 1 / (1-d);
        elseif y(1)<2 && x(1)>block_size/2
            d = (1-x(1)/block_size);
            factor = 1 / (1-d);
        elseif y(1)<2 && x(1)<block_size/2
            d = x(1)/block_size;
            factor = 1 / (1-d);
        end
    end
% Normalize the candidate points
    x2w = x/block_size;
    y2w = y/block_size;
% Calculate the distance between the candidate point and the upper left
% corner
    d1 = sqrt(x2w.^2+y2w.^2);
% The point whose frequency distance from the upper left corner is between
% 0.09 and 1.01 is selected as the candidate point.
    x_upleft = x(d1<1+dev&d1>1-dev);
    y_upleft = y(d1<1+dev&d1>1-dev);
% The candidate points with frequency distance between 0.09 and 1.01 in 
% the upper left corner were processed
    if ~isempty(x_upleft)
    % Select the most likely rotation candidate
        [x1,y1] = getOnly(x_upleft,y_upleft,Txx);
    % Get P_2
        P_2 = Txx(1:block_size/2,1:block_size/2);
    % Calculate the possible rotation angle
        angle = getangle(x1,y1,Txx,1,block_size);
    % Get at most 5 scaling candidate points with peak value greater than
    % 50 in the P_2
        [x2,y2] = getdp(P_2,5,block_size,50);
    % Judge whether the angle between the rotation candidate and the  
    % possible scaling candidate is consistent
        a = length(x_upleft);
        temp_angle = angle;
        while(a>1)
            angle_temp = abs(atan((y2-1)./(x2-1)).*180./pi-angle);
            if ~(min(angle_temp)<1)
                x_upleft(x_upleft==x1)='';
                y_upleft(y_upleft==y1)='';
                a = length(x_upleft);
                [x1,y1] = getOnly(x_upleft,y_upleft,Txx);
                if ~isempty(x1) && ~isempty(y1)
                    angle = getangle(x1,y1,Txx,1,block_size);
                end
            else
                a=0;
            end
        end
        angle_temp = abs(atan((y2-1)./(x2-1)).*180./pi-angle);
        % There is no relationship between rotation candidate points and 
        % scaling candidate points. It is considered that only rotation is 
        % used and the rotation angle uses the maximum rotation peak
        if ~(min(angle_temp)<1)
            angle = temp_angle;
        end
        % Calculate scaling factor 
        if factor==1
            factor = getfactor(angle,P_2,block_size); 
        end
    else
        angle = 0;
    end
end
function factor = getfactor(angle,Txx,block_size)
    % The scaling factor is calculated by rotation angle, frequency domain 
    % and block
    factor = 1;
    h = block_size/2; 
    [x,y] = getdp(Txx,10,h,50);
    x1 = x(sqrt(x.^2+y.^2)<h);
    y1 = y(sqrt(x.^2+y.^2)<h);
    if length(x)<3
        factor = 1;
    else
        angleXY = abs(atan((y1-1)./(x1-1)).*180./pi-angle);
        if min(angleXY)<1
            need_x = x1(angleXY==min(angleXY));
            need_y = y1(angleXY==min(angleXY));
            if length(need_x)>1
                eng = zeros(1,length(need_x));
                for i = 1:length(need_x)
                    eng(i) = Txx(need_x(i),need_y(i));
                end
                need_x = need_x(eng==max(eng));
                need_y = need_y(eng==max(eng));
                while(abs(sqrt((1-need_x/256)^2+(1-need_y/256)^2)-1)<0.01)
                    eng(eng==max(eng))=0;
                    need_x = x1(angleXY==min(angleXY));
                    need_y = y1(angleXY==min(angleXY));
                    need_x = need_x(eng==max(eng));
                    need_y = need_y(eng==max(eng));
                end
            end
            d = sqrt(need_x(1)^2+need_y(1)^2)/block_size;
            factor = 1/(1-d);
        end
    end
end
function [only_x,only_y]=getOnly(x,y,Txx)
% The possible candidate points are selected by the peak size
    if length(x)>1 && length(y)>1
        max_x = x(1);
        max_y = y(1);
        for j = 2:min(length(x),length(y))
            if Txx(x(j),y(j))>Txx(max_x,max_y)
                max_x = x(j);
                max_y = y(j);
            end
        end
        only_x = max_x;
        only_y = max_y;
    else
        only_x = x;
        only_y = y;
        
    end
end
function angle = getangle(x,y,Txx,type,block_size)
% Calculate rotation angle
    max_point = [x(1),y(1)];
    if length(x)>1 && length(y)>1
        for m = 2:length(x)
            if Txx(x(m),y(m))>Txx(max_point(1),max_point(2))
                 max_point = [x(m),y(m)];
            end
        end
    end
    if type==1
        angle = atan(max_point(2)/max_point(1))*180/pi;
    elseif type==2
        angle = atan((block_size-max_point(2))/max_point(1))*180/pi;
    elseif type==3
        angle = atan(max_point(2)/(block_size-max_point(1)))*180/pi;
    elseif type==4
        angle = atan((block_size-max_point(2))/(block_size-max_point(1)))*180/pi;
    end
end
function [x,y] = getdp(Txx,num,h,threshold)
% Calculate at most num values greater than threshold in Txx
    i = 1;
    while(i<=num)
       [maxx,maxy] = find(Txx == max(max(Txx)));
       if length(maxx)~=1
            x(i:i+length(maxx)-1) = maxx;
            y(i:i+length(maxx)-1) = maxy;
            i = i + length(maxx);  
            for j = 1:length(maxx)
                if sum(maxx~=1)&&sum(maxx~=h)&&sum(maxy~=1)&&sum(maxy~=h)
                    Txx(maxx(j)-1+(maxx(j)==1):maxx(j)+1-(maxx(j)==h),maxy(j)-1+(maxy(j)==1):maxy(j)+1-(maxy(j)==h))=0;
                else
                    Txx(maxx(j):maxx(j),maxy(j):maxy(j))=0;
                end
            end
       else
            x(i) = maxx;
            y(i) = maxy;
            i = i + 1;
            if maxx~=1&&maxx~=h&&maxy~=1&&maxy~=h
                Txx(maxx-1:maxx+1,maxy-1:maxy+1)=0;
            else
                Txx(maxx:maxx,maxy:maxy)=0;
            end
       end
       if max(max(Txx)) < threshold
           break;
       end
    end
end
