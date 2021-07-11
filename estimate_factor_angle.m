function [angle,factor]=estimate_factor_angle(Txx,dev,block_size)
    %标签初始化
    factor = 1;
    h = block_size;
    % 判断是否只经过缩放
    % 获取图像中峰值超过80的至多10个点的坐标
    [x,y] = getdp(Txx,40,h,80);
    if x(1)<2 || y(1)<2
        %% if the brightest peak locate at axis, then no rotation, i.e.,
        % angle = 0
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
        angle = 0;
    else
        %% else estimate first then scaling factor
        % normalize the peak location, 
        % (x-1, y-1) is better than (x,y)
        x2w = (x-1)/block_size;
        y2w = (y-1)/block_size;
        %计算候选点与左上角的距离
        d1 = sqrt(x2w.^2+y2w.^2);
        % we take those pixels whose distance to (1,1) between (1-dev, 1+dev)
        % as the curve for searching peak of rotation
        x_upleft = x(d1<1+dev&d1>1-dev);
        y_upleft = y(d1<1+dev&d1>1-dev);
        %处理左上部分角点频率距离在0.09-1.01之间的备选点
        if ~isempty(x_upleft)
            %选取最可能的旋转候选点
            [x1,y1] = getOnly(x_upleft,y_upleft,Txx);
            %取左上角块
            Txx_r = Txx(1:block_size/2,1:block_size/2);
            %计算可能的旋转角度
            angle = getangle(x1,y1,Txx,1,block_size);
            %获取左上角块中峰值大于50的至多5个缩放候选点
            % [x2,y2] = getdp(Txx_r,5,h,50);
            [x2,y2] = getdp_excludeaxis(Txx_r, 3, h/2, 50);
            %判断旋转候选点与可能的缩放候选点间是否存在角度一致
            a = length(x_upleft);
            temp_angle = angle;
            while(a>1)
                % 计算Txx_r中峰值对应角度和当前angle 差异
                angle_temp = abs(atan((y2-1)./(x2-1)).*180./pi-angle);
                if ~(min(angle_temp)<1)
                    % 如果在左上区域不能找到对应两点，去掉候选角度
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
            %不存在旋转候选点与缩放候选点的关系，认为只经过旋转，旋转角度使用最大旋转峰值的结果
            if ~(min(angle_temp)<1)
                % Almost fail if reach here
                angle = temp_angle;
            end
            %计算缩放因子
            factor = getfactor(angle,Txx_r,block_size); 
        else
            % No rescaling or rotation
            angle = 0;
            factor = 1;
        end
    end
end

function factor = getfactor(angle,Txx,block_size)
    %通过旋转角度、频域和块计算缩放因子
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

% get the brightest peak from the candidates
function [only_x,only_y]=getOnly(x,y,Txx)
%通过峰值大小选择可能的候选点
    if min(length(x),length(y))>1
        max_x = x(1);
        max_y = y(1);
        for j = 2: min(length(x),length(y))
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
%计算旋转角度
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
%计算Txx中大于threshold的至多num个值
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

%%% 20210702
function [x,y] = getdp_excludeaxis(Txx,num,h,threshold)
    x=0;
    y=0;
    meshXY = meshgrid(1:h, 1:h);
    meshYX = meshXY';
    R2 = (meshXY.^2+meshYX.^2).^0.5;
    interest = (R2 < h);
    interest(1,:)=0;
    interest(:,1)=0;
    %%% have a basic idea of interest region
    % figure, imshow(interest)
    %%%
    Txx = Txx.*interest;
    %计算Txx中大于threshold的至多num个值
    i = 1;
    while(i<=num)
       [maxx,maxy] = find(Txx == max(Txx(:)));
       if max(Txx(:)) < threshold
           break;
       end
       if length(maxx)~=1
            x(i:i+length(maxx)-1) = maxx;
            y(i:i+length(maxx)-1) = maxy;
            i = i + length(maxx);  
            for j = 1:length(maxx)
                Txx(maxx(j)-1: maxx(j)+1, maxy(j)-1: maxy(j)+1)=0;
            end
       else
            x(i) = maxx;
            y(i) = maxy;
            i = i + 1;
            Txx(maxx-1:maxx+1,maxy-1:maxy+1)=0;
       end
    end
%     if(isempty(x)||isempty(y))
%         x=0;
%         y=0;
%     end
end