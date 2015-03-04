function getBoundary(Initial_3DMaskWM, Initial_3DMaskGM, InputImages_3D)
%load ~/Downloads/CS269_PROJECT/Snakes' Input - MAT files'/SnakesInput_2.mat
inp = InputImages_3D(:,:,60);
inp = rangeIt(inp);
inp = uint8(inp);

im2 = Initial_3DMaskWM(1,:,:,60);
im = zeros(size(inp));
im(:,:) = im2(1,:,:);
ed = edge(im);
ed = imdilate( ed, strel('line', 2,0) );
[x,y] = find(ed == 1);

dx = 0;
dy = 0;
for i=-1:1
    for j=-1:1
        if(ed(x(1)+i,y(1)+j) == 1)
            dx = i;
            dy = j;
        end
    end
end
dir = '';
if(dx == 1)
    if(dy == 1)
        dir = 'SE';
    elseif(dy == -1)
        dir = 'NE';
    else
        dir = 'E';
    end
elseif(dx == 0)
    if(dy == 1)
        dir = 'N';
    elseif(dy == -1)
        dir = 'S';
    end
else
    if(dy == 1)
        dir = 'SW';
    elseif(dy == -1)
        dir = 'NW';
    else
        dir = 'W';
    end
end

P = bwtraceboundary(ed, [x(1) y(1)],dir,8);

x = P(:,1);
y = P(:,2);
flag = zeros(size(inp));
for i=1:length(x)
    flag(x(i),y(i)) = 1;
end

x1 = x(1);
y1 = y(1);

xx = x(2);
yy = y(2);

Q = [x1 y1];

breakFlag = 0;
visited = zeros(size(flag));
visited(x1,y1) = 1;
flag(357,95) = 1;
[sx sy] = size(flag);
%disp = zeros(sx,sy,3);
%disp(:,:,1) = flag;
%imshow(disp);hold on;
%pause;
while(xx ~= x1 || yy ~= y1)
    breakFlag = 0;
    for j=1:-1:-1
        for i=1:-1:-1
            if((flag(xx+i,yy+j) == 1 && visited(xx+i,yy+j) == 0) || (xx+i == x1 && yy+j == y1))
                visited(xx+i,yy+j) = 1;
                %  disp(xx+i,yy+j,1) = 0;
                %  disp(xx+i,yy+j,2) = 1;
                %  imshow(disp);
                %plot(yy+j,xx+i);
                %pause;
                Q = [Q;[xx,yy]];
                xx = xx + i;
                yy = yy + j;
                breakFlag = 1;
                break;
            end
        end
        if(breakFlag == 1)
            break;
        elseif(breakFlag == 0 && j == -1)
            xx = Q(end,1);
            yy = Q(end,2);
        end
    end
end

% % tmpim = zeros(size(flag));
% % x = Q(:,1);
% % y = Q(:,2);
% % for i=1:length(x)
% %     tmpim(x(i),y(i)) = 1;
% % end
% % imshow(tmpim);hold on;
snake(inp,Q(:,2),Q(:,1),50,1);
end
