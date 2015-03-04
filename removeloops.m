function [x,y] = removeloops(x,y,im)
%%removeloops(x,y,im) returns the new coordinates of the snake after
%%removing the loops in the snake
%%Input Arguments:
%%x: X coordinates of the contour nodes
%%y: Y coordinates of the contour nodes
%%im: Input image
%%Output:
%%x: X coordinates of the contour nodes after removing the loop
%%y: Y coordinates of the contour nodes after removing the loop

%%Get the integer coordinates for the snake contour
x = round(x);
y = round(y);

for i=1:length(x)
    if x(i)<=0
        x(i) = 1;
    elseif x(i) > size(im,2)
        x(i) = size(im,2);
    end
    
    if y(i)<=0
        y(i) = 1;
    elseif y(i) > size(im,1)
        y(i) = size(im,1);
    end
end

%%Create a matrix to check if more than one contour node appears on a
%%distinct cell
nodes = zeros(size(im));

%%Maintains the intersection points
intersections = [];

for i=1:length(x)
    %%Check if node already has a point on it
    %%If yes then add it to the intersection points
    if nodes(y(i),x(i)) ~= 0
        intersections = [intersections;[nodes(y(i),x(i)),i]];
    end
    %%Update the node with new point
    nodes(y(i),x(i)) = i;  
    %nodes = mark(nodes,x(i),y(i),i);
end

%%logical array to check which points to keep and which to remove
select = ones(length(x),1);
%%Set the logical array
for i=1:size(intersections,1)
    select(intersections(i,1)+1:intersections(i,2)) = 0;
end
%%Convert to logical
select = logical(select);
%%Select from the X and Y coordinates
x = x(select);
y = y(select);

end

function nodes = mark(nodes,x,y,val)
%%mark(nodes,x,y,val) returns the nodes matrix after marking the
%%9-neighbors of the point

[Xsize,Ysize] = size(nodes);

for i=-1:1
    for j=-1:1
        if(x+i <= Xsize && y+j <= Ysize && x+i >= 1 && y+j >= 1)
            nodes(y+j,x+i) = val;
        end
    end
end

end