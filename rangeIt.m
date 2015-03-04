function E = rangeIt(en)
%%rangeIt(en) returns a new image with values in the range 0-255
%%Input Arguments:
%%en: Input image
%%Output:
%%E: Ranged Image
    
    %%Minimum value in the input image
    mi = min(min(en));
    
    %%If the minimum value of negative, make all the values positive
    if  mi < 0
        en = en - mi;
    end
    
    %%Maximum value in the input image
    ma = max(max(en));
    
    %%Normalize and scale all the values from 0-255
    E = (en./ma)*255;        
end