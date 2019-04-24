function anaglyphImg = makeAnaglyph( L_img, R_img, type )
%MAKEANAGLYPH Converts a stereo pair into anaglyph. See website below
%             for available types.
%             http://www.3dtv.at/Knowhow/AnaglyphComparison_en.aspx

    % Assumes side-by-side pair if R_img is not provided or empty.
    if(nargin<2 || isempty(R_img) )
        R_img = L_img(:,ceil(end/2)+1:end,:);
        L_img = L_img(:,1:floor(end/2),:);
    end

    %Error checking and type conversion
    [h,w,d] = size(L_img);
    assert( all([h,w,d]==size(R_img)), 'Left and right images must be the same dimension.');
    if(nargin<3)
        type = 'optimized';
    end
    if( ~isfloat(L_img) )
        L_img = double(L_img) ./ 255;
    end
    if( ~isfloat(R_img) )
        R_img = double(R_img) ./ 255;
    end
    
    %Create conversion transformation
    if( strcmp(type,'true') )
        L_mat = [0.299, 0.587, 0.114; 0,0,0; 0,0,0];
        R_mat = [0,0,0; 0,0,0; 0.299, 0.587, 0.114];
    elseif( strcmp(type,'gray') )
        L_mat = [0.299, 0.587, 0.114; 0,0,0; 0,0,0];
        R_mat = [0,0,0; 0.299, 0.587, 0.114; 0.299, 0.587, 0.114];
    elseif( strcmp(type,'color') )
        L_mat = [1,0,0; 0,0,0; 0,0,0];
        R_mat = [0,0,0; 0,1,0; 0,0,1];
    elseif( strcmp(type,'halfcolor') )
        L_mat = [0.299, 0.587, 0.114; 0,0,0; 0,0,0];
        R_mat = [0,0,0; 0,1,0; 0,0,1];
    elseif( strcmp(type,'optimized') )
        L_mat = [0, 0.7, 0.3; 0,0,0; 0,0,0];
        R_mat = [0,0,0; 0,1,0; 0,0,1];
    else
        error('Unreconized type. Available types: true, gray, color, halfcolor, optimized.\n');
    end
    
    %Make anaglyph
    anaglyphImg = anaglyphConvert(L_img, L_mat, R_img, R_mat);
    
    %Apply gamma correction to red channel?
    anaglyphImg(:,:,1) = anaglyphImg(:,:,1).^(1/1.5);
end

function anaglyph = anaglyphConvert(L_img, Lm, R_img, Rm)
    h = size(L_img,1);
    w = size(L_img,2);
    
    L_mult_r = cat(3, Lm(1,1)*ones(h,w), Lm(1,2)*ones(h,w), Lm(1,3)*ones(h,w));
    L_mult_g = cat(3, Lm(2,1)*ones(h,w), Lm(2,2)*ones(h,w), Lm(2,3)*ones(h,w));
    L_mult_b = cat(3, Lm(3,1)*ones(h,w), Lm(3,2)*ones(h,w), Lm(3,3)*ones(h,w));
    
    R_mult_r = cat(3, Rm(1,1)*ones(h,w), Rm(1,2)*ones(h,w), Rm(1,3)*ones(h,w));
    R_mult_g = cat(3, Rm(2,1)*ones(h,w), Rm(2,2)*ones(h,w), Rm(2,3)*ones(h,w));
    R_mult_b = cat(3, Rm(3,1)*ones(h,w), Rm(3,2)*ones(h,w), Rm(3,3)*ones(h,w));
    
    anaglyph = cat(3, sum(L_img.*L_mult_r,3) + sum(R_img.*R_mult_r,3), ...
                      sum(L_img.*L_mult_g,3) + sum(R_img.*R_mult_g,3), ...
                      sum(L_img.*L_mult_b,3) + sum(R_img.*R_mult_b,3) );
end