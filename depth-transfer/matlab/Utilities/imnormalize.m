function [ I_norm ] = imnormalize( img, low_prc, up_prc )
%IMNORMALIZE normalize an ND vector with some error checking
	I = img;
    if ~isfloat(I)
        I = double(I);
    end
	if( abs(max(I(:))-min(I(:)))<1e-8 )
        I_norm = img;
    else
        if(nargin==1)
            I_norm = (I - min(I(:))) / (max(I(:)) - min(I(:)));
        else
            m = prctile(I(:),low_prc);
            M = prctile(I(:),up_prc);
            I_norm = (I-m)/(M-m);
        end
	end
end

