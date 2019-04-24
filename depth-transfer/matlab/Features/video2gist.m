function gist = video2gist( V )
%VIDEOGIST Computes gist after temporal median filtering the video sequence
% V is an NxMx{1,3}xK video sequence (K=#frames)
    if(nargin==0) %Return number of features if no args given
        gist = im2gist();
        return;
    end
    %gist = im2gist(V(:,:,:,1)); %Using just a single frame might suffice
    gist = im2gist(median(V,4));
end

