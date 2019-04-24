function cv = clip(v, bds)
    if(nargin<2)
        bds = [0,1];
    end
	cv = v;
    cv(v(:)<bds(1)) = bds(1);
    cv(v(:)>bds(2)) = bds(2);
end