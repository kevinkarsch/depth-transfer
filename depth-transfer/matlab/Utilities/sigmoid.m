function sx = sigmoid(x, center, scale)
%SIGMOID passes x through a sigmoidal function with desired center and scale
    sx = 1./(1+exp(scale.*(center-x)));
end
