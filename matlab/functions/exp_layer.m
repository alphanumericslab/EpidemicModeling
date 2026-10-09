classdef exp_layer < nnet.layer.Layer
    % EXP_LAYER Learnable exponential activation: Z = exp(Alpha .* X).
    % layer = exp_layer(num_channels,name). Requires Deep Learning Toolbox.
    % Alpha initializes with shape 1-by-1-by-num_channels. Assign supplied
    % weights explicitly for MATLAB/Python parity; random streams differ.
    % Author: Reza Sameni | Emory University
    % Reference: Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
    properties (Learnable)
        Alpha
    end
    methods
        function layer = exp_layer(num_channels, name)
            % EXP_LAYER Construct the activation with learnable channel scales.
            layer.Name = name;
            layer.Description = "Exponential activation with " + num_channels + " channels";
            layer.Alpha = randn([1 1 num_channels]);
        end
        function z = predict(layer, x)
            % PREDICT Evaluate the elementwise exponential forward pass.
            z = exp(layer.Alpha .* x);
        end
    end
end
