classdef my_tanh_layer < nnet.layer.Layer
    % MY_TANH_LAYER Learnable scaled tanh: Z = Alpha .* tanh(X ./ Alpha).
    % layer = my_tanh_layer(num_channels,name,initial_param_std).
    % Requires Deep Learning Toolbox. The original square num_channels-by-
    % num_channels scale initialization is preserved. Assign explicit nonzero
    % Alpha weights for reproducible cross-language forward passes.
    % Author: Reza Sameni | Emory University
    % Reference: Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
    properties (Learnable)
        Alpha
    end

    methods

        function layer = my_tanh_layer(num_channels, name, initial_param_std)
            % MY_TANH_LAYER Construct a scaled tanh with learnable parameters.
            layer.Name = name;
            layer.Description = "Scaled tanh with " + num_channels + " channels";
            layer.Alpha = initial_param_std * randn(num_channels);
        end

        function z = predict(layer, x)
            % PREDICT Evaluate the scaled-tanh forward pass.
            z = layer.Alpha .* tanh(x ./ layer.Alpha);
        end

    end
end
