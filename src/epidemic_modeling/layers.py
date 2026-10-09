"""NumPy and optional PyTorch ports of the original learnable activations.

Author: Reza Sameni, Emory University. Stored weights are required for
cross-language parity because MATLAB/NumPy/PyTorch random streams differ.
"""

import numpy as np


def exp_layer(x, alpha):
    """Evaluate the learnable exponential activation ``exp(alpha*x)``.

    x and alpha must be broadcast-compatible; supply MATLAB Alpha weights
    explicitly. The optional torch_exp_layer factory adds automatic gradients.

    Returns a NumPy array with the broadcast shape of x and alpha. Large
    positive alpha*x can overflow; choose scales appropriate for the input.
    """

    return np.exp(np.asarray(alpha) * np.asarray(x))


def my_tanh_layer(x, alpha):
    """Evaluate ``alpha*tanh(x/alpha)`` with broadcast-compatible nonzero weights.

    Pass the original MATLAB Alpha array to reproduce an existing layer.
    A zero scale is rejected because the historical expression is undefined.

    Returns a NumPy array with the broadcast shape of x and alpha.
    Each output approaches +/-abs(alpha) as the input magnitude increases.
    """
    alpha = np.asarray(alpha)

    if np.any(alpha == 0):
        raise ValueError("alpha scales must be nonzero")

    return alpha * np.tanh(np.asarray(x) / alpha)


def torch_exp_layer(alpha):
    """Create a trainable PyTorch module initialized with explicit exponential weights.

    Requires the optional ``neural`` dependency. Input tensor and alpha must
    be broadcast-compatible. No GPU, network, or random initialization is used.
    """
    import torch

    class exponential_layer(torch.nn.Module):
        """Trainable exponential activation initialized from supplied weights."""

        def __init__(self, weights):
            """Register the learnable scale tensor."""
            super().__init__()
            self.alpha = torch.nn.Parameter(
                torch.as_tensor(weights, dtype=torch.float64).clone()
            )

        def forward(self, x):
            """Evaluate the activation with automatic differentiation."""

            return torch.exp(self.alpha * x)

    return exponential_layer(alpha)


def torch_my_tanh_layer(alpha):
    """Create a trainable PyTorch scaled-tanh module with explicit nonzero weights.

    Requires the optional neural dependency. alpha supplies the initial
    nonzero scale weights. Returns a torch.nn.Module with float64 learnable
    weights; inputs must have compatible shape, dtype, and device.
    """
    import torch

    if np.any(np.asarray(alpha) == 0):
        raise ValueError("alpha scales must be nonzero")

    class scaled_tanh_layer(torch.nn.Module):
        """Trainable scaled hyperbolic tangent initialized from supplied weights."""

        def __init__(self, weights):
            """Register the learnable scale tensor."""
            super().__init__()
            self.alpha = torch.nn.Parameter(
                torch.as_tensor(weights, dtype=torch.float64).clone()
            )

        def forward(self, x):
            """Evaluate the activation with automatic differentiation."""

            return self.alpha * torch.tanh(x / self.alpha)

    return scaled_tanh_layer(alpha)
