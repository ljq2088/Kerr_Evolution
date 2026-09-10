class ContractError(ValueError):
    """An input violates a frozen C0 formula or array contract."""


class OpenBoundaryError(NotImplementedError):
    """A calculation requires a boundary closure left open for V1.2."""

