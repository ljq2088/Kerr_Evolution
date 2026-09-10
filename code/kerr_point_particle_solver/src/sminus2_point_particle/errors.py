from kerr_waveform_tools.contracts import ContractError


class OpenBoundaryError(NotImplementedError):
    pass


__all__ = ["ContractError", "OpenBoundaryError"]
