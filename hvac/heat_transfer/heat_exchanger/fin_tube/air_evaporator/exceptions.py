class EvaporatorError(Exception):
    pass


class SuperheatingError(EvaporatorError):
    pass


class BoilingError(EvaporatorError):
    pass

