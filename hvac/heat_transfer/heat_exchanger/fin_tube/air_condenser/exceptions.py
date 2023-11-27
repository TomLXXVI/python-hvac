class CondenserError(Exception):
    pass


class DesuperheatingError(CondenserError):
    pass


class CondensingError(CondenserError):
    pass


class SubcoolingError(CondenserError):
    pass
