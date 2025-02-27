"""
Module with miscellaneous functions. 
"""
import inspect


def print_doc_string(obj: object) -> None:
    """Prints a docstring with all lines justified to the left."""
    info = None
    if obj.__doc__ is not None:
        info = inspect.getdoc(obj)
    if callable(obj):
        signature = inspect.signature(obj)
        if signature.parameters:
            info += (
                "\n\n"
                "Parameter Types & Defaults\n"
                "--------------------------\n"
            )
            n = len(signature.parameters)
            for i, param in enumerate(signature.parameters.values()):
                if param.name != 'self':
                    info += str(param)
                if i < (n - 1):
                    info += '\n'
    print(info)
