"""
This file contains pretty-print functions for logging
"""

def message(tag: str = "info", message: str = "", /):
    """
    Add a colored tag before a printed message
    """
    colors = {
        "info": "\033[1;36m@INFO:\033[0m)",
        "cmd": "\033[1;32m@CMD:\033[0m",
        "error": "\033[41m@ERROR:\033[0m",
        "doc": "\033[0;33m@DOC:\033[0m",
        "warning": "\033[1;33m@WARNING:\033[0m"
    }
    print(colors.get(tag.lower(), "") + message)
