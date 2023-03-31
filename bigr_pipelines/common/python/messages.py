#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file contains pretty-print functions for logging
"""

import logging


def message(tag: str = "info", message: str = ""):
    """
    Add a colored tag before a printed message
    """
    colors = {
        "info": "\033[1;36m@INFO:\033[0m)",
        "cmd": "\033[1;32m@CMD:\033[0m",
        "error": "\033[41m@ERROR:\033[0m",
        "doc": "\033[0;33m@DOC:\033[0m",
        "warning": "\033[1;33m@WARNING:\033[0m",
    }
    print(colors.get(tag.lower(), "") + message)


class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""

    grey = "\x1b[38;21m"
    yellow = "\x1b[33;21m"
    red = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "(%(filename)s:%(lineno)d) - %(name)s - %(levelname)s - %(message)s "

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
