import logging
import sys
from acpype.params import tmpLogFile


class LogFormatter(logging.Formatter):

    """
    Define log formatter
    """

    err_fmt = "ERROR: %(msg)s"
    warn_fmt = "WARNING: %(msg)s"
    dbg_fmt = "DEBUG: %(msg)s"
    info_fmt = "%(msg)s"

    def __init__(self):
        super().__init__(fmt="%(levelno)d: %(msg)s", datefmt=None, style="%")

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = LogFormatter.dbg_fmt
        elif record.levelno == logging.INFO:
            self._style._fmt = LogFormatter.info_fmt
        elif record.levelno == logging.ERROR:
            self._style._fmt = LogFormatter.err_fmt
        elif record.levelno == logging.WARNING:
            self._style._fmt = LogFormatter.warn_fmt
        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)
        # Restore the original format configured by the user
        self._style._fmt = format_orig
        return result


def set_logging_conf(level=20):
    # Setting logging configurations
    logging.root.setLevel(0)
    logger = logging.getLogger(__name__)
    if logger.handlers:
        logger.handlers.pop()

    fmt = LogFormatter()
    file_handler = logging.FileHandler(filename=tmpLogFile)
    stdout_handler = logging.StreamHandler(sys.stdout)
    file_handler.setLevel(logging.DEBUG)
    stdout_handler.setLevel(level)
    stdout_handler.setFormatter(fmt)
    file_handler.setFormatter(fmt)
    if not logger.handlers:
        logger.addHandler(file_handler)
    logger.addHandler(stdout_handler)
    # logger.setLevel(logging.DEBUG)
    return logger
