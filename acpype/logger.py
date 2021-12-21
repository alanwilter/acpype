import logging


class MyFileFormatter(logging.Formatter):

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
            self._style._fmt = MyFileFormatter.dbg_fmt

        elif record.levelno == logging.INFO:
            self._style._fmt = MyFileFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._style._fmt = MyFileFormatter.err_fmt

        elif record.levelno == logging.WARNING:
            self._style._fmt = MyFileFormatter.warn_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result


def set_logging_conf():
    # Setting logging configurations
    fmt = MyFileFormatter()
    file_handler = logging.FileHandler(filename="/tmp/acpype_run.log")
    file_handler.setFormatter(fmt)
    handlers = [file_handler]
    logging.basicConfig(
        level=logging.DEBUG, handlers=handlers,
    )
