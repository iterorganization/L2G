import logging

log = logging.getLogger(__name__)
_file_handler = None
_stream_handler = None

fmt = '%(asctime)s %(levelname)-8s %(name)27s : %(message)s'
fmt_date = '%Y-%m-%d %T'


def addFileHandler(output_file='l2g.log'):
    """Add a handler to redirect logs to a file.

    Arguments:
        output_file (str): Location of the log file.

    """
    global _file_handler, log
    if _file_handler is not None:
        return

    _file_handler = logging.FileHandler(output_file, mode='a')
    _file_handler.setFormatter(logging.Formatter(fmt, fmt_date))
    log.addHandler(_file_handler)

def addStreamHandler():
    """Adds a stream handler for redirecting the log output to the terminal.
    There is no streamhandler added by default, for better integration into
    GUI.
    """
    global _stream_handler, log
    if _stream_handler is not None:
        return
    _stream_handler = logging.StreamHandler()
    _stream_handler.setFormatter(logging.Formatter(fmt, fmt_date))
    log.addHandler(_stream_handler)

def enableLogging():
    """Raises the l2g logger level to INFO.
    """
    global log
    log.setLevel(logging.INFO)

def enableDebugging():
    """Raises the l2g logger level to DEBUG.
    """
    global log
    log.setLevel(logging.DEBUG)