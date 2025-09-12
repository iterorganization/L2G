import logging
import sys
import os

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

old_factory = logging.getLogRecordFactory()

def custom_log_record_factory(*args, **kwargs) -> logging.LogRecord:
    """Custom log record in order to correctly obtain information which module
    has called logging.
    """
    record = old_factory(*args, **kwargs)
    module_name = record.name

    # Check if the module is a Cython module
    module = sys.modules.get(module_name)
    if module:
        # If the module has a `__file__` attribute, use it to resolve the filename
        if hasattr(module, '__file__') and not module.__file__ is None:
            file_name: str = module.__file__

            if file_name.endswith("so"):
                # Extract the name of the pyx file by
                basename = os.path.basename(file_name)
                split_basename = basename.split(".")
                if len(split_basename) <= 1:
                    pass
                else:
                    file_name = split_basename[0]
            record.filename = file_name
        else:
            # Fallback for Cython modules without `__file__`. The name of the
            # *.so files contains cpython and other identifiers so only include
            # the string up to the first '.'
            record.filename = module_name + ".pyx"
    return record

logging.setLogRecordFactory(custom_log_record_factory)

_file_handler = None
_stream_handler = None

fmt = '%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] : %(message)s'
fmt_date = '%Y-%m-%d %T'

def addFileHandler(output_file: str='l2g.log') -> None:
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

def addStreamHandler() -> None:
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

def enableLogging() -> None:
    """Raises the l2g logger level to INFO.
    """
    global log
    log.setLevel(logging.INFO)

def enableDebugging() -> None:
    """Raises the l2g logger level to DEBUG.
    """
    global log
    log.setLevel(logging.DEBUG)