import logging
from logging import StreamHandler, Logger


class ModuleLogger:
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL

    # specify the formatting of the log records
    FORMATTER = logging.Formatter(
        '[%(name)s | %(levelname)s] %(message)s'
    )

    @classmethod
    def _create_console_handler(cls, log_level: int) -> StreamHandler:
        # create a log handler that logs records to the console
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(cls.FORMATTER)
        console_handler.setLevel(log_level)
        return console_handler

    @classmethod
    def get_logger(cls, logger_name: str, log_level: int = logging.DEBUG) -> Logger:
        logger = logging.getLogger(logger_name)
        if not logger.handlers:
            # if a logger with the same `logger_name` is called multiple times,
            # don't add its handlers again
            console_handler = cls._create_console_handler(log_level)
            logger.addHandler(console_handler)
            logger.setLevel(log_level)
        return logger
