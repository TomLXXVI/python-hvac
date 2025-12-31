import logging
from logging import StreamHandler, FileHandler, Logger
from logging.handlers import QueueHandler
from multiprocessing import Queue
from pathlib import Path


class ModuleLogger:
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL

    # specify the formatting of the log records
    FORMATTER = logging.Formatter(
        '[%(process)s | %(name)s | %(levelname)s] %(message)s'
    )

    @classmethod
    def create_console_handler(cls, log_level: int | None = None) -> StreamHandler:
        # Create a log handler that logs records to the console
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(cls.FORMATTER)
        console_handler.setLevel(log_level or cls.DEBUG)
        return console_handler

    @classmethod
    def create_file_handler(cls, file_path: Path | str, log_level: int | None = None) -> FileHandler:
        # Create a log handler that logs records to a file.
        file_handler = logging.FileHandler(file_path, mode='a', encoding='utf-8')
        file_handler.setFormatter(cls.FORMATTER)
        file_handler.setLevel(log_level or cls.DEBUG)
        return file_handler

    @classmethod
    def create_queue_handler(cls, queue: Queue) -> QueueHandler:
        queue_handler = QueueHandler(queue)
        return queue_handler

    @classmethod
    def get_logger(
        cls,
        logger_name: str,
        file_path: Path | str | None = None,
        log_level: int = logging.DEBUG
    ) -> Logger:
        """Returns a logger with `logger_name`. If a logger with `logger_name`
        already exists, this same logger will be returned. Otherwise, a new
        logger is returned. If a new logger and `file_path` is not None, a
        logger will be returned with a console and file handler, that will log
        the same messages both to the console and to the file. Only messages
        with a priority equal or higher than `log_level` will be logged.
        """
        logger = logging.getLogger(logger_name)
        if not logger.handlers:
            # if a logger with the same `logger_name` is called multiple times,
            # don't add its handlers again
            console_handler = cls.create_console_handler(log_level)
            logger.addHandler(console_handler)
            if file_path is not None:
                file_handler = cls.create_file_handler(file_path, log_level)
                logger.addHandler(file_handler)
            logger.setLevel(log_level)
        return logger

    @classmethod
    def get_logger_to_queue(
        cls,
        logger_name: str,
        queue: Queue,
        log_level: int = logging.DEBUG,
    ) -> Logger:
        """Returns a logger that puts log messages on a queue."""
        logger = logging.getLogger(logger_name)
        if not logger.handlers:
            queue_handler = cls.create_queue_handler(queue)
            logger.addHandler(queue_handler)
            logger.setLevel(log_level)
        return logger

    @staticmethod
    def sort_log_by_process_id(file_path: Path | str) -> list[str]:
        """Sorts the log entries in a log file at `file_path` by the process-ID
        mentioned in each log entry and saves the sorted file back to disk
        with the suffix '_sorted' added to the original filename.
        """
        log_entries = []
        traceback_lines = []
        with open(file_path, 'r', encoding='utf-8') as fh:
            for line in fh.readlines():
                if line.startswith('['):
                    # We have a new log entry.
                    # If previous lines came from a traceback, this traceback
                    # belongs to the previous log entry. In that case, we first
                    # add the traceback to the last log entry in `log_entries`:
                    if traceback_lines:
                        traceback_str = ''.join(traceback_lines)
                        log_entries[-1] += traceback_str
                        traceback_lines = []
                    # Append the new log entry to `log_entries`:
                    log_entries.append(line)
                else:
                    # Lines that don't start with '[' belong to a traceback.
                    # We collect the traceback lines in a separate list that
                    # we will add to the last log entry, when a new log entry is
                    # encountered:
                    traceback_lines.append(line)

        log_entries = [entry.strip() for entry in log_entries]

        def _get_pid(entry: str) -> str:
            *_, after = entry.partition('[')
            before, *_ = after.partition(']')
            pid = before.split('|')[0]
            pid = pid.strip()
            return pid

        sorted_log_entries = sorted(
            log_entries,
            key=lambda entry: _get_pid(entry)
        )

        if isinstance(file_path, str):
            file_path = Path(file_path)

        dir_path = file_path.parent
        file_name = file_path.stem
        extension = file_path.suffix
        new_filepath = dir_path / Path(file_name + "_sorted" + extension)

        with open(new_filepath, 'w', encoding='utf-8') as fh:
            fh.writelines(line + '\n' for line in sorted_log_entries)

        return sorted_log_entries
