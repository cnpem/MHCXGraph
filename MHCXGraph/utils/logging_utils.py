from __future__ import annotations

import logging
from pathlib import Path
import coloredlogs

_LOG: "VerboseLoggerAdapter | None" = None


class ConsoleFilter(logging.Filter):
    def __init__(self, *, debug: bool, verbose: bool):
        super().__init__()
        self.debug = debug
        self.verbose = verbose

    def filter(self, record: logging.LogRecord) -> bool:
        if record.levelno == logging.DEBUG and not self.debug:
            return False
        if getattr(record, "verbose_only", False) and not self.verbose:
            return False
        return True


class VerboseLoggerAdapter(logging.LoggerAdapter):
    def vinfo(self, msg, *args, **kwargs):
        extra = dict(kwargs.pop("extra", {}))
        extra["verbose_only"] = True
        kwargs["extra"] = extra
        self.logger.info(msg, *args, **kwargs)


def setup_logging(
    *,
    outdir: Path,
    debug: bool,
    verbose: bool,
    log_file_name: str = "debug.log",
) -> VerboseLoggerAdapter:
    global _LOG

    root = logging.getLogger("MHCXGraph")
    root.handlers.clear()
    root.propagate = False
    root.setLevel(logging.DEBUG)

    outdir.mkdir(parents=True, exist_ok=True)
    file_path = outdir / log_file_name

    console_fmt = "%(asctime)s %(message)s"
    console_date = "%Y-%m-%d %H:%M:%S"

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG if debug else logging.INFO)
    ch.addFilter(ConsoleFilter(debug=debug, verbose=verbose))
    ch.setFormatter(coloredlogs.ColoredFormatter(fmt=console_fmt, datefmt=console_date))
    root.addHandler(ch)

    fh = logging.FileHandler(file_path, mode="a", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s %(filename)s:%(lineno)d | %(message)s",
            "%Y-%m-%d %H:%M:%S",
        )
    )
    root.addHandler(fh)

    _LOG = VerboseLoggerAdapter(root, {})
    _LOG.info(f"Logger initialized | debug={debug} verbose={verbose}")
    _LOG.debug(f"Debug log file at: {file_path}")
    return _LOG


def get_log() -> VerboseLoggerAdapter:
    global _LOG
    if _LOG is None:
        base = logging.getLogger("MHCXGraph")
        _LOG = VerboseLoggerAdapter(base, {})
    return _LOG
