# core/tracking.py
from __future__ import annotations
import os, pickle, json, time, random, string, threading
from pathlib import Path
from dataclasses import dataclass
from typing import Any, Optional, Callable, Union
import numpy as np

__all__ = [
    "init_tracker", "get_current", "set_enabled",
    "save", "mark", "track", "tracked", "Tracker"
]

# ---------- Config ----------

@dataclass
class TrackerConfig:
    root: Union[str, Path] = "CrossSteps"     # diretório raiz p/ todas as execuções
    outdir: Optional[str] = None              # subdir da execução (ex: run_name). Se None, gera aleatório.
    enabled: bool = True                      # liga/desliga persistência
    prefer_npy_for_ndarray: bool = False      # se True, arrays vão como .npy ao invés de .pkl
    add_timestamp_prefix: bool = False        # se True, inclui epoch no prefixo (além do contador)

# ---------- Núcleo ----------

class Tracker:
    def __init__(self, cfg: TrackerConfig):
        self.cfg = cfg
        self._lock = threading.Lock()
        self._step = 0

        run_dir = cfg.outdir or _rand_id(6)
        self.base = Path(cfg.root) / run_dir
        self.base.mkdir(parents=True, exist_ok=True)

        # metadados mínimos da execução
        meta = {
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "run_dir": str(self.base),
            "enabled": cfg.enabled,
        }
        (self.base / "_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    # contador incremental seguro (thread-safe)
    def _next_index(self) -> int:
        with self._lock:
            self._step += 1
            return self._step

    def _prefix(self, step_label: str) -> str:
        idx = self._next_index()
        ts = f"{int(time.time())}_" if self.cfg.add_timestamp_prefix else ""
        # prefixo: 001_[ts]stepLabel_
        return f"{idx:03d}_{ts}{step_label}_" if step_label else f"{idx:03d}_"

    def _ensure_dir(self, rel: Optional[Union[str, Path]] = None) -> Path:
        d = self.base if rel is None else (self.base / rel)
        d.mkdir(parents=True, exist_ok=True)
        return d

    def save(
        self,
        step_label: str,
        key: str,
        obj: Any,
        *,
        subdir: Optional[str] = None,
        filename: Optional[str] = None
    ) -> Optional[Path]:
        """Salva qualquer objeto. Nome: 'NNN_stepLabel_key.pkl' (ou .npy p/ arrays)."""
        if not self.cfg.enabled:
            return None

        folder = self._ensure_dir(subdir)
        prefix = self._prefix(step_label)
        if filename is None:
            stem = f"{prefix}{_sanitize(key)}"
        else:
            stem = _sanitize(filename)

        # arrays podem ir em .npy se preferir
        if isinstance(obj, np.ndarray) and self.cfg.prefer_npy_for_ndarray:
            fpath = folder / f"{stem}.npy"
            np.save(fpath, obj)
            return fpath

        # fallback geral: pickle
        fpath = folder / f"{stem}.pkl"
        with open(fpath, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
        return fpath

    def mark(
        self,
        step_label: str,
        key: str,
        info: Union[str, dict, list, int, float, bool, None],
        *,
        subdir: Optional[str] = None
    ) -> Optional[Path]:
        """Marca um checkpoint leve em .txt (se str) ou .json (se não for str)."""
        if not self.cfg.enabled:
            return None

        folder = self._ensure_dir(subdir)
        prefix = self._prefix(step_label)
        stem = f"{prefix}{_sanitize(key)}"

        if isinstance(info, str):
            p = folder / f"{stem}.txt"
            p.write_text(info, encoding="utf-8")
            return p
        else:
            p = folder / f"{stem}.json"
            p.write_text(json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
            return p

    def track(
        self,
        step_label: str,
        key: str,
        func: Callable[..., Any],
        *args, subdir: Optional[str] = None, **kwargs
    ) -> Any:
        """Executa a função, salva o resultado e retorna o valor."""
        result = func(*args, **kwargs)
        self.save(step_label, key, result, subdir=subdir)
        return result

    # Context manager opcional (para trocar temporariamente o tracker global)
    def __enter__(self):
        self._prev = get_current(_default_none=True)
        _set_global(self)
        return self

    def __exit__(self, exc_type, exc, tb):
        _set_global(self._prev)

# ---------- Singleton global ----------

_global_tracker: Optional[Tracker] = None

def _rand_id(n: int) -> str:
    return "".join(random.choices(string.ascii_lowercase + string.digits, k=n))

def _sanitize(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in str(s))

def _set_global(t: Optional[Tracker]) -> None:
    global _global_tracker
    _global_tracker = t

def init_tracker(
    *,
    root: Union[str, Path] = "CrossSteps",
    outdir: Optional[str] = None,
    enabled: bool = True,
    prefer_npy_for_ndarray: bool = False,
    add_timestamp_prefix: bool = False
) -> Tracker:
    """Inicializa/atualiza o tracker global com sua configuração."""
    cfg = TrackerConfig(
        root=root,
        outdir=outdir,
        enabled=enabled,
        prefer_npy_for_ndarray=prefer_npy_for_ndarray,
        add_timestamp_prefix=add_timestamp_prefix,
    )
    tracker = Tracker(cfg)
    _set_global(tracker)
    return tracker

def get_current(_default_none: bool = False) -> Optional[Tracker]:
    if _global_tracker is None:
        if _default_none:
            return None  # usado apenas internamente
        # default “preguiçoso” para não quebrar quem esquecer de inicializar
        return init_tracker()  # CrossSteps/<rand>
    return _global_tracker

def set_enabled(flag: bool) -> None:
    t = get_current()
    t.cfg.enabled = flag

# ---------- atalhos globais (back‑compat) ----------

def save(step: str, key: str, obj: Any, **opts) -> Optional[Path]:
    """Compatível com seu uso: save('create_graphs','Subgraph', obj)."""
    return get_current().save(step, key, obj, **opts)

def mark(step: str, key: str, info: Any, **opts) -> Optional[Path]:
    return get_current().mark(step, key, info, **opts)

def track(step: str, key: str, func: Callable[..., Any], *args, **kwargs) -> Any:
    return get_current().track(step, key, func, *args, **kwargs)

# Decorator opcional para funções – salva retorno automaticamente
def tracked(step: str, key: str, *, subdir: Optional[str] = None):
    def deco(fn: Callable[..., Any]):
        def wrapper(*args, **kwargs):
            res = fn(*args, **kwargs)
            get_current().save(step, key, res, subdir=subdir)
            return res
        return wrapper
    return deco
