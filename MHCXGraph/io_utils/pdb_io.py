import os
from pathlib import Path
from typing import List

def _match_ext(name: str, extensions: List[str]) -> bool:
    """Checagem simples por sufixo. Ex.: ['.pdb', '.pdb.gz', '.cif']"""
    return any(name.endswith(ext) for ext in extensions)

def list_pdb_files(pdb_dir: str, extensions: List[str] = None) -> List[str]:
    """
    Lista SOMENTE os arquivos do diretório (sem recursão) que batem com extensions
    e imprime o menu de seleção. Retorna apenas os NOMES dos arquivos (como antes).
    """
    extensions = extensions or [".pdb", ".cif"]
    dir_path = Path(pdb_dir).expanduser().resolve()

    if not dir_path.exists() or not dir_path.is_dir():
        raise Exception(f"Directory not found: {pdb_dir}")

    files = sorted(f for f in os.listdir(dir_path) if _match_ext(f, extensions))

    if not files:
        exts = ", ".join(extensions)
        raise Exception(f"No structure files ({exts}) found in: {pdb_dir}")

    print("Select structure files to associate:")
    print("1 - All")
    for i, fname in enumerate(files, start=2):
        print(f"{i} - {fname}")

    return files

def get_user_selection(pdb_files: List[str], pdb_dir: str):
    """
    Pergunta no terminal e retorna [[full_path, file_name], ...]
    Mantém a semântica antiga: '1' significa todos; caso contrário números separados por vírgula.
    """
    raw = input("\nEnter the numbers (comma-separated) or '1' for All: ").strip()

    if raw.lower() in {"1", "all", "*"}:
        return [[os.path.join(pdb_dir, f), f] for f in pdb_files]

    selected_files = []
    tokens = [t.strip() for t in raw.split(",") if t.strip()]

    for tok in tokens:
        if not tok.isdigit():
            print(f"Ignoring invalid token: {tok!r}")
            continue
        idx = int(tok)
        if idx == 1:
            # prioridade para All se aparecer misturado
            return [[os.path.join(pdb_dir, f), f] for f in pdb_files]
        mapped = idx - 2
        if 0 <= mapped < len(pdb_files):
            fname = pdb_files[mapped]
            selected_files.append([os.path.join(pdb_dir, fname), fname])
        else:
            print(f"Number {idx} is out of range. Skipping.")

    if not selected_files:
        print("No valid selections. Selecting none.")
    return selected_files


