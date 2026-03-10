import argparse
import copy
import os
from pathlib import Path

import pandas as pd
from Bio import pairwise2
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa

# =========================
# Basic amino-acid mapping
# =========================

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
    "MSE": "M",
}


def residue_to_one_letter(residue) -> str:
    resname = residue.get_resname().strip().upper()
    return THREE_TO_ONE.get(resname, "X")


def parse_imgt_label(label: str) -> tuple[int, str]:
    label = str(label).strip()
    if not label:
        raise ValueError("Empty IMGT label")

    i = 0
    while i < len(label) and label[i].isdigit():
        i += 1

    if i == 0:
        raise ValueError(f"Invalid IMGT label: {label}")

    resseq = int(label[:i])
    icode = label[i:] if i < len(label) else " "
    if len(icode) > 1:
        raise ValueError(f"Unsupported insertion code in label: {label}")

    return resseq, icode


# =========================
# Structure I/O
# =========================

def load_structure(path: str, structure_id: str = "structure"):
    ext = os.path.splitext(path)[1].lower()
    if ext in [".cif", ".mmcif"]:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(structure_id, path)


def save_structure(structure, output_path: str):
    ext = os.path.splitext(output_path)[1].lower()
    if ext in [".cif", ".mmcif"]:
        io = MMCIFIO()
    else:
        io = PDBIO()
    io.set_structure(structure)
    io.save(output_path)


# =========================
# Sequence extraction
# =========================

def get_chain_polymer_residues(chain) -> list:
    residues = []
    for residue in chain:
        if is_aa(residue, standard=False):
            residues.append(residue)
    return residues


def get_chain_sequence_and_residues(chain) -> tuple[str, list]:
    residues = get_chain_polymer_residues(chain)
    seq = "".join(residue_to_one_letter(r) for r in residues)
    return seq, residues


# =========================
# Template loading
# =========================

def load_mhcii_templates(display_csv: Path, numbering_csv: Path):
    display_df = pd.read_csv(display_csv)
    numbering_df = pd.read_csv(numbering_csv)

    required_display_cols = {"seq", "type"}
    required_number_cols = {"MHCII_a", "MHCII_b"}

    if not required_display_cols.issubset(display_df.columns):
        raise ValueError(
            f"imgt_display_all.csv must contain columns {required_display_cols}, "
            f"but has {display_df.columns.tolist()}"
        )

    if not required_number_cols.issubset(numbering_df.columns):
        raise ValueError(
            f"imgt_numbering_mapping_all.csv must contain columns {required_number_cols}, "
            f"but has {numbering_df.columns.tolist()}"
        )

    templates = {"MHCII_a": [], "MHCII_b": []}
    for mhc_type in ["MHCII_a", "MHCII_b"]:
        sub = display_df[display_df["type"] == mhc_type].copy()
        numbering = numbering_df[mhc_type].astype(str).tolist()

        for _, row in sub.iterrows():
            gapped_seq = str(row["seq"]).strip()
            if len(gapped_seq) != len(numbering):
                raise ValueError(
                    f"Template length != numbering length for {mhc_type}. "
                    f"Template length={len(gapped_seq)}, numbering length={len(numbering)}"
                )

            ungapped_seq = gapped_seq.replace("-", "")
            template_name = f"{row.get('class', 'NA')}|{row.get('allele', 'NA')}|{mhc_type}"

            templates[mhc_type].append({
                "name": template_name,
                "type": mhc_type,
                "gapped_seq": gapped_seq,
                "ungapped_seq": ungapped_seq,
                "numbering": numbering,
            })

    return templates


# =========================
# Alignment and scoring
# =========================

def align_query_to_template(query_seq: str, template_ungapped: str):
    alns = pairwise2.align.globalms(
        template_ungapped,
        query_seq,
        2.0,
        -1.0,
        -10.0,
        -1.0,
        penalize_end_gaps=(False, False),
        one_alignment_only=True
    )
    if not alns:
        return None
    return alns[0]


def compute_alignment_stats(aln) -> dict[str, float]:
    a = aln.seqA
    b = aln.seqB

    matches = 0
    aligned_pairs = 0
    template_non_gap = 0
    query_non_gap = 0

    for x, y in zip(a, b):
        if x != "-":
            template_non_gap += 1
        if y != "-":
            query_non_gap += 1
        if x != "-" and y != "-":
            aligned_pairs += 1
            if x == y:
                matches += 1

    identity_over_aligned = matches / aligned_pairs if aligned_pairs > 0 else 0.0
    coverage_template = aligned_pairs / template_non_gap if template_non_gap > 0 else 0.0
    coverage_query = aligned_pairs / query_non_gap if query_non_gap > 0 else 0.0

    return {
        "score": aln.score,
        "matches": matches,
        "aligned_pairs": aligned_pairs,
        "template_len": template_non_gap,
        "query_len": query_non_gap,
        "identity_over_aligned": identity_over_aligned,
        "coverage_template": coverage_template,
        "coverage_query": coverage_query,
    }


def choose_best_template_for_chain(query_seq: str, templates_for_type: list[dict]) -> dict:
    best = None

    for tpl in templates_for_type:
        aln = align_query_to_template(query_seq, tpl["ungapped_seq"])
        if aln is None:
            continue

        stats = compute_alignment_stats(aln)

        min_cov = stats["coverage_template"]
        if min_cov < 0.50:
            continue

        rank_tuple = (
            min_cov,
            stats["identity_over_aligned"],
            stats["score"]
        )

        if best is None or rank_tuple > best["rank_tuple"]:
            best = {
                "template": tpl,
                "alignment": aln,
                "stats": stats,
                "rank_tuple": rank_tuple,
            }

    return best


# =========================
# Projection
# =========================

def project_query_onto_gapped_template(
    gapped_template: str,
    numbering: list[str],
    aligned_template_ungapped: str,
    aligned_query: str,
    chain_residues: list
):
    if len(gapped_template) != len(numbering):
        raise ValueError("gapped_template and numbering must have same length")

    template_pos_to_query_pos = {}
    insertions_before_template_pos = {}

    tpos = 0
    qpos = 0

    insertions_before_template_pos[0] = []

    for t_char, q_char in zip(aligned_template_ungapped, aligned_query):
        if t_char != "-" and q_char != "-":
            template_pos_to_query_pos[tpos] = qpos
            tpos += 1
            qpos += 1
            if tpos not in insertions_before_template_pos:
                insertions_before_template_pos[tpos] = []

        elif t_char != "-" and q_char == "-":
            template_pos_to_query_pos[tpos] = None
            tpos += 1
            if tpos not in insertions_before_template_pos:
                insertions_before_template_pos[tpos] = []

        elif t_char == "-" and q_char != "-":
            if tpos not in insertions_before_template_pos:
                insertions_before_template_pos[tpos] = []
            insertions_before_template_pos[tpos].append(qpos)
            qpos += 1

    n_template = tpos
    if n_template not in insertions_before_template_pos:
        insertions_before_template_pos[n_template] = []

    projected_query_chars = []
    projected_residues = []
    numbering_map = []
    still_unplaced_query_positions = []

    gapped_pos = 0
    current_template_pos = 0

    while gapped_pos < len(gapped_template):
        tpl_char = gapped_template[gapped_pos]

        if tpl_char == "-":
            start = gapped_pos
            while gapped_pos < len(gapped_template) and gapped_template[gapped_pos] == "-":
                gapped_pos += 1
            end = gapped_pos

            insertion_bucket = insertions_before_template_pos.get(current_template_pos, [])

            gap_block_len = end - start
            n_to_place = min(gap_block_len, len(insertion_bucket))

            for i in range(gap_block_len):
                col_num = numbering[start + i]

                if i < n_to_place:
                    q_idx = insertion_bucket[i]
                    aa = residue_to_one_letter(chain_residues[q_idx])
                    projected_query_chars.append(aa)
                    projected_residues.append(chain_residues[q_idx])
                    numbering_map.append(col_num)
                else:
                    projected_query_chars.append("-")
                    projected_residues.append(None)
                    numbering_map.append(col_num)

            if len(insertion_bucket) > gap_block_len:
                still_unplaced_query_positions.extend(insertion_bucket[gap_block_len:])

            insertions_before_template_pos[current_template_pos] = []

        else:
            num_label = numbering[gapped_pos]
            q_idx = template_pos_to_query_pos.get(current_template_pos, None)

            if q_idx is None:
                projected_query_chars.append("-")
                projected_residues.append(None)
                numbering_map.append(num_label)
            else:
                aa = residue_to_one_letter(chain_residues[q_idx])
                projected_query_chars.append(aa)
                projected_residues.append(chain_residues[q_idx])
                numbering_map.append(num_label)

            current_template_pos += 1
            gapped_pos += 1

    for bucket in sorted(insertions_before_template_pos):
        if insertions_before_template_pos[bucket]:
            still_unplaced_query_positions.extend(insertions_before_template_pos[bucket])

    return projected_query_chars, projected_residues, numbering_map, still_unplaced_query_positions


# =========================
# Debug printing
# =========================

def print_alignment_blocks(label: str, seq: str, width: int = 80):
    for i in range(0, len(seq), width):
        print(f"{label:<12} {seq[i:i+width]}")


def print_debug_alignment(chain_id: str, mhc_type: str, template_name: str,
                          aln, gapped_template: str,
                          projected_query_chars: list[str], numbering: list[str]):
    print("\n" + "=" * 120)
    print(f"DEBUG | chain={chain_id} | assigned_type={mhc_type} | template={template_name}")
    print("-" * 120)
    print("Raw alignment against UNGAPPED template:")
    print_alignment_blocks("TEMPLATE", aln.seqA)
    print_alignment_blocks("QUERY", aln.seqB)

    match_line = "".join("|" if (a == b and a != "-" and b != "-") else " " for a, b in zip(aln.seqA, aln.seqB))
    print_alignment_blocks("MATCH", match_line)

    print("-" * 120)
    print("Projection onto FIXED GAPPED template:")
    projected_query = "".join(projected_query_chars)
    print_alignment_blocks("TEMPLATE", gapped_template)
    print_alignment_blocks("QUERY", projected_query)

    print("-" * 120)
    print("Projected numbering map (only positions where query residue exists):")
    pairs = []
    for aa, num in zip(projected_query_chars, numbering):
        if aa != "-" and num is not None:
            pairs.append(f"{num}:{aa}")
    chunk = 18
    for i in range(0, len(pairs), chunk):
        print("MAP         " + "  ".join(pairs[i:i+chunk]))
    print("=" * 120 + "\n")


# =========================
# Chain detection
# =========================

def detect_mhcII_chains(model, templates, debug: bool = False):
    chain_candidates = []

    for chain in model:
        seq, residues = get_chain_sequence_and_residues(chain)
        if len(seq) < 50:
            continue

        best_a = choose_best_template_for_chain(seq, templates["MHCII_a"])
        best_b = choose_best_template_for_chain(seq, templates["MHCII_b"])

        if best_a is None and best_b is None:
            continue

        chain_candidates.append({
            "chain": chain,
            "chain_id": chain.id,
            "seq": seq,
            "residues": residues,
            "best_a": best_a,
            "best_b": best_b,
        })

    if not chain_candidates:
        raise RuntimeError("Could not find any suitable protein chain candidates for MHC-II.")

    alpha_sorted = sorted(
        [c for c in chain_candidates if c["best_a"] is not None],
        key=lambda x: x["best_a"]["rank_tuple"],
        reverse=True
    )
    beta_sorted = sorted(
        [c for c in chain_candidates if c["best_b"] is not None],
        key=lambda x: x["best_b"]["rank_tuple"],
        reverse=True
    )

    if not alpha_sorted:
        raise RuntimeError("Could not detect an MHCII_a chain.")
    if not beta_sorted:
        raise RuntimeError("Could not detect an MHCII_b chain.")

    alpha_choice = alpha_sorted[0]
    beta_choice = None
    for c in beta_sorted:
        if c["chain_id"] != alpha_choice["chain_id"]:
            beta_choice = c
            break

    if beta_choice is None:
        raise RuntimeError("Could not detect a distinct MHCII_b chain.")

    result = {
        "MHCII_a": {
            "chain": alpha_choice["chain"],
            "chain_id": alpha_choice["chain_id"],
            "seq": alpha_choice["seq"],
            "residues": alpha_choice["residues"],
            "best": alpha_choice["best_a"],
        },
        "MHCII_b": {
            "chain": beta_choice["chain"],
            "chain_id": beta_choice["chain_id"],
            "seq": beta_choice["seq"],
            "residues": beta_choice["residues"],
            "best": beta_choice["best_b"],
        }
    }

    if debug:
        print("\nCandidate chain summary:")
        for c in chain_candidates:
            a_txt = "NA"
            b_txt = "NA"
            if c["best_a"] is not None:
                a = c["best_a"]["stats"]
                a_txt = (
                    f"score={a['score']:.1f}, "
                    f"id={a['identity_over_aligned']:.3f}, "
                    f"covQ={a['coverage_query']:.3f}, "
                    f"covT={a['coverage_template']:.3f}, "
                    f"tpl={c['best_a']['template']['name']}"
                )
            if c["best_b"] is not None:
                b = c["best_b"]["stats"]
                b_txt = (
                    f"score={b['score']:.1f}, "
                    f"id={b['identity_over_aligned']:.3f}, "
                    f"covQ={b['coverage_query']:.3f}, "
                    f"covT={b['coverage_template']:.3f}, "
                    f"tpl={c['best_b']['template']['name']}"
                )

            print(f"  Chain {c['chain_id']}: len={len(c['seq'])}")
            print(f"    as MHCII_a -> {a_txt}")
            print(f"    as MHCII_b -> {b_txt}")

        print("\nSelected chains:")
        print(f"  MHCII_a -> chain {result['MHCII_a']['chain_id']} using {result['MHCII_a']['best']['template']['name']}")
        print(f"  MHCII_b -> chain {result['MHCII_b']['chain_id']} using {result['MHCII_b']['best']['template']['name']}")

    return result


# =========================
# Renumbering
# =========================

def renumber_chain_from_projection(chain, projected_residues, numbering_map, debug=False):
    mapping = []
    seen_new_ids = set()

    for residue, num_label in zip(projected_residues, numbering_map):
        if residue is None or num_label is None:
            continue

        resseq, icode = parse_imgt_label(num_label)
        new_id = (" ", resseq, icode)

        if new_id in seen_new_ids:
            raise RuntimeError(
                f"Duplicate target numbering {new_id} inside chain {chain.id}. "
                f"This usually means the projection produced duplicated labels."
            )

        old_id = residue.id
        mapping.append((old_id, new_id, residue.get_resname()))
        seen_new_ids.add(new_id)

    for old_id, new_id, _ in mapping:
        residue = chain[old_id]
        residue.id = new_id

    if debug:
        print(f"\nRenumbering chain {chain.id}:")
        for old_id, new_id, resname in mapping:
            print(f"  {resname:>3}  {old_id}  ->  {new_id}")

    return mapping


def continue_numbering_from_unmatched(
    aa_residues_in_order,
    numbering_map,
    projected_residues,
    unmatched_query_positions,
    debug=False
):
    if not unmatched_query_positions:
        return

    assigned_numbers = []
    for residue, num_label in zip(projected_residues, numbering_map):
        if residue is None or num_label is None:
            continue
        resseq, _icode = parse_imgt_label(num_label)
        assigned_numbers.append(resseq)

    if not assigned_numbers:
        return

    next_resseq = max(assigned_numbers) + 1

    for q_idx in sorted(unmatched_query_positions):
        residue = aa_residues_in_order[q_idx]
        old_id = residue.id
        residue.id = (" ", next_resseq, " ")
        if debug:
            print(f"  CONT  {residue.get_resname():>3} {old_id} -> {residue.id}")
        next_resseq += 1


# =========================
# Per-file processing
# =========================

def process_structure_file_mhcii(
    input_path: str,
    output_path: str,
    templates,
    alpha_chain: str = None,
    beta_chain: str = None,
    debug: bool = False,
    warn_score: float = 50.0
):
    structure = load_structure(input_path, structure_id=os.path.basename(input_path))
    model = next(structure.get_models())

    detected = detect_mhcII_chains(model, templates, debug=debug)

    if alpha_chain is not None:
        if alpha_chain not in model:
            raise RuntimeError(f"alpha-chain '{alpha_chain}' not found in structure")
        seq, residues = get_chain_sequence_and_residues(model[alpha_chain])
        best = choose_best_template_for_chain(seq, templates["MHCII_a"])
        if best is None:
            raise RuntimeError(f"Could not assign chain {alpha_chain} as MHCII_a")
        detected["MHCII_a"] = {
            "chain": model[alpha_chain],
            "chain_id": alpha_chain,
            "seq": seq,
            "residues": residues,
            "best": best,
        }

    if beta_chain is not None:
        if beta_chain not in model:
            raise RuntimeError(f"beta-chain '{beta_chain}' not found in structure")
        seq, residues = get_chain_sequence_and_residues(model[beta_chain])
        best = choose_best_template_for_chain(seq, templates["MHCII_b"])
        if best is None:
            raise RuntimeError(f"Could not assign chain {beta_chain} as MHCII_b")
        detected["MHCII_b"] = {
            "chain": model[beta_chain],
            "chain_id": beta_chain,
            "seq": seq,
            "residues": residues,
            "best": best,
        }

    if detected["MHCII_a"]["chain_id"] == detected["MHCII_b"]["chain_id"]:
        raise RuntimeError("MHCII_a and MHCII_b resolved to the same chain.")

    out_structure = copy.deepcopy(structure)
    out_model = next(out_structure.get_models())

    chain_info = {}
    for mhc_type in ["MHCII_a", "MHCII_b"]:
        original = detected[mhc_type]
        out_chain = out_model[original["chain_id"]]
        _out_seq, out_residues = get_chain_sequence_and_residues(out_chain)

        best = original["best"]
        tpl = best["template"]
        aln = best["alignment"]

        score = best["stats"]["score"]
        if score < warn_score:
            print(
                f"WARNING: low alignment score in {os.path.basename(input_path)} | "
                f"{mhc_type} chain {original['chain_id']} | "
                f"template {tpl['name']} | score={score:.1f}"
            )

        projected_query_chars, projected_residues, numbering_map, unmatched_query_positions = (
            project_query_onto_gapped_template(
                tpl["gapped_seq"],
                tpl["numbering"],
                aln.seqA,
                aln.seqB,
                out_residues
            )
        )

        if debug and unmatched_query_positions:
            print(
                f"WARNING: chain {original['chain_id']} has "
                f"{len(unmatched_query_positions)} query residues that could not be fitted into "
                f"available numbered template columns."
            )
            print("Unplaced query residue indices:", unmatched_query_positions)

        if debug:
            print_debug_alignment(
                chain_id=original["chain_id"],
                mhc_type=mhc_type,
                template_name=tpl["name"],
                aln=aln,
                gapped_template=tpl["gapped_seq"],
                projected_query_chars=projected_query_chars,
                numbering=tpl["numbering"]
            )

        chain_info[mhc_type] = {
            "chain": out_chain,
            "aa_residues_in_order": out_residues,
            "projected_residues": projected_residues,
            "numbering_map": numbering_map,
            "unmatched_query_positions": unmatched_query_positions
        }

    for mhc_type in ["MHCII_a", "MHCII_b"]:
        renumber_chain_from_projection(
            chain_info[mhc_type]["chain"],
            chain_info[mhc_type]["projected_residues"],
            chain_info[mhc_type]["numbering_map"],
            debug=debug
        )

    for mhc_type in ["MHCII_a", "MHCII_b"]:
        continue_numbering_from_unmatched(
            chain_info[mhc_type]["aa_residues_in_order"],
            chain_info[mhc_type]["numbering_map"],
            chain_info[mhc_type]["projected_residues"],
            chain_info[mhc_type]["unmatched_query_positions"],
            debug=debug
        )

    save_structure(out_structure, output_path)


# =========================
# Main
# =========================

def main():
    parser = argparse.ArgumentParser(
        description="Renumber all MHC-II structures in a directory using IMGT templates and numbering."
    )
    parser.add_argument("-i", "--input-dir", required=True, help="Input directory with .pdb/.cif/.mmcif files")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--numbering-csv", required=True, help="Path to imgt_numbering_mapping_all.csv")
    parser.add_argument("--display-csv", required=True, help="Path to imgt_display_all.csv")
    parser.add_argument("--alpha-chain", default=None, help="Optional fixed chain ID for MHCII_a")
    parser.add_argument("--beta-chain", default=None, help="Optional fixed chain ID for MHCII_b")
    parser.add_argument("--warn-score",type=float,default=50.0,help="Warn if the best alignment score is below this value")
    parser.add_argument("--debug", action="store_true", help="Print debug output")
    parser.add_argument("--suffix", default="", help="Optional suffix added before file extension in output")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    templates = load_mhcii_templates(args.display_csv, args.numbering_csv)

    valid_ext = {".pdb", ".cif", ".mmcif"}
    files = sorted(
        f for f in os.listdir(args.input_dir)
        if os.path.isfile(os.path.join(args.input_dir, f))
        and os.path.splitext(f)[1].lower() in valid_ext
    )

    if not files:
        raise RuntimeError("No .pdb, .cif, or .mmcif files found in input directory.")

    n_ok = 0
    n_fail = 0

    for fname in files:
        input_path = os.path.join(args.input_dir, fname)
        stem, ext = os.path.splitext(fname)
        out_name = f"{stem}{args.suffix}{ext}" if args.suffix else fname
        output_path = os.path.join(args.output_dir, out_name)

        print(f"\nProcessing: {fname}")
        try:
            process_structure_file_mhcii(
                input_path=input_path,
                output_path=output_path,
                templates=templates,
                alpha_chain=args.alpha_chain,
                beta_chain=args.beta_chain,
                debug=args.debug,
                warn_score=args.warn_score
            )
            print(f"  OK -> {output_path}")
            n_ok += 1
        except Exception as e:
            print(f"  FAILED -> {fname}: {e}")
            n_fail += 1

    print("\nFinished.")
    print(f"  Success: {n_ok}")
    print(f"  Failed : {n_fail}")


if __name__ == "__main__":
    main()
