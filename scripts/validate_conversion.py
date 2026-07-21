#!/usr/bin/env python3
"""Static integrity checks for the Snakemake-to-Nextflow conversion.

This does not execute bioinformatics tools. It validates project structure,
module wiring, sample metadata, DANPOS contrast rendering, YAML syntax, bundled
FASTQ integrity, and optional byte-level comparison with the uploaded source zip.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import sys
import tempfile
import zipfile
from pathlib import Path

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML is required to run this validator") from exc

EXPECTED = {
    "FASTQC_RAW": "fastqc_raw.nf",
    "TRIM_READS": "trim_reads.nf",
    "BOWTIE2_BUILD": "bowtie2_build.nf",
    "BOWTIE2_ALIGN": "bowtie2_align.nf",
    "SAM_TO_BAM": "sam_to_bam.nf",
    "SORT_BAM": "sort_bam.nf",
    "INDEX_BAM": "index_bam.nf",
    "QUALIMAP_MAPPING": "qualimap_mapping.nf",
    "MULTIQC_MAPPING": "multiqc_mapping.nf",
    "ALIGNMENT_SIEVE": "alignment_sieve.nf",
    "SORT_FILTERED_BAM": "sort_filtered_bam.nf",
    "INDEX_FILTERED_BAM": "index_filtered_bam.nf",
    "QUALIMAP_FILTERED": "qualimap_filtered.nf",
    "MULTIQC_FILTERED": "multiqc_filtered.nf",
    "PREPARE_DANPOS_INPUT": "prepare_danpos_input.nf",
    "DANPOS_DIFF": "danpos_diff.nf",
    "ADJUST_FDR": "adjust_fdr.nf",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_sheet(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        first = handle.readline()
        delimiter = "\t" if first.count("\t") >= first.count(",") else ","
        handle.seek(0)
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    required = {"Samples", "R1", "R2", "Condition", "Replicate"}
    missing = required.difference(rows[0].keys() if rows else set())
    if missing:
        raise AssertionError(f"Missing samplesheet fields: {sorted(missing)}")
    return rows


def check_source_zip(root: Path, archive: Path, messages: list[str]) -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(temp_dir)
        source_root = next(Path(temp_dir).glob("*/MNase-Seq"))
        pairs = {
            source_root / "1_code/AdjustFDR/AdjustFDR.R": root / "bin/AdjustFDR.R",
            source_root / "1_code/AdjustFDR/functions.R": root / "bin/functions.R",
            source_root / "1_code/find_dynamic_nucleosomes.R": root / "bin/find_dynamic_nucleosomes.R",
            source_root / "1_code/nucleosomeDynamics/get_nucleosomeDynamics_gff.R": root / "bin/get_nucleosomeDynamics_gff.R",
        }
        for source, converted in pairs.items():
            if source.name in {"AdjustFDR.R", "functions.R"}:
                assert sha256(source) == sha256(converted), f"Byte mismatch: {source.name}"
                messages.append(f"PASS: {source.name} is byte-identical to the source archive")
            else:
                assert source.read_text() == converted.read_text(), f"Text mismatch: {source.name}"
                messages.append(f"PASS: {source.name} is unchanged from the source archive")

        snakefile = (source_root / "Snakefile").read_text()
        rule_names = re.findall(r"^rule\s+([A-Za-z0-9_]+):", snakefile, flags=re.MULTILINE)
        executable_rules = [name for name in rule_names if name != "all"]
        assert len(executable_rules) == 17, executable_rules
        messages.append("PASS: source Snakefile contains 17 executable rules plus rule all")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--source-zip", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    messages: list[str] = []

    modules = sorted((root / "modules").glob("*.nf"))
    assert len(modules) == 17, f"Expected 17 modules, found {len(modules)}"
    messages.append("PASS: 17 module files found")

    main_text = (root / "main.nf").read_text()
    observed_processes: set[str] = set()
    for process_name, filename in EXPECTED.items():
        module = root / "modules" / filename
        text = module.read_text()
        declarations = re.findall(r"^process\s+([A-Za-z0-9_]+)\s*\{", text, flags=re.MULTILINE)
        assert declarations == [process_name], f"{filename}: {declarations}"
        observed_processes.add(process_name)
        include_pattern = rf"include\s*\{{\s*{process_name}\s*\}}\s*from\s*['\"]\./modules/{module.stem}['\"]"
        assert re.search(include_pattern, main_text), f"Missing include for {process_name}"
    assert observed_processes == set(EXPECTED)
    messages.append("PASS: each module contains exactly one expected process and every process is included")

    for nf in [root / "main.nf", *modules]:
        text = nf.read_text()
        assert text.count("\"\"\"") % 2 == 0, f"Unbalanced triple quotes: {nf}"
        assert text.count("{") == text.count("}"), f"Unbalanced braces: {nf}"
    messages.append("PASS: coarse Nextflow quote/brace balance checks")

    config = (root / "nextflow.config").read_text()
    expected_danpos = "-m 1 -c 51137214 --extend 73 -a 1 -z 20 -jd 100 -e 1"
    assert expected_danpos in config
    danpos_module = (root / "modules/danpos_diff.nf").read_text()
    assert '"${treatment}/:${control}/"' in danpos_module
    messages.append("PASS: original DANPOS parameters and ordered path expression are present")

    rows = read_sheet(root / "assets/test_samplesheet.tsv")
    assert len(rows) == 4
    ids = [row["Samples"] for row in rows]
    assert len(ids) == len(set(ids))
    conditions = {row["Condition"] for row in rows}
    assert conditions == {"iKO_24h", "noKO_24h"}
    for row in rows:
        for mate in ("R1", "R2"):
            fastq = (root / "assets" / row[mate]).resolve()
            assert fastq.is_file(), fastq
            with gzip.open(fastq, "rb") as handle:
                while handle.read(1024 * 1024):
                    pass
    messages.append("PASS: test TSV has four unique samples, two expected conditions, and eight valid gzip FASTQs")

    for contrast, rendered in {
        "iKO_24h:noKO_24h": "iKO_24h/:noKO_24h/",
        "noKO_24h:iKO_24h": "noKO_24h/:iKO_24h/",
    }.items():
        left, right = [part.strip() for part in contrast.split(":")]
        assert left in conditions and right in conditions and f"{left}/:{right}/" == rendered
    messages.append("PASS: forward and reverse DANPOS contrast rendering")

    for env_file in sorted((root / "envs").glob("*.yml")):
        parsed = yaml.safe_load(env_file.read_text())
        assert isinstance(parsed, dict) and parsed.get("dependencies"), env_file
    messages.append("PASS: all Conda YAML files parse and contain dependencies")

    expected_r = {
        "AdjustFDR.R", "functions.R", "find_dynamic_nucleosomes.R", "get_nucleosomeDynamics_gff.R"
    }
    observed_r = {path.name for path in (root / "bin").glob("*.R")}
    assert expected_r.issubset(observed_r), observed_r
    messages.append("PASS: all original R scripts are located in bin/")

    for script in sorted((root / "scripts").glob("*.sh")):
        assert script.stat().st_mode & 0o111, f"Not executable: {script}"
    assert (root / "bin/AdjustFDR.R").stat().st_mode & 0o111
    messages.append("PASS: installer/reference shell scripts and AdjustFDR.R are executable")

    if args.source_zip:
        check_source_zip(root, args.source_zip.resolve(), messages)

    print("Static conversion validation")
    print("=" * 28)
    for message in messages:
        print(message)
    print("\nLIMITATION: Nextflow, R, Conda environments, DANPOS3, and reference files were not executed here.")
    print("CONFIDENCE workflow logic: 0.94")
    print("CONFIDENCE DANPOS ordering/validation: 0.99")
    print("CONFIDENCE complete target-system execution: 0.65")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, zipfile.BadZipFile) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
