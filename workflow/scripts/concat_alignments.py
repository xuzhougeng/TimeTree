#!/usr/bin/env python3
"""
Concatenate multiple sequence alignments into a supermatrix.
Generate partition file in NEXUS format for IQ-TREE.
"""

import tempfile
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Set
from Bio import SeqIO


@dataclass
class SupermatrixResult:
    """Result of supermatrix concatenation."""
    supermatrix: Dict[str, str]  # species -> concatenated sequence
    partitions: List[Tuple[str, int, int]]  # (og_id, start, end)
    species: List[str]  # sorted species list
    total_length: int
    n_partitions: int


def parse_alignment(aln_file: Path) -> Dict[str, str]:
    """Parse FASTA alignment into dict of species -> sequence."""
    seqs = {}
    for record in SeqIO.parse(aln_file, "fasta"):
        seqs[record.id] = str(record.seq)
    return seqs


def load_alignments(aln_files: List[Path]) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """
    Load all alignment files.

    Args:
        aln_files: List of alignment file paths

    Returns:
        Tuple of (alignments dict, sorted species list)
        alignments: {og_id: {species: sequence}}
    """
    all_species: Set[str] = set()
    alignments: Dict[str, Dict[str, str]] = {}

    for aln_file in sorted(aln_files):
        og_id = aln_file.stem.replace('.aln', '')
        seqs = parse_alignment(aln_file)
        alignments[og_id] = seqs
        all_species.update(seqs.keys())

    return alignments, sorted(all_species)


def build_supermatrix(
    alignments: Dict[str, Dict[str, str]],
    species_list: List[str]
) -> SupermatrixResult:
    """
    Build supermatrix by concatenating alignments.

    Args:
        alignments: {og_id: {species: sequence}}
        species_list: Sorted list of all species

    Returns:
        SupermatrixResult with concatenated sequences and partition info
    """
    supermatrix: Dict[str, str] = defaultdict(str)
    partitions: List[Tuple[str, int, int]] = []
    pos = 1

    # Use consistent ordering
    sorted_og_ids = sorted(alignments.keys())

    for og_id in sorted_og_ids:
        seqs = alignments[og_id]
        # Get alignment length from first sequence
        aln_len = len(next(iter(seqs.values())))

        for sp in species_list:
            if sp in seqs:
                supermatrix[sp] += seqs[sp]
            else:
                # Missing species: fill with gaps
                supermatrix[sp] += '-' * aln_len

        # Record partition
        partitions.append((og_id, pos, pos + aln_len - 1))
        pos += aln_len

    total_len = len(next(iter(supermatrix.values()))) if supermatrix else 0

    # Validate partition integrity
    _validate_partitions(partitions, total_len, sorted_og_ids, alignments)

    return SupermatrixResult(
        supermatrix=dict(supermatrix),
        partitions=partitions,
        species=species_list,
        total_length=total_len,
        n_partitions=len(partitions)
    )


def _validate_partitions(
    partitions: List[Tuple[str, int, int]],
    total_len: int,
    sorted_og_ids: List[str],
    alignments: Dict[str, Dict[str, str]]
) -> None:
    """
    Validate partition integrity.

    Raises:
        ValueError: If partition validation fails
    """
    # Check 1: Partition count matches gene count
    if len(partitions) != len(sorted_og_ids):
        raise ValueError(
            f"Partition count ({len(partitions)}) != gene count ({len(sorted_og_ids)})"
        )

    # Check 2: Partitions are contiguous and cover entire length
    expected_pos = 1
    for i, (og_id, start, end) in enumerate(partitions):
        # Check start position
        if start != expected_pos:
            raise ValueError(
                f"Partition {og_id}: expected start {expected_pos}, got {start}"
            )

        # Check partition length matches alignment length
        partition_len = end - start + 1
        actual_len = len(next(iter(alignments[og_id].values())))
        if partition_len != actual_len:
            raise ValueError(
                f"Partition {og_id}: length mismatch. "
                f"Partition says {partition_len}, alignment has {actual_len}"
            )

        # Check order matches sorted order
        if og_id != sorted_og_ids[i]:
            raise ValueError(
                f"Partition order mismatch at index {i}: "
                f"expected {sorted_og_ids[i]}, got {og_id}"
            )

        expected_pos = end + 1

    # Check 3: Last partition ends at total length
    if partitions and partitions[-1][2] != total_len:
        raise ValueError(
            f"Last partition ends at {partitions[-1][2]}, "
            f"but total length is {total_len}"
        )

    # Check 4: Sum of partition lengths equals total length
    sum_partition_lens = sum(end - start + 1 for _, start, end in partitions)
    if sum_partition_lens != total_len:
        raise ValueError(
            f"Sum of partition lengths ({sum_partition_lens}) != "
            f"total length ({total_len})"
        )


def write_phylip(result: SupermatrixResult, out_path: Path) -> None:
    """Write supermatrix in PHYLIP format."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        f.write(f"{len(result.species)} {result.total_length}\n")
        for sp in result.species:
            f.write(f"{sp}  {result.supermatrix[sp]}\n")


def write_fasta(result: SupermatrixResult, out_path: Path) -> None:
    """Write supermatrix in FASTA format."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        for sp in result.species:
            f.write(f">{sp}\n{result.supermatrix[sp]}\n")


def write_partitions(result: SupermatrixResult, out_path: Path) -> None:
    """Write partition file in NEXUS format for IQ-TREE."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        f.write("#nexus\n")
        f.write("begin sets;\n")
        for og_id, start, end in result.partitions:
            f.write(f"  charset {og_id} = {start}-{end};\n")
        f.write("end;\n")


def concat_alignments(
    aln_files: List[Path],
    out_phy: Path,
    out_fasta: Path,
    out_part: Path
) -> SupermatrixResult:
    """
    Main concatenation logic.

    Args:
        aln_files: List of input alignment files
        out_phy: Output PHYLIP file path
        out_fasta: Output FASTA file path
        out_part: Output partition file path

    Returns:
        SupermatrixResult with concatenation details
    """
    # Load all alignments
    alignments, species_list = load_alignments(aln_files)

    # Build supermatrix
    result = build_supermatrix(alignments, species_list)

    # Write outputs
    write_phylip(result, out_phy)
    write_fasta(result, out_fasta)
    write_partitions(result, out_part)

    print(f"Supermatrix: {len(result.species)} taxa, {result.total_length} sites, "
          f"{result.n_partitions} partitions")

    return result


def main():
    """Snakemake entry point."""
    aln_files = [Path(f) for f in snakemake.input.alignments]
    out_phy = Path(snakemake.output.supermatrix)
    out_part = Path(snakemake.output.partitions)
    out_fasta = Path(snakemake.output.fasta)

    concat_alignments(aln_files, out_phy, out_fasta, out_part)


# =============================================================================
# Test Functions
# =============================================================================

def test_concat_alignments():
    """Test the concatenation logic with mock data."""
    import os

    # Create temporary directory for test files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create mock alignment files
        aln1 = tmpdir / "OG0001.aln.faa"
        aln2 = tmpdir / "OG0002.aln.faa"

        # Gene 1: all 3 species present
        aln1.write_text(
            ">Species_A\nMKVLFFVA\n"
            ">Species_B\nMKVLFFIA\n"
            ">Species_C\nMKILFFVA\n"
        )

        # Gene 2: Species_B missing
        aln2.write_text(
            ">Species_A\nACDEFG\n"
            ">Species_C\nACDEYG\n"
        )

        # Output paths
        out_phy = tmpdir / "output" / "supermatrix.phy"
        out_fasta = tmpdir / "output" / "supermatrix.faa"
        out_part = tmpdir / "output" / "partitions.nex"

        # Run concatenation
        result = concat_alignments(
            aln_files=[aln1, aln2],
            out_phy=out_phy,
            out_fasta=out_fasta,
            out_part=out_part
        )

        # Assertions
        print("\n=== Test Results ===")

        # Test 1: Check species count
        assert len(result.species) == 3, f"Expected 3 species, got {len(result.species)}"
        print("✓ Species count correct (3)")

        # Test 2: Check total length (8 + 6 = 14)
        assert result.total_length == 14, f"Expected length 14, got {result.total_length}"
        print("✓ Total length correct (14)")

        # Test 3: Check partition count
        assert result.n_partitions == 2, f"Expected 2 partitions, got {result.n_partitions}"
        print("✓ Partition count correct (2)")

        # Test 4: Check gap filling for missing species
        assert result.supermatrix["Species_B"].endswith("------"), \
            f"Expected gaps for missing Species_B in gene 2, got: {result.supermatrix['Species_B']}"
        print("✓ Gap filling correct for missing species")

        # Test 5: Check partition positions
        assert result.partitions[0] == ("OG0001", 1, 8), \
            f"Partition 1 incorrect: {result.partitions[0]}"
        assert result.partitions[1] == ("OG0002", 9, 14), \
            f"Partition 2 incorrect: {result.partitions[1]}"
        print("✓ Partition positions correct")

        # Test 6: Check output files exist
        assert out_phy.exists(), "PHYLIP file not created"
        assert out_fasta.exists(), "FASTA file not created"
        assert out_part.exists(), "Partition file not created"
        print("✓ All output files created")

        # Test 7: Verify PHYLIP format
        phy_content = out_phy.read_text()
        assert phy_content.startswith("3 14"), f"PHYLIP header incorrect: {phy_content[:20]}"
        print("✓ PHYLIP format correct")

        # Test 8: Verify partition file format
        part_content = out_part.read_text()
        assert "#nexus" in part_content.lower(), "Missing NEXUS header"
        assert "charset OG0001 = 1-8;" in part_content, "Partition 1 definition incorrect"
        assert "charset OG0002 = 9-14;" in part_content, "Partition 2 definition incorrect"
        print("✓ Partition file format correct")

        # Test 9: Verify sequence integrity
        expected_A = "MKVLFFVA" + "ACDEFG"  # 8 + 6 = 14
        expected_B = "MKVLFFIA" + "------"  # 8 + 6 gaps = 14
        expected_C = "MKILFFVA" + "ACDEYG"  # 8 + 6 = 14

        assert result.supermatrix["Species_A"] == expected_A, \
            f"Species_A sequence incorrect:\n  Expected: {expected_A}\n  Got: {result.supermatrix['Species_A']}"
        assert result.supermatrix["Species_B"] == expected_B, \
            f"Species_B sequence incorrect:\n  Expected: {expected_B}\n  Got: {result.supermatrix['Species_B']}"
        assert result.supermatrix["Species_C"] == expected_C, \
            f"Species_C sequence incorrect:\n  Expected: {expected_C}\n  Got: {result.supermatrix['Species_C']}"
        print("✓ Sequence integrity verified")

        print("\n=== All tests passed! ===\n")

        # Print example output
        print("Example supermatrix:")
        for sp, seq in result.supermatrix.items():
            print(f"  {sp}: {seq}")

        return True


def test_edge_cases():
    """Test edge cases."""
    import os

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        print("\n=== Edge Case Tests ===")

        # Edge case 1: Single gene
        aln1 = tmpdir / "OG0001.aln.faa"
        aln1.write_text(">Sp_A\nMKVL\n>Sp_B\nMKIL\n")

        out_phy = tmpdir / "single" / "super.phy"
        out_fasta = tmpdir / "single" / "super.faa"
        out_part = tmpdir / "single" / "parts.nex"

        result = concat_alignments([aln1], out_phy, out_fasta, out_part)
        assert result.n_partitions == 1, "Single gene should have 1 partition"
        assert result.total_length == 4, "Single gene length incorrect"
        print("✓ Single gene case handled correctly")

        # Edge case 2: Species with special characters in name
        aln2 = tmpdir / "OG0002.aln.faa"
        aln2.write_text(">Species_name_with_underscores\nACDE\n>Another-species.v1\nACDF\n")

        out_phy2 = tmpdir / "special" / "super.phy"
        out_fasta2 = tmpdir / "special" / "super.faa"
        out_part2 = tmpdir / "special" / "parts.nex"

        result2 = concat_alignments([aln2], out_phy2, out_fasta2, out_part2)
        assert "Species_name_with_underscores" in result2.species
        assert "Another-species.v1" in result2.species
        print("✓ Special characters in species names handled correctly")

        print("\n=== Edge case tests passed! ===\n")
        return True


def test_partition_validation():
    """Test that partition validation catches errors."""
    print("\n=== Partition Validation Tests ===")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create multiple genes with known lengths
        genes = [
            ("OG0001", "MKVL", 4),      # length 4
            ("OG0002", "ACDEFG", 6),    # length 6
            ("OG0003", "HIJ", 3),       # length 3
        ]

        aln_files = []
        for og_id, seq, _ in genes:
            aln_file = tmpdir / f"{og_id}.aln.faa"
            aln_file.write_text(f">Species_A\n{seq}\n>Species_B\n{seq}\n")
            aln_files.append(aln_file)

        out_phy = tmpdir / "out" / "super.phy"
        out_fasta = tmpdir / "out" / "super.faa"
        out_part = tmpdir / "out" / "parts.nex"

        result = concat_alignments(aln_files, out_phy, out_fasta, out_part)

        # Verify partition positions
        # OG0001: 1-4, OG0002: 5-10, OG0003: 11-13
        expected_partitions = [
            ("OG0001", 1, 4),
            ("OG0002", 5, 10),
            ("OG0003", 11, 13),
        ]

        for i, (expected, actual) in enumerate(zip(expected_partitions, result.partitions)):
            assert expected == actual, \
                f"Partition {i} mismatch: expected {expected}, got {actual}"

        # Verify total length
        assert result.total_length == 13, f"Expected total length 13, got {result.total_length}"

        # Verify sequences are correctly concatenated
        assert len(result.supermatrix["Species_A"]) == 13
        assert result.supermatrix["Species_A"] == "MKVLACDEFGHIJ"

        print("✓ Partition positions verified: OG0001(1-4), OG0002(5-10), OG0003(11-13)")
        print("✓ Total length verified: 13")
        print("✓ Sequence concatenation verified")

        # Read partition file and verify content
        part_content = out_part.read_text()
        assert "charset OG0001 = 1-4;" in part_content
        assert "charset OG0002 = 5-10;" in part_content
        assert "charset OG0003 = 11-13;" in part_content
        print("✓ Partition file content verified")

        print("\n=== Partition validation tests passed! ===\n")
        return True


if __name__ == "__main__":
    # Check if running under Snakemake
    try:
        snakemake
        main()
    except NameError:
        # Running standalone - execute tests
        print("Running in standalone mode - executing tests\n")
        test_concat_alignments()
        test_edge_cases()
        test_partition_validation()
