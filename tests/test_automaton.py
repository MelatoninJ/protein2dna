"""Unit tests for the Aho-Corasick automaton on its own (no codons, no host)."""

import random

import pytest

from oracle import hits
from p2d.optimize.automaton import MAX_EXPANSION, Automaton, expand_pattern, reverse_complement


def brute_force(text: str, patterns: list[str]) -> set[tuple[int, int]]:
    """(end, length) of the longest forbidden match ending at each position, both strands."""
    best: dict[int, int] = {}
    for start, end in hits(text, patterns):
        best[end] = max(best.get(end, 0), end - start)
    return set(best.items())


def test_reverse_complement():
    assert reverse_complement("GGTCTC") == "GAGACC"
    assert reverse_complement("CATATG") == "CATATG"  # palindromic


def test_expand_plain_and_ambiguous():
    assert expand_pattern("GAATTC") == ["GAATTC"]
    assert sorted(expand_pattern("GRCT")) == ["GACT", "GGCT"]
    assert len(expand_pattern("NN")) == 16


def test_expand_rejects_blowup_and_garbage():
    with pytest.raises(ValueError, match="expands to more than"):
        expand_pattern("N" * 7)
    assert len(expand_pattern("N" * 6)) == 4096 == MAX_EXPANSION
    with pytest.raises(ValueError, match="invalid IUPAC"):
        expand_pattern("GAXTC")
    with pytest.raises(ValueError, match="empty"):
        expand_pattern("")


def test_finds_pattern_on_both_strands():
    a = Automaton(["GGTCTC"])  # BsaI: not a palindrome
    assert [(m.start, m.length) for m in a.find_all("AAGGTCTCAA")] == [(2, 6)]
    assert [(m.start, m.length) for m in a.find_all("AAGAGACCAA")] == [(2, 6)]
    assert a.find_all("AAGGTCTAA") == []


def test_reports_original_pattern_not_expansion():
    a = Automaton(["GRATC"])
    (m,) = a.find_all("TTGAATCTT")
    assert m.pattern == "GRATC"


def test_overlapping_matches_are_all_seen():
    a = Automaton(["AAA"])
    assert [m.end for m in a.find_all("AAAAA")] == [3, 4, 5]


def test_longest_match_wins_when_patterns_nest():
    """A pattern that is a suffix of another must not hide the longer one."""
    a = Automaton(["CATATG", "ATG"])
    hit = {(m.end, m.length) for m in a.find_all("GCATATGC")}
    assert (7, 6) in hit


def test_non_acgt_base_breaks_a_match():
    a = Automaton(["GAATTC"])
    assert a.find_all("GAANTTC") == []
    assert len(a.find_all("GAATTCGAATTC")) == 2


def test_state_carries_across_calls():
    """Fixed context is fed first, then the optimised region continues from its state."""
    a = Automaton(["GAATTC"])
    state = a.feed("TTGAA")
    _, matches = a.scan("TTCGG", state)
    assert [(m.end, m.length) for m in matches] == [(3, 6)]  # 6-long match ends 3 bases in


def test_empty_pattern_set_is_a_valid_noop():
    a = Automaton([])
    assert a.find_all("ACGT" * 10) == []
    assert a.max_len == 0


@pytest.mark.parametrize(
    "patterns",
    [
        ["GAATTC", "CTCGAG"],
        ["GGTCTC", "CATATG"],
        ["GRATCY", "TTTT", "ACN"],
        ["AAAAAA", "AAAAAT"],
    ],
)
def test_matches_regex_oracle_on_random_text(patterns):
    rng = random.Random(1234)
    a = Automaton(patterns)
    for _ in range(200):
        text = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 80)))
        got = {(m.end, m.length) for m in a.find_all(text)}
        assert got == brute_force(text, patterns), text
