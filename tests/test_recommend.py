"""Pair suggestions: completeness, correctness of what each pair leaves on the protein, grouping."""

import os
from pathlib import Path

import pytest
from Bio.Seq import Seq

import p2d
from p2d.enzymes import list_cut_sites, plan_cut
from p2d.order import make_order
from p2d.recommend import group_by_outcome, rank_pairs
from p2d.vector import Vector

PROTEIN = "KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE"
SAMPLES = Path(os.environ["P2D_VECTOR_SAMPLES"]) if os.environ.get("P2D_VECTOR_SAMPLES") else None


@pytest.fixture(scope="module")
def vector():
    return p2d.load_bundled()


@pytest.fixture(scope="module")
def host():
    return p2d.load_host("ecoli_bl21")


def brute_force(vector, host, protein, cterm, window=None):
    """Every pair that plan_cut accepts, designs cleanly and orders cleanly, found the slow way."""
    sites = list_cut_sites(vector, window)
    found = {}
    for up in sites:
        for down in sites:
            if down.start < up.end:
                continue
            try:
                site = plan_cut(vector, up, down)
                plan = p2d.ConstructPlan.single_chain(protein)
                a = p2d.assemble(
                    p2d.design_insert(
                        plan,
                        vector,
                        site,
                        host,
                        strategy=p2d.CodonStrategy.MAX_CAI,
                        cterm=cterm,
                    )
                )
            except ValueError:
                continue
            if a.ok and make_order(a).ok:
                found[(up.id, down.id)] = a
    return found


def real_ends(vector, a, protein):
    """What the finished plasmid really makes, read from the vector's own start codon."""
    seq = a.seq[vector.expression_start :]
    translated = str(Seq(seq[: len(seq) - len(seq) % 3]).translate(to_stop=True))
    at = translated.index(protein)
    return translated[:at], translated[at + len(protein) :]


@pytest.mark.parametrize("cterm", [p2d.CTerm.STOP, p2d.CTerm.TAG, p2d.CTerm.VECTOR_FRAME])
def test_every_workable_pair_is_found_and_no_other(vector, host, cterm):
    options, tried = rank_pairs(vector, PROTEIN, host, cterm=cterm)
    expected = brute_force(vector, host, PROTEIN, cterm)
    assert {(o.up.id, o.down.id) for o in options} == set(expected)
    assert tried >= len(options) > 0


def test_what_each_pair_leaves_matches_independent_translation(vector, host):
    for cterm in (p2d.CTerm.STOP, p2d.CTerm.TAG):
        for option in rank_pairs(vector, PROTEIN, host, cterm=cterm)[0]:
            a = brute_force(vector, host, PROTEIN, cterm)[(option.up.id, option.down.id)]
            before, after = real_ends(vector, a, PROTEIN)
            assert (option.before, option.after) == (before, after), (
                option.up.label,
                option.down.label,
            )
            assert option.extra_residues == len(before) + len(after)
            assert option.plasmid_bp == len(a.seq) and option.ligation_ok


def test_n_terminal_tag_flag_means_vectors_own_residues_remain(vector, host):
    options, _ = rank_pairs(vector, PROTEIN, host)
    assert any(o.n_tag for o in options) and any(not o.n_tag for o in options)
    for o in options:
        assert o.n_tag == o.before.startswith("MGSSHHHHHH"), (o.up.label, o.before)


def test_tag_mode_reads_his_tag_for_every_pair(vector, host):
    for o in rank_pairs(vector, PROTEIN, host, cterm=p2d.CTerm.TAG)[0]:
        assert o.after.endswith("HHHHHH")
    for o in rank_pairs(vector, PROTEIN, host, cterm=p2d.CTerm.STOP)[0]:
        assert o.after == ""


def test_options_are_ordered_by_what_they_add(vector, host):
    options, _ = rank_pairs(vector, PROTEIN, host)
    keys = [o.sort_key for o in options]
    assert keys == sorted(keys)
    assert options[0].extra_residues == min(o.extra_residues for o in options)
    assert options[0].before == "MG" and not options[0].n_tag, "NcoI adds the least here"


def test_ranking_is_deterministic(vector, host):
    first = [(o.up.id, o.down.id) for o in rank_pairs(vector, PROTEIN, host)[0]]
    again = [(o.up.id, o.down.id) for o in rank_pairs(vector, PROTEIN, host)[0]]
    assert first == again


def test_limit_and_window_narrow_search(vector, host):
    everything, tried_all = rank_pairs(vector, PROTEIN, host)
    top, _ = rank_pairs(vector, PROTEIN, host, limit=3)
    assert [(o.up.id, o.down.id) for o in top] == [(o.up.id, o.down.id) for o in everything[:3]]
    narrow, tried = rank_pairs(vector, PROTEIN, host, window=(144, 180))
    assert tried < tried_all and all(o.up.start >= 144 and o.down.end <= 180 for o in narrow)


def test_n_tag_filter_partitions_options(vector, host):
    both, _ = rank_pairs(vector, PROTEIN, host)
    keep, _ = rank_pairs(vector, PROTEIN, host, n_tag=True)
    drop, _ = rank_pairs(vector, PROTEIN, host, n_tag=False)
    assert all(o.n_tag for o in keep) and not any(o.n_tag for o in drop)
    assert {(o.up.id, o.down.id) for o in keep} | {(o.up.id, o.down.id) for o in drop} == {
        (o.up.id, o.down.id) for o in both
    }


def test_groups_partition_options_and_each_is_one_outcome(vector, host):
    options, _ = rank_pairs(vector, PROTEIN, host)
    groups = group_by_outcome(options)
    members = [(o.up.id, o.down.id) for g in groups for o in g.pairs]
    assert sorted(members) == sorted((o.up.id, o.down.id) for o in options), (
        "no pair lost or doubled"
    )
    assert len({(g.before, g.after, g.n_tag, g.pad_bases) for g in groups}) == len(groups)
    for g in groups:
        assert {o.outcome for o in g.pairs} == {(g.before, g.after, g.n_tag, g.pad_bases)}
        assert g.best.sort_key == min(o.sort_key for o in g.pairs)
    bests = [g.best.sort_key for g in groups]
    assert bests == sorted(bests), "groups follow the ranking of their best pair"


def test_grouping_by_outcome(vector, host):
    options, _ = rank_pairs(vector, PROTEIN, host, cterm=p2d.CTerm.STOP)
    groups = group_by_outcome(options)
    assert len(groups) < len(options), "with a stop the 3' site does not change the protein"
    assert max(len(g.pairs) for g in groups) > 1


def test_empty_input_gives_no_groups():
    assert group_by_outcome([]) == []


def test_start_codon_is_required(vector, host):
    bare = Vector.from_sequence(vector.seq, name="bare")
    with pytest.raises(ValueError, match="start codon"):
        rank_pairs(bare, PROTEIN, host)


def test_window_with_too_many_pairs_must_be_narrowed(host):
    import random

    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(30_000))
    big = Vector.from_sequence(seq[:80] + "ATG" + seq[83:], name="big")
    big.expression_start = 80
    with pytest.raises(ValueError, match="narrow the window"):
        rank_pairs(big, PROTEIN, host)


def test_window_with_one_site_has_nothing_to_pair(vector, host):
    assert rank_pairs(vector, PROTEIN, host, window=(144, 150)) == ([], 0)


def test_real_pet28a_results_are_correct_and_distinct(host):
    from p2d.analysis import analyze_vector, confirm
    from p2d.vectorio import read_vector_file

    path = SAMPLES / "pET-28a_plus.dna" if SAMPLES else None
    if not path or not path.is_file():
        pytest.skip("pET-28a_plus.dna not available (set P2D_VECTOR_SAMPLES)")
    v = confirm(analyze_vector(read_vector_file(path)))
    window = (5060, 5260)
    options, tried = rank_pairs(v, PROTEIN, host, window=window)
    expected = brute_force(v, host, PROTEIN, p2d.CTerm.STOP, window)
    assert {(o.up.id, o.down.id) for o in options} == set(expected) and len(options) >= 50
    for o in options:
        before, after = real_ends(v, expected[(o.up.id, o.down.id)], PROTEIN)
        assert (o.before, o.after) == (before, after)
    groups = group_by_outcome(options)
    assert len(groups) >= 8 and len(groups) < len(options) / 3
    top = groups[0]
    assert top.before == "MG" and not top.n_tag and top.best.up.enzyme == "NcoI", (
        "NcoI drops the tag"
    )
    ndei = next(g for g in groups if g.best.up.enzyme == "NdeI")
    assert ndei.before == "MGSSHHHHHHSSGLVPRGSHM" and ndei.n_tag

    tagged, _ = rank_pairs(v, PROTEIN, host, window=window, cterm=p2d.CTerm.TAG)
    best = min(tagged, key=lambda o: (o.extra_residues, o.sort_key))
    assert best.up.enzyme == "NcoI" and best.down.enzyme == "XhoI" and best.after == "LEHHHHHH"


def test_added_parts_leave_vector_residues_alone(vector, host):
    from p2d.model import LINKERS, TAGS

    n, c = (TAGS["flag"], LINKERS["gs"]), (LINKERS["g4s"], TAGS["strep2"])
    plain, _ = rank_pairs(vector, PROTEIN, host)
    built, _ = rank_pairs(vector, PROTEIN, host, n_parts=n, c_parts=c)
    key = lambda o: (o.up.id, o.down.id)  # noqa: E731
    assert {key(o): o.outcome for o in plain} == {key(o): o.outcome for o in built}
    option = built[0]
    site = plan_cut(vector, option.up, option.down)
    chain = p2d.Chain("payload", PROTEIN, n_terminal=list(n), c_terminal=list(c))
    a = p2d.assemble(
        p2d.design_insert(p2d.ConstructPlan([chain]), vector, site, host, cterm=p2d.CTerm.STOP)
    )
    assert a.ok and "DYKDDDDKGS" + PROTEIN + "GGGGSWSHPQFEK" in a.protein
