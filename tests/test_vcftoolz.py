# -*- coding: utf-8 -*-

"""
test_vcftoolz
----------------------------------

Tests for `vcftoolz` module.
"""

from vcftoolz import vcftoolz


def test_get_unique_set_elements():
    """Test: Given a list of sets, return a new list of sets with the overlapping
    elements removed, leaving only the set elements that are unique to each
    set.
    """
    unique_sets = vcftoolz.get_unique_set_elements([set()])
    assert(unique_sets == [set()])

    unique_sets = vcftoolz.get_unique_set_elements([set(), set()])
    assert(unique_sets == [set(), set()])

    unique_sets = vcftoolz.get_unique_set_elements([set(), set(), set()])
    assert(unique_sets == [set(), set(), set()])

    unique_sets = vcftoolz.get_unique_set_elements([{100, 101, 1, 2}, {100, 101, 3, 4}, {100, 101, 5, 6}])
    assert(unique_sets == [{1, 2}, {3, 4}, {5, 6}])

    unique_sets = vcftoolz.get_unique_set_elements([{100, 101}, {100, 101}, {100, 101}])
    assert(unique_sets == [set(), set(), set()])

    unique_sets = vcftoolz.get_unique_set_elements([{"1"}])
    assert(unique_sets == [{"1"}])

    unique_sets = vcftoolz.get_unique_set_elements([{"10", "11"}, {"01", "11"}])
    assert(unique_sets == [{"10"}, {"01"}])

    unique_sets = vcftoolz.get_unique_set_elements([{"100", "110", "101", "111"}, {"110", "010", "111", "011"}, {"101", "111", "011", "001"}])
    assert(unique_sets == [{"100"}, {"010"}, {"001"}])


def test_get_missing_set_elements():
    """Test: Given a list of sets, return a new list of sets containing the elements that are
    missing from each set but are present in more than one other set.
    """
    missing_sets = vcftoolz.get_missing_set_elements([set()])
    assert(missing_sets == [set()])

    missing_sets = vcftoolz.get_missing_set_elements([set(), set()])
    assert(missing_sets == [set(), set()])

    missing_sets = vcftoolz.get_missing_set_elements([set(), set(), set()])
    assert(missing_sets == [set(), set(), set()])

    missing_sets = vcftoolz.get_missing_set_elements([{100, 101}, {100, 101}, {100, 101}])
    assert(missing_sets == [set(), set(), set()])

    missing_sets = vcftoolz.get_missing_set_elements([{"1"}])
    assert(missing_sets == [set()])

    missing_sets = vcftoolz.get_missing_set_elements([{"10", "11"}, {"01", "11"}])
    assert(missing_sets == [set(), set()])

    missing_sets = vcftoolz.get_missing_set_elements([{"100", "110", "101", "111"}, {"110", "010", "111", "011"}, {"101", "111", "011", "001"}])
    assert(missing_sets == [{"011"}, {"101"}, {"110"}])


def test_get_set_intersections():
    """Test: Given a list of sets, return a new list of sets with all the possible
    mutually exclusive overlapping combinations of those sets.  Another way
    to think of this is the mutually exclusive sections of a venn diagram
    of the sets.  If the original list has N sets, the returned list will
    have (2**N)-1 sets.
    """
    intersections = vcftoolz.get_set_intersections([{"100", "110", "101", "111"}, {"110", "010", "111", "011"}, {"101", "111", "011", "001"}])
    assert(intersections == [('111', {'111'}), ('011', {'011'}), ('101', {'101'}), ('001', {'001'}), ('110', {'110'}), ('010', {'010'}), ('100', {'100'})])
