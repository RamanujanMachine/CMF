from sympy.abc import x, y, z

from ramanujantools import Position


def test_add():
    p1 = Position({x: 17, y: 3})
    p2 = Position({y: 4, z: 3})
    assert {x: 17, y: 7, z: 3} == p1 + p2


def test_iadd():
    p = Position({x: 17, y: 3})
    p += Position({y: 4, z: 3})
    assert {x: 17, y: 7, z: 3} == p


def test_isub():
    p = Position({x: 17, y: 3})
    p -= Position({y: 4, z: 3})
    assert {x: 17, y: -1, z: -3} == p


def test_sub():
    p1 = Position({x: 17, y: 3})
    p2 = Position({y: 4, z: 3})
    assert {x: 17, y: -1, z: -3} == p1 - p2


def test_mul():
    p = Position({x: 1, y: -2, z: 3})
    expected = {key: -7 * value for key, value in p.items()}
    assert expected == -7 * p
    assert expected == p * -7


def test_neg():
    p = Position({x: 1, y: -2, z: 3})
    assert {key: -value for key, value in p.items()} == -p


def test_longest():
    p = Position({x: 1, y: 2, z: -3})
    assert 3 == p.longest()


def test_longest_empty():
    p = Position()
    assert 0 == p.longest()


def test_shortest():
    p = Position({x: 1, y: 2, z: -3})
    assert 1 == p.shortest()


def test_shortest_empty():
    p = Position()
    assert 0 == p.longest()


def test_signs():
    p = Position({x: 1, y: 2, z: -3})
    assert Position({x: 1, y: 1, z: -1}) == p.signs()


def test_signs_empty():
    assert Position() == Position().signs()