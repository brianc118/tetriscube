# TetrisCube

Solution to find all 9839 solutions of the 4x4 Tetris Cube. Extended to work with lxmxn volumes with arbitrarily defined pieces.

## Solution breakdown

1. 
	a) Reduce problem to variant of set-cover with the additional constraint that each piece needs to be used exactly once.

	b) Note we can weaken the constraint to that each piece needs to be used at least once since the solved is trying to minimize number of pieces (as in set-cover)

	c) The universe set, U, is the set of all points in the cube unioned with the set {a_i}^n for the n vertices. Then we create a collection S with the ith piece corresponding to the set { points } union {a_i}. Thus to cover U we require each piece to be used at least once and the pieces to cover all points of the cube.

2. To solve the set-cover problem we use linear programming (in fact ILP). There is no problem objective; everything is defined by equality constraints.