# libtropicana

`libtropicana` is a library designed to find regular simplicial subdivision 
of lattice convex polytopes and also compute normalized volume as a byproduct.
It is written completely in C++ with optional interface for leveraging 
`spBLAS` (Sparse `BLAS`) routines.
`libtropicana` is open source software. 
Users may freely distribute its source under the terms of the LGPL license.

## Background

_Convex polytopes_ are higher dimensional generalizations of two dimensional
convex polygons or three dimensional polyhedrons that school children are 
usually familiar with.
They are geometric objects with flat sides that have the property that
the line segment connecting any two points in this object must be entirely
in this object as well (convex).
The mathematical definition is more technical --- A convex polytopes is a
finite intersection of half spaces that is bounded.
Equivalently, a convex polytope can also be defined as the convex hull
of a finite set of points (the smallest convex set containing them).
It is called a _lattice_ convex polytope if it can be realized as the
convex hull of points that only have integer coordinates.

In the study of convex polytopes it is often required to divide a 
convex polytope into smaller pieces that are easier to understand.
For many purpose (e.g. volume computation) it is convenient to divide
a convex polytope into a collection of simplicies which would be called
a _simplicial subdivision_.

## How it works

`libtropicana` is based on a pivoting algorithm similar to the core
of `lrs`.
But unlike `lrs`, which puts a special emphasis on memory efficiency
and accuracy, `libtropicana` focuses on speed (potentially at the 
expense of higher memory consumption).



