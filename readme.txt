Conformal Blocks Software Package Prototype V0.2
------------------------------------------------

IMPORTANT:
*Please do not distribute without my permision
*This software is provided without any guarantees that the results it produces are
correct.  You should check all computations using Swinarski's program as a reference
before accepting them as true
*This software will have bugs, and I cannot guarantee they will not harm your computer

Setup
-----
The program has been converted to a Cython module and needs to be compiled before use.  Follow these
steps to setup the program.  You will need Sage installed on your computer.

1.  Extract the folder "ConformalBlocks" wherever you like on your computer.
2.  Open a terminal and change directory to the "ConformalBlocks" folder.
3.  Execute the following command:

$ sage setup.py build_ext --inplace

This builds a Python extension which you can now use with Sage.

Using the program in Sage
-------------------------
To use the program in Sage, run Sage in the "ConformalBlocks" folder.  Once started
import the conformal blocks package by execute the following command:

sage: import fusion_prod.cbd as cbd

To see the list of objects in the module, type "cbd." then press [TAB].  Computations are done
using these objects.

Example Usage
-------------
Two basic things the package can calculate are the tensor product and fusion product
of two representations.  The rep's are inputted in terms of the fundamental weight coordinates
of its highest weight.

sage: import fusion_prod.cbd as cbd
sage: liealg = cbd.TypeALieAlgebra(3)
sage: liealg.tensor([1,2,0],[3,1,1])
{[0, 2, 3]: 1,
 [0, 3, 1]: 1,
 [1, 0, 4]: 1,
 [1, 1, 2]: 2,
 [1, 2, 0]: 1,
 [1, 3, 2]: 1,
 [1, 4, 0]: 1,
 [2, 0, 1]: 1,
 [2, 1, 3]: 2,
 [2, 2, 1]: 3,
 [2, 4, 1]: 1,
 [3, 0, 2]: 2,
 [3, 1, 0]: 2,
 [3, 2, 2]: 2,
 [3, 3, 0]: 2,
 [4, 0, 3]: 1,
 [4, 1, 1]: 3,
 [4, 3, 1]: 1,
 [5, 0, 0]: 1,
 [5, 1, 2]: 1,
 [5, 2, 0]: 1,
 [6, 0, 1]: 1}
sage: liealg.fusion([1,2,0],[3,1,1],5)
{[0, 2, 3]: 1,
 [0, 3, 1]: 1,
 [1, 0, 4]: 1,
 [1, 1, 2]: 2,
 [1, 2, 0]: 1,
 [1, 4, 0]: 0,
 [2, 0, 1]: 1,
 [2, 2, 1]: 1,
 [3, 0, 2]: 1,
 [3, 1, 0]: 1,
 [5, 0, 0]: 0}

The program can also calculate ranks and divisors of conformal blocks.  The syntax to
create a symmetric bundle is the following:

SymmetricConformalBlocksBundle(<Lie Algebra>, <Weight>, <Number of points>, <Level>)

For example:

sage: V = cbd.SymmetricConformalBlocksBundle(liealg, [1,2,1], 9, 4)
sage: V.getRank()
41412
sage: V.get_norm_sym_divisor_ray()
[151081L, 189833L, 229020L]

(Ignore the L's)


Running Scripts
---------------

Two example scripts are included.  Scripts allow you to do more complicated calculations.
To run them follow these steps:

1.  Open a terminal and change directory to "ConformalBlocks".
2.  Run sage.
3.  Enter the command ``load("experiments/<Script>.py")" to run the script.


Change Log
----------
V0.1: Initial version
V0.2: Reorganization and Cython conversion
- Code has been significantly reorganized to allow easier interactive use.  In particular the program is
  now contained in a single module, and most methods now accept tuples or lists of integers in place of
  Weight or IrrRep objects
- The module is now a Cython module

