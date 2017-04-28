Conformal Blocks Software Package Prototype
-------------------------------------------

IMPORTANT:
*This software is in development and is provided without any guarantees that the results
it produces are correct.  You should check all computations using David Swinarski's program as a
reference before accepting them as true.

Using the program in Sage
-------------------------
To use the program in Sage, run Sage in the "conformal-blocks-master" folder.  Once started
import the conformal blocks package by execute the following command:

sage: import conformal_blocks.cbbundle as cbbundle

To see the list of objects in the module, type "cbbundle." then press [TAB].  Computations are done
using these objects.

Using the program with IPython
------------------------------
An alternative environment is to use IPython running on top of a Python2 runtime of your
choice.  Running on a recent version of the PyPy JIT runtime will result in a 3-8 times
speed increase (depending on the particular computations).  Usage in IPython is currently
identical to usage in Sage, but future features might not be available outside of the
Sage environment.

Example Usage
-------------
Two basic things the package can calculate are the tensor product and fusion product
of two representations.  The rep's are inputted in terms of the fundamental weight coordinates
of its highest weight.

sage: import conformal_blocks.cbbundle as cbbundle
sage: liealg = cbbundle.TypeALieAlgebra(3)
sage: liealg.tensor((1,2,0),(3,1,1))
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
sage: liealg.fusion((1,2,0),(3,1,1),5)
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

sage: V = cbbundle.SymmetricConformalBlocksBundle(liealg, (1,2,1), 9, 4)
sage: V.get_rank()
41412
sage: V.get_norm_sym_divisor_ray()
[151081L, 189833L, 229020L]


Running Scripts
---------------

Three example scripts are included.  Scripts allow you to do more complicated calculations.
To run them follow these steps:

1.  Open a terminal and change directory to "conformal-blocks-master".
2.  Start Sage or Ipython
3.  In Sage, enter the command

    sage: load("experiments/<Script>.py")

    to run the script.  In IPython enter the command:

    : run "experiments/<Script>.py"

