Affine arithmetic
=================


Affine arithmetic (AA) is a computational model for numerical analysis based on the notion of **affine forms** - combinations of primitive variables which stand for sources of undertainty during the computation.

Affine arithmetic is meant to alleviate one of the major problems in `interval arithmetic <https://en.wikipedia.org/wiki/Interval_arithmetic>`_, namely the **dependency problem**: 

    *if an interval occurs several times in a calculation using parameters, and each occurrence is taken independently, then this can lead to an unwanted expansion of the resulting intervals*.

For example, the function :math:`f = x^x + x` is bounded by :math:`[-\frac{1}{4}, 2]` for :math:`x \in [-1, 1]`, but interval arithmetic ignoring the dependency between :math:`x^2` and :math:`x` effectively calculates the *infumum* and *supremum* of the function :math:`h(x,y) = x^2 + y` over :math:`x,y \in [-1, 1`.
