Affine arithmetic
=================


Affine arithmetic (AA) is a computational model for numerical analysis based on the notion of **affine forms** - combinations of primitive variables which stand for sources of undertainty during the computation.

Affine arithmetic is meant to alleviate one of the major problems in `interval arithmetic <https://en.wikipedia.org/wiki/Interval_arithmetic>`_, namely the **dependency problem**: 

    *if an interval occurs several times in a calculation using parameters, and each occurrence is taken independently, then this can lead to an unwanted expansion of the resulting intervals*.

For example, the function :math:`f = x^x + x` is bounded by :math:`[-\frac{1}{4}, 2]` for :math:`x \in [-1, 1]`, but interval arithmetic ignoring the dependency between :math:`x^2` and :math:`x` effectively calculates the *infumum* and *supremum* of the function :math:`h(x,y) = x^2 + y` over :math:`x,y \in [-1, 1`.

Affine model
------------

In affine arithmetic, a quantity :math:`x` is represented as the following affine form:

.. math::

    x = x_0 + x_1 \epsilon_1 + ... + x_n \epsilon_n

where :math:`\epsilon_1, ..., \epsilon_n` are symbolic real variables whose values are unknown but assumed to lie in :math:`[-1, 1]`. Note that the number :math:`n` changes during the calculation.

In the case of a multivariate function :math:`f=(x_1,...,x_m)` the following affine forms are initialized:

.. math::
    \begin{align*}
        x_1 &= \frac{\overline{x}_1 + \underline{x}_1}{2} +  \frac{\overline{x}_1 - \underline{x}_1}{2} \epsilon_1\\
        \vdots\\
        x_m &= \frac{\overline{x}_m + \underline{x}_m}{2} +  \frac{\overline{x}_m - \underline{x}_m}{2} \epsilon_m\\
    \end{align*}

where :math:`[\underline{x}_k, \overline{x}_k]` is the domain of variable :math:`x_k`.

An affine form can be converted to an interval using the formula:

.. math::

    [x_0 - \Delta, x_0 + \Delta], \Delta = \sum_{i=1}^n |x_i|

Let :math:`I(x)` be the conversion of an affine form to a corresponding standard interval.

For two affine forms, :math:`x = x_0 + \sum_{i=1}^n x_i \epsilon_i` and :math:`y = y_0 + \sum_{i=1}^n y_i \epsilon_i` the following linear operations are defined:

.. math::
    
    \begin{align}
        x \pm y &= (x_0 \pm y_0) + \sum_{i=1}^n (x_i \pm y_i) \epsilon_i\\
        x \pm \alpha &= (x_0 \pm \alpha) + \sum_{i=1}^n x_i \epsilon_i\\
        \alpha x &= (\alpha x_0) + \sum_{i=1}^n (\alpha x_i)
    \end{align}

A nonlinear function :math:`f(x)` of an affine form is generally not able to be represented directly as an affine form. We must therefore consider a linear approximation of :math:`f` and a representation of the approximation error by introducing a new noise symbol :math:`\epsilon_{n+1}`.

Let :math:`X = I(x)` be the range of :math:`x`. For a nonlinear function :math:`f(x)`, a linear approximation in the form :math:`ax + b` will have a maximum approximation error :math:`\delta`:

.. math::

    \delta = \max_{x \in X} | f(x) - (ax + b) |

The result of the nonlinear operation can then be represented as follows:

.. math::
    \begin{align}
        f(x) &= ax + b + \delta \epsilon_{n+1}\\
             &= a(x_0 + x_1 \epsilon_1 + ... + x_n \epsilon_n) + b + \delta \epsilon_{n+1}
    \end{align}

Nonlinear binomial operations are calculated similarly.

Minima and maxima of multivariate functions
-------------------------------------------

Without loss of generality, we consider finding maxima of a two-dimensional function :math:`f(x_1, x_2)` in the box :math:`X^{(0)} = (X_1^{(0)}, X_2^{(0)}) = ([\underline{X_1^{(0)}}, \overline{X_1^{(0)}}], [\underline{X_2^{(0)}}, \overline{X_2^{(0)}}])`.

Notation
""""""""

For an interval :math:`J`, let the center and the width of :math:`J` be :math:`c(J)` and :math:`w(J)`, respectively.

For a box :math:`X`, let :math:`F_A(X)` be the range boundary of :math:`f` in :math:`X` obtained by applying *AA* and let the upper bound of :math:`I(F_A(X))` be :math:`\overline{F_A(X)}`.


Algorithm
"""""""""

1. Initialize two lists :math:`\mathcal{S}` and :math:`\mathcal{T}` in which to store boxes and range boundaries of :math:`f`:

    .. math::

        \mathcal{S} = \mathcal{T} = \emptyset

   and set stopping criteria :math:`\epsilon_r, \epsilon_b`.

2. **If** :math:`w(X_1^{(0)}) < w(X_2^{(0)})`, divide :math:`X(0)` into:
   
   .. math::

        \begin{align}
            X^{(1)} &=([\underline{X_1^{(0)}}, \overline{X_1^{(0)}}], [\underline{X_2^{(0)}}, c(X_2^{(0)})])\text{ and}\\
            X^{(2)} &=([\underline{X_1^{(0)}}, \overline{X_1^{(0)}}], [c(X_2^{(0)}), \overline{X_2^{(0)}}])
        \end{align}

   **otherwise** divide :math:`X(0)` into:

   .. math::

        \begin{align}
            X^{(1)} &=([\underline{X_1^{(0)}}, c(X_1^{(0)})], [\underline{X_2^{(0)}}, \overline{X_2^{(0)}}])\text{ and}\\
            X^{(2)} &=([c(X_1^{(0)}), \overline{X_1^{(0)}}], [\underline{X_2^{(0)}}, \overline{X_2^{(0)}}])
        \end{align}

3. Calculate :math:`F_A(X^{(1)})` and :math:`F_A(X^{(2)})`, then calculate :math:`\underline{f_{\max}^{(1)}}` and :math:`\underline{f_{\max}^{(2)}}`. The lower bound of the maxima is then given as :math:`\underline{f_{\max}} = \max(\underline{f^{(1)}_{\max}}, \underline{f^{(2)}_{\max})}` (use **Algorithm 1**).

4. **If** :math:`\overline{F_A(X^{(1)})} < \underline{f_{\max}}`, insert :math:`X^{(2)}` and :math:`F_A(X^{(2)})` into :math:`\mathcal{S}`. **Discard** :math:`X^{(1)}`.

   **If** :math:`\overline{F_A(X^{(2)})} < \underline{f_{\max}}`, insert :math:`X^{(1)}` and :math:`F_A(X^{(1)})` into :math:`\mathcal{S}`. **Discard** :math:`X^{(2)}`.

   **Otherwise** insert :math:`X^{(1)}`, :math:`F_A(X^{(1)})`, :math:`X^{(2)}`, :math:`F_A(X^{(2)})` into :math:`\mathcal{S}`.


5. **Repeat**:

   5.1. **If** :math:`S=\emptyset` go to **step 6**.

        **Otherwise**, find the box :math:`X^{(i)} \in \mathcal{S}` for which :math:`F_A(X^{(i)})` is largest. Select :math:`X^{(i)}` and :math:`F_A(X^{(i)})` as the box and the range boundary to be processed and remove them from :math:`\mathcal{S}`.

   5.2. Calculate :math:`\underline{f_{\max}^{(i)}}` (the candidates of :math:`f_{\max}`) by utilizing :math:`X^{(i)}` and :math:`F_A(X^{(i)})` and by applying **Algorithm 1**. Update :math:`\underline{f_{\max}}=\max\{{\underline{f_{\max}^{(i)}}}\}`

   5.3. Discard any box :math:`X` and range boundary :math:`F_A(X)` from :math:`\mathcal{S}` and :math:`\mathcal{T}` for which :math:`\overline{F_A(X)} < \underline{f_{\max}}`.

   5.4. Narrow :math:`X^{(i)}` down by utilizing :math:`X^{(i)}, F_A(X^{(i)})` and :math:`\underline{f_{\max}}` using **Algorithm 2**.

   5.5. **If** :math:`w(X_1^{(i)}) < w(X_2^{(i)})`, divide :math:`X^{(i)}` into:

       .. math::

            \begin{align}
                X^{(j)} &=([\underline{X_1^{(i)}}, \overline{X_1^{(i)}}], [\underline{X_2^{(i)}}, c(X_2^{(i)})])\text{ and}\\
                X^{(k)} &=([\underline{X_1^{(i)}}, \overline{X_1^{(i)}}], [c(X_2^{(i)}), \overline{X_2^{(i)}}])
            \end{align}


    **otherwise** divide :math:`X^{(i)}` into:

       .. math::

            \begin{align}
                X^{(j)} &=([\underline{X_1^{(i)}}, c(X_1^{(i)})], [\underline{X_2^{(i)}}, \overline{X_2^{(i)}}])\text{ and}\\
                X^{(k)} &=([c(X_1^{(i)}), \overline{X_1^{(i)}}], [\underline{X_2^{(i)}}, \overline{X_2^{(i)}}])
            \end{align}

    
  5.6. Calculate :math:`F_A(X^{(j)})` and :math:`F_A(X^{(k)})`.

         **If** :math:`\max_{1 \leq h \leq m} w(X_h^{(j)}) < \epsilon_r` **and** :math:`w(I(F_A(X^{(j)}))) < \epsilon_b`, insert :math:`X^{(j)}` and :math:`F_A(X^{(j)})` into :math:`\mathcal{T}`. **Otherwise**, insert :math:`X^{(j)}` and :math:`F_A(X^{(j)})` into :math:`\mathcal{S}`. Repeat the procedure for :math:`X^{(k)}` and :math:`F_A(X^{(k)})`.

6. Group the boxes left in :math:`\mathcal{T}` that have a common point with each other. Let the boxes that belong to the same group be :math:`Y^{(1)},...,Y^{(l)}`. The maximum value in the group is calculated as

.. math::

    \cup_{h=1}^l I(F_A(Y^{(h)}))

The point corresponding to the maximum value is calculated as:

.. math::

    \cup_{h=1}^l} Y^{(h)}

Repeat for all groups.
