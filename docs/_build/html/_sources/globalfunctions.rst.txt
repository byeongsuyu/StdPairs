Global functions
============================


Decompose a difference between submonoids as submonoids
-----------------------------------------------------------

If we have two translations of submonoids :math:`b+\mathbb{N}G` and :math:`b+\mathbb{N}G'` where :math:`G' \supseteq G` as a face, then there is an algorithm to represent the difference :math:`(b+\mathbb{N}G) \setminus (b+\mathbb{N}G)'` as a union of finitely many translations of submonoids :math:`\bigcup_{i=1}^{m}a_{i}+\mathbb{N}F_{i}` where :math:`G \subseteq F_{i}` for any :math:`i \in \{1,2, \cdots, m\}`. This algorithm was implemented as ``pair_difference()``. More details on this algorithm can be found in `Yu, 2020 <https://arxiv.org/abs/2010.08903>`__ or `Matusevich and Yu, 2020 <https://arxiv.org/abs/2005.10968>`__.


.. currentmodule:: stdpairs

.. autofunction:: pair_difference


Find a prime ideal corresponding to a face
-----------------------------------------------------------

.. currentmodule:: stdpairs

.. autofunction:: prime_ideal


Division of pairs
-----------------------------------------------------------

.. currentmodule:: stdpairs

.. autofunction:: div_pairs



Interface using ``Macaulay2``
-------------------------------------

Using `MonomialSubalgebra <https://faculty.math.illinois.edu/Macaulay2/doc/Macaulay2-1.15/share/doc/Macaulay2/MonomialAlgebras/html/index.html>`__ class in `Macualay2 <https://faculty.math.illinois.edu/Macaulay2/>`__, one can translate ``MonomialIdeal`` object in this package into `MonomialSubalgebra <https://faculty.math.illinois.edu/Macaulay2/doc/Macaulay2-1.15/share/doc/Macaulay2/MonomialAlgebras/html/index.html>`__ object in `Macualay2 <https://faculty.math.illinois.edu/Macaulay2/>`__ via methods below.


.. currentmodule:: stdpairs

.. autofunction:: from_macaulay2
.. autofunction:: to_macaulay2


Save and load via ``string``
--------------------------------------------

``SageMath`` provide a global function writing and reading objects as a binary file by  `pickling <https://doc.sagemath.org/html/en/reference/misc/sage/misc/persist.html>`__. However, it may not work in case Sage code for the package is dramatically changed. To avoid such a catastrophy, we provides some global functions which may save and load your calculation about the monomial ideal objects easily.



.. currentmodule:: stdpairs

.. autofunction:: txt_to_affinemonoid
.. autofunction:: txt_to_monomialideal











