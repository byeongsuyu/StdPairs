stdPairs
================

stdPairs is a library of SageMath that gives symbolic computation over a monomial ideal of an affine (non-normal) semigroup ring. In particular, it provides associate prime ideals, multiplicity, and irredundant irreducible primary decomposition of a monomial ideal. It implements algorithms in [Matusevich and Yu, 2020](https://arxiv.org/abs/2005.10968). 

The source file can be found on the [github pages](https://github.com/byeongsuyu/StdPairs/). 

Minimum requirement
--------------
- [SageMath](https://www.sagemath.org) version 9.1 (Release Date: 2020-05-20) or higher.
- [4ti2 Package](https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/four_ti_2.html) in [SageMath](https://www.sagemath.org)
- [Macaulay2](http://www2.macaulay2.com/Macaulay2/) version 1.14 or higher.

Load Package
-------------
On command line
```
$ ./sage
┌────────────────────────────────────────────────────────────────────┐
│ SageMath version 9.1, Release Date: 2020-05-20                     │
│ Using Python 3.7.3. Type "help()" for help.                        │
└────────────────────────────────────────────────────────────────────┘
sage:load("[path]/stdPairs.spyx")
```
On jupyter attached to Sage
```
load("[path]/stdPairs.spyx")
```

User Manuel
-------------
See the section 2 of [Yu, 2020](https://arxiv.org/abs/2010.08903). 


Copyright & License
-

StdPairs is Copyright (c) 2020-2020 Byeongsu Yu. All Rights Reserved.

Permission to modify and redistribute is granted under the terms of the GPL 3.0 license. See the [LICENSE.txt](https://github.com/byeongsuyu/StdPairs/blob/main/LICENSE) file for the full license.
