Test StdPairs package
==========================

To be certain on functioning of the implemented methods and attributes of the package, one can test all examples in the documens directly using the command below.

.. code-block:: shell

   sage -t $SAGE_ROOT/local/lib/python3.9/site-packages/stdpairs/affinemonoid.py
   sage -t --long $SAGE_ROOT/local/lib/python3.9/site-packages/stdpairs/monomialideal.py
   sage -t $SAGE_ROOT/local/lib/python3.9/site-packages/stdpairs/properpair.py
   sage -t $SAGE_ROOT/local/lib/python3.9/site-packages/stdpairs/macaulay2_interface.py
   sage -t $SAGE_ROOT/local/lib/python3.9/site-packages/stdpairs/txt_to_object.py

where ``$SAGE_ROOT`` is the folder your sage was compiled or installed. Make sure that testing ``monomialideal.py`` file takes 2~3 minutes, unlike others. If one wants to test for more serious examples, we provide two test examples as below.

1. `Test overall <https://github.com/byeongsuyu/StdPairs/blob/551979519a41afeab3a12cc6723b928dacedce0f/tests/test_overall.sage>`_ executes all methods and attributes of the package and compares their return with pre-defined true computational result. It took a 15 minutes depending on your machine.

2. `Test indifference between Macaulay2 and StdPairs <https://github.com/byeongsuyu/StdPairs/blob/551979519a41afeab3a12cc6723b928dacedce0f/tests/test_standard_cover.sage>`_ executes ``MonomialIdeal.standard_cover()`` method over a (randomly selected) monomial ideal over a polynomial ring and compares its return with the same calculation from `Macaulay2 <http://www2.macaulay2.com/Macaulay2/>`_.

To test your installed package, please download above testing files, and type

.. code-block:: shell

   sage -t --long test_overall.sage

on the folder having tests files downloaded. For example, above one gives you messages as below.

.. code-block:: shell

   too many failed tests, not using stored timings
   Running doctests with ID 2021-03-14-19-03-06-6813f135.
   Git branch: develop
   Using --optional=4ti2,build,dochtml,homebrew,pip,sage,sage_spkg
   Doctesting 1 file.
   This is test suits for stdpairs.
   running AffineMonoid and its methods . . .  pass
   running txt_to_affinemonoid() . . . pass
   sage -t --long --random-seed=0 test_overall.sage
   running txt_to_monomialideal() . . . pass
   running MonomialIdeal and its methods . . .Calculate the standard cover of an ideal
   It takes a few minutes, depending on the system.
   Cover for 1  generator was calculated.  2  generators remaining. 
   Cover for 2  generators was calculated.  1  generators remaining. 
   Cover for 3  generators was calculated.  0  generators remaining. 
   Calculate the standard cover of an ideal
   It takes a few minutes, depending on the system.
   Cover for 1  generator was calculated.  1  generators remaining. 
   Cover for 2  generators was calculated.  0  generators remaining. 
   Calculate the standard cover of an ideal
   It takes a few minutes, depending on the system.
   Cover for 1  generator was calculated.  1  generators remaining. 
   Cover for 2  generators was calculated.  0  generators remaining. 
   Calculate the standard cover of an ideal
   It takes a few minutes, depending on the system.
   Cover for 1  generator was calculated.  1  generators remaining. 
   Cover for 2  generators was calculated.  0  generators remaining. 
   Calculate the standard cover of an ideal
   It takes a few minutes, depending on the system.
   Cover for 1  generator was calculated.  1  generators remaining. 
   Cover for 2  generators was calculated.  0  generators remaining. 
   Calculate the standard cover of an ideal
   It takes a few minutes, depending on the system.
   Cover for 1  generator was calculated.  2  generators remaining. 
   Cover for 2  generators was calculated.  1  generators remaining. 
   Cover for 3  generators was calculated.  0  generators remaining. 
   Calculate the standard cover of an ideal
   It takes a few minutes, depending on the system.
   Cover for 1  generator was calculated.  4  generators remaining. 
   Cover for 2  generators was calculated.  3  generators remaining. 
   Cover for 3  generators was calculated.  2  generators remaining. 
   Cover for 4  generators was calculated.  1  generators remaining. 
   Cover for 5  generators was calculated.  0  generators remaining. 
   Calculate the standard cover of an ideal
   It takes a few minutes, depending on the system.
   Cover for 1  generator was calculated.  1  generators remaining. 
   Cover for 2  generators was calculated.  0  generators remaining. 
   Calculate the standard cover of an ideal
   It takes a few minutes, depending on the system.
   Cover for 1  generator was calculated.  1  generators remaining. 
   Cover for 2  generators was calculated.  0  generators remaining. 
   Calculate the standard cover of an ideal
   It takes a few minutes, depending on the system.
   Cover for 1  generator was calculated.  2  generators remaining. 
   Cover for 2  generators was calculated.  1  generators remaining. 
   Cover for 3  generators was calculated.  0  generators remaining. 
   pass
   running prime_ideal() . . . pass
   running ProperPair class and its methods . . . pass
   running div_pairs() . . . pass
   running macaulay2() . . . pass
   running from_macaulay2() . . . pass
   running to_macaulay2() . . . pass
   StdPair package test was done! Everything works well.
   [0 tests, 0.00 s]
   ----------------------------------------------------------------------
   All tests passed!
   ----------------------------------------------------------------------
   Total time for all tests: 732.3 seconds
   cpu time: 0.0 seconds
   cumulative wall time: 0.0 seconds

