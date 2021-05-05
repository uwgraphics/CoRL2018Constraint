Sympy Cython Module Creator
---------------------------

**symcymod** allows easy conversion of multiple sympy expressions
into Cython modules and c code. It makes use of sympy's codegen
utility to create cython files and corresponding code and header
files.

Requirements
------------
**sympy**, **cython**


How to use symcymod?
--------------------

Provide a list of tuples (**name**,**sympy expression**,**argument order**)
to create a module with the c library.

Note: The generated methods are suffixed with a **_c**.

For an example, please see the end of the file **symcymod.py**.
