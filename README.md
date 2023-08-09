# **IMPALIB**
=============
Table of contents
-----------------

* [Introduction](#introduction)
* [Application](#application)
* [Supported Constraints](#supported-constraints)
* [Code Parameters](#code-parameters)
* [Usage](#usage)
* [Requirements and Installation](#requirements-and-installation)
* [License](#license)
### **Introduction**
------------

Would you like to solve optimization problems using message-passing algorithms? And would you like perform Belief Propagation in C++ and python? Then IMPALIB (Iterative Message Passing Algorithm Library) is exactly for you.

IMPALIB:
- can be used in three different ways
    1. A small header-only library written in modern and pure C++
    2. A pure python code
    3. A C++ code with a python wrapper
- is very easy to integrate and use
- is self contained and depends only [cnpy](https://github.com/rogersce/cnpy) library for unit testing - also header-only library which is used for reading and writing numpy files
- supports inference for optimization problems with various possible constraints
- results in a much faster solution than Google OR-Tools
