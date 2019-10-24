This is the consequent biological sequences tool package.

Command-line Tools
******************

**align_check.py** -- sequence alignment with statiscal significance
**align_check_mp.py** -- same as above with multiprocessing support
**dot_plot.py** -- makes comparison dot plots of two protein sequences


Modules
*******

**matrix.py** -- read and write scoring matrices
**sequence.py** -- read, write and retrieve biological sequences
**alignment.py** -- routines implementing various alignment algorithms


Installation
------------

1. Create a python virtual environment::

   $ python3 -m venv <name>

2. Activate the virtual environment::

   $ . <name>/bin/activate

3. Clone the repo::

   $ git clone https://github.com/pjmartel/consequent

4. Install with pip::

   $ cd consequent ; pip install .


Usage
----

The command-line tools are available on the virtual environent exec path::
   $ dot_plot.py --help

Import the modules from Python codo::
   import consequent.matrix

After using, use::
   deactivate
to exit the virtual environment.
