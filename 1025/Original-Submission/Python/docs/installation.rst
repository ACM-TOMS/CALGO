Installation
============

PARyOpt requires Python 3.5 or above, NumPy, and SciPy for basic functionality.
Paramiko is required for cost functions evaluated on remote machines
(HPC clusters). Matplotlib is used for visualization for these examples,
but is not required.

These all will be installed when you do a

.. code-block:: bash

   pip install paryopt

(if you are on Ubuntu, you may need to do ``sudo apt-get install python3-pip``
 and use ``pip3`` here instead!)

Or, if you prefer an Anaconda environment:

.. code-block:: bash

   conda create -n paryopt python=3.5 numpy scipy matplotlib paramiko
   activate paryopt

Or, if you are using a manual download:

.. code-block:: bash

    tar -xvf paryopt-1.0.1.tar.gz
    cd paryopt/
    python3.5 setup.py install

If you want a virtual environment (preferred), do:

.. code-block:: bash

    tar -xvf paryopt-1.0.1.tar.gz
    cd paryopt/
    python3.5 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
    pip install setup.py
