
.. _gquest_install:


Installation
====================================================

Git
----------------------------------------------------

The repository lives at  https://git.mccullerlab.com/gquest/gquest-design

Clone it from git. Use the https link if you don't intend to upload. Use the git+ssh link if you want to be more series.

.. code-block:: bash

  git clone https://git.mccullerlab.com/gquest/gquest-design.git

or, if you have your ssh keys set up and want push powers.

.. code-block:: bash

  git clone ssh://git@git.mccullerlab.com:2224/gquest/gquest-design.git


Conda
----------------------------------------------------

There are two stages of installing the gquest design repository. First is setting up a virtual environment (generally conda) and installing the dependencies. The second is installing the repository itself in developer mode, so that its code is available to the tests and notebooks.

Install conda or miniconda and make sure it is activated. The command

.. code-block:: bash

  conda

should exist in your bash if conda is functioning. I would then recommend creating a dedicated environment. If you have an existing environment that you like with all of your favorite packages, make a clone of it by replacing base with your environment name. If not, then either omit the clone argument or keep base as the target.

.. code-block:: bash 

   conda create --name gquest --clone base
   conda activate gquest

Some conda environments don't ship with pip, but we will want the conda-ified pip so that pip installs to the conda envrionment rather than the system or user environment, so make sure that it is installed. Additionally, ensure that `conda-forge <https://conda-forge.org/docs/user/introduction.html>`_ is active, to get the latest packages .

.. code-block:: bash 

  conda config --add channels conda-forge
  conda install mamba pip

In the above, the "mamba" package is also included. After installing it, you can change most commands from conda to mamba and they will be much faster. As far as I have seen, mamba is a drop-in replacement and does not cause problems.

Dependencies through conda
++++++++++++++++++++++++++++++++++++++++++++++++++++

Now install as many of the dependencies using conda as possible. The remainder will be done through pip in the next step.

.. code-block:: bash 

  conda install numpy scipy matplotlib pytest pytest-watch pytest-html jupyter jupytext h5py tabulate setuptools_scm networkx

You might want to include "git" in that set if your system doesn't already have it.

Development Environment
----------------------------------------------------

There are a number of in-development libraries from LM's set of "wield" packages. These should be installed through git. They simplest way to get started is to install them using the requirements file.

.. code-block:: bash 

   cd gquest-design
   pip install -r requirements.txt

.. code-block:: bash 

   cd gquest-design
   pip install -e .


now it should be working

testing
-------------------------------------------------------

You can check that everything worked correctly by running the tests with pytest.


.. code-block:: bash 

   cd gquest-design
   pytest

This should collect and run all of the tests. You can explore then outputs in the test_results folder that appears. Near each test file, there will also be a test_results folder that is local to the test file. In this way, the tests behave slightly like notebooks.

If you run

.. code-block:: bash 

   cd gquest-design
   ptw -- path/to/test_thing.py -k test_name 

Then you can run a single test, and the "ptw" pytest-watch command will run the test any time you modify the test file, or any other code. In this way you can edit and watch tests very rapidly, without needing to click or use a mouse. Open a pdf viewer to watch plots update as the test runs, or stare at the console output.

Installing wield Dependencies (optional)
-------------------------------------------------------

If you look at the install targets in the requirements.txt you'll find severall packages installed through git. Those are installed to some conda folder. If instead you want to install them locally, so that you can test and edit them, then you can do that. In this case, let's install wield-model, since it can do the beam propagation and matching calculations.


.. code-block:: bash 

   git clone https://github.com/wieldphysics/wield-model.git
   cd wield-model
   pip install -e .

and now that specific folder is available for python to use. If you modify it, then you can see the changes.

We can do similarly for gwinc

.. code-block:: bash 

   git clone https://git.ligo.org/lee-mcculler/pygwinc.git
   cd pygwinc
   git checkout superQKall
   pip install -e .
   python -m gwinc Aplus

These are the https git interfaces. You should go to the respective repository sites and use the git+ssh interface and fork the repository if you want to edit and submit pull requests.



