==================================
Deterministic Arrival Time Monitor  
==================================

Repo for the deterministic approach for the arrival time monitor.

Requirements
============

If direct access to ``psana`` is needed, then ``python 2.7`` will be needed,
otherwise, use ``python >3.4``.

There some additional packages that will need to be installed listed in the
``requirements.txt`` file. If installing using pip: ::

  $ pip install -r requirements.txt

If using conda: ::

  $ conda install `cat requirements.txt` -c conda-forge

Development
===========

For development on the deterministic algorithm, first clone the repo into the
desired directory: ::

  $ git clone https://github.com/ryancoffee/TimeTool_deterministic.git

Install the requirements using one of the methods listed above, and then add the
project to your python path by running the ``setup.py`` file: ::

  $ python setup.py develop

This will install the pcakage in development mode. To remove it, simply run the
same command with an additional ``--uninstall`` flag: ::

  $ python setup.py develop --uninstall
