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


Docker
======

The dockerfile included in the repo contains all the necessary components to run
``datm``. If running this from ``daq-tst-dev05``, a slight modification needs to
be made to the Dockerfile. Simply remove the colon in line 4 so that it reads as
the following: ::

  3| #
  4| FROM centos7
  5| MAINTAINER ...

If the docker image has not been built already, ``cd`` into this directory and
then run the following command: ::

  $ docker build --rm -t datm/base .

Once it finishes (this can take some time), you can launch the container
normally with the following: ::

  $ docker run --rm -ti datm/base


Development
===========

For development on the deterministic algorithm, first clone the repo into the
desired directory: ::

  $ git clone https://github.com/slaclab/datm.git

Install the requirements using one of the methods listed above, and then add the
project to your python path by running the ``setup.py`` file: ::

  $ python setup.py develop

This will install the pcakage in development mode. To remove it, simply run the
same command with an additional ``--uninstall`` flag: ::

  $ python setup.py develop --uninstall
