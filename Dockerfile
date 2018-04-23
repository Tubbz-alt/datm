#
# Deterministic Arrival Time Monitor Dockerfile
#
FROM centos:7
MAINTAINER Abdullah P Rashed Ahmed <apra@slac.stanford.edu>

# Install yum packages
RUN yum clean all && yum -y install bzip2.x86_64 libgomp.x86_64 telnet.x86_64 gcc-c++ strace curl

# Install Miniconda2 and add it to the path
RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
RUN chmod +x /tmp/miniconda.sh
RUN bash /tmp/miniconda.sh -bfp /usr/local
RUN rm -rf /tmp/miniconda.sh
RUN conda update --yes conda
RUN conda install --yes python=2
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh

# Conda install psana
RUN conda install --yes --channel lcls-rhel7 psana-conda
RUN conda install --yes --channel conda-forge "mpich>=3" mpi4py h5py pytables libtiff=4.0.6 
RUN rm -rf /opt/conda/lib/python2.7/site-packages/numexpr-2.6.2-py2.7.egg-info

# Install datm requirements from the master branch
RUN conda install --yes --channel conda-forge `curl https://raw.githubusercontent.com/slaclab/datm/master/requirements.txt`

# Install additional packages
RUN conda install --yes --channel conda-forge ipython ipdb 

# Clean up
RUN conda clean --all --yes
