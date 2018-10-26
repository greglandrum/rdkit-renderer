FROM debian:jessie
MAINTAINER Greg Landrum <greg.landrum@t5informatics.com>

# adapted from continuumio/miniconda3
RUN apt-get update --fix-missing && apt-get install -y \
  wget \
  bzip2 \
  ca-certificates \
  libc6

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
     wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh && \
     /bin/bash /Miniconda3-4.5.4-Linux-x86_64.sh -b -p /opt/conda && \
     rm Miniconda3-4.5.4-Linux-x86_64.sh

ENV PATH /opt/conda/bin:$PATH
ENV LANG C

# actually do the conda install
#RUN conda config --add channels  https://conda.anaconda.org/greglandrum
RUN conda install -y gunicorn flask nomkl
RUN conda config --append channels conda-forge
RUN conda install -y flasgger
RUN conda config --add channels rdkit
RUN conda install -y nox cairo

ARG FORCE_REBUILD=unknown

RUN FORCE_REBUILD=${FORCE_REBUILD} conda install -y rdkit
COPY . /src

EXPOSE 8000
WORKDIR "/src"
ENTRYPOINT ["/opt/conda/bin/gunicorn", "--bind=0.0.0.0:8000"]
CMD ["wsgi:app"]
