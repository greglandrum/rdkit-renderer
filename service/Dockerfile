FROM debian:jessie
LABEL maintainer="Greg Landrum <greg.landrum@t5informatics.com>"

# adapted from continuumio/miniconda3
RUN apt-get update --fix-missing && apt-get install -y \
  wget \
  bzip2 \
  ca-certificates \
  libc6 \
  xorg-dev

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
  wget --quiet -O installer.sh https://github.com/conda-forge/miniforge/releases/download/4.10.1-0/Mambaforge-4.10.1-0-Linux-x86_64.sh && \
  /bin/bash /installer.sh -b -p /opt/conda && \
  rm installer.sh

ENV PATH /opt/conda/bin:$PATH
ENV LANG C

# actually do the conda install
RUN conda config --set auto_update_conda false
RUN mamba install -y gunicorn flask nomkl flasgger
RUN conda install -y rdkit

ARG FORCE_REBUILD=unknown

RUN FORCE_REBUILD=${FORCE_REBUILD} mamba install -y rdkit
COPY . /src

EXPOSE 8000
WORKDIR "/src"
ENTRYPOINT ["/opt/conda/bin/gunicorn", "--bind=0.0.0.0:8000"]
CMD ["wsgi:app"]
