FROM debian:jessie
MAINTAINER Greg Landrum <greg.landrum@gmail.com>

# adapted from continuumio/miniconda3
RUN apt-get update --fix-missing && apt-get install -y \
  wget \
  bzip2 \
  ca-certificates \
  libc6 \
  libglib2.0 \
  libxext6 \
  libsm6 \
  libxrender1 \
  git \
  curl \
  grep \
  sed \
  dpkg \
  libcairo2 

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-3.19.0-Linux-x86_64.sh && \
    /bin/bash /Miniconda3-3.19.0-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-3.19.0-Linux-x86_64.sh

RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean


ENV PATH /opt/conda/bin:$PATH
ENV LANG C.UTF-8

RUN conda config --add channels  https://conda.binstar.org/greglandrum
RUN conda install -y rdkit gunicorn flask cairo cffi

COPY . /src

EXPOSE 8000
WORKDIR "/src"
ENTRYPOINT ["/opt/conda/bin/gunicorn", "--bind=0.0.0.0:8000"]
CMD ["wsgi:app"]
