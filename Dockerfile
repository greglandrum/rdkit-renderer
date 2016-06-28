FROM debian:jessie
MAINTAINER Greg Landrum <greg.landrum@gmail.com>

# adapted from continuumio/miniconda3
RUN apt-get update --fix-missing && apt-get install -y \
  wget \
  bzip2 \
  ca-certificates \
  libc6

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.0.5-Linux-x86_64.sh && \
    /bin/bash /Miniconda3-4.0.5-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-4.0.5-Linux-x86_64.sh

ENV PATH /opt/conda/bin:$PATH
ENV LANG C

# actually do the conda install
RUN conda config --add channels  https://conda.anaconda.org/greglandrum
RUN conda install -y gunicorn flask cairo_nox nomkl
RUN conda install -y rdkit boost=1.56
COPY . /src

EXPOSE 8000
WORKDIR "/src"
ENTRYPOINT ["/opt/conda/bin/gunicorn", "--bind=0.0.0.0:8000"]
CMD ["wsgi:app"]
