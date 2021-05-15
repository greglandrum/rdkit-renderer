FROM condaforge/mambaforge
LABEL maintainer="Greg Landrum <greg.landrum@t5informatics.com>"

ENV LANG C

# actually do the conda install
RUN conda config --set auto_update_conda false
RUN mamba install -y gunicorn flask nomkl flasgger
ARG FORCE_REBUILD=unknown

RUN FORCE_REBUILD=${FORCE_REBUILD} mamba install -y rdkit
COPY . /src

EXPOSE 8000
WORKDIR "/src"
ENTRYPOINT ["/opt/conda/bin/gunicorn", "--bind=0.0.0.0:8000"]
CMD ["wsgi:app"]
