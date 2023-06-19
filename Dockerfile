###################
# STAGE 1: builder
###################

# Build currently doesn't work on > Java 11 (i18n utils are busted) so build on 8 until we fix this
FROM adoptopenjdk/openjdk8:x86_64-debianslim-jre8u345-b01 as builder

WORKDIR /app/source

ENV PATH="$PATH:/opt/conda/bin"
ENV FC_LANG en-US
ENV LC_CTYPE en_US.UTF-8

# bash:    various shell scripts
# wget:    installing lein
# git:     ./bin/version
# make:    backend building
# gettext: translations
RUN apt-get update && apt-get install -y coreutils bash git wget make gettext

# lein:    backend dependencies and building
ADD ./bin/lein /usr/local/bin/lein
RUN chmod 744 /usr/local/bin/lein

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_22.11.1-1-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p /opt/conda
RUN /opt/conda/bin/conda install -c conda-forge -c bioconda mamba
RUN /opt/conda/bin/mamba install -c conda-forge -c bioconda -y python=3.9 r-base=3.6.3 r-renv blas lapack cxx-compiler

ADD ./bin/quartet-protqc-report /opt/conda/bin/quartet-prot-report
# Install report locally instead of remote to install the latest version.
ADD report /report
RUN /opt/conda/envs/venv/bin/pip install /report

ADD ./resources/bin/protqc.sh /opt/conda/bin/protqc.sh
ADD ./resources/renv /opt/conda/renv
ADD ./resources/renv.lock /opt/conda/renv.lock
# Install protqc locally instead of remote to install the latest version.
# ADD ./build/Rprofile /opt/conda/etc/Rprofile
# RUN /opt/conda/envs/venv/bin/Rscript /opt/conda/etc/Rprofile

# Disable cache to install all packages into the conda environment.
COPY protqc /protqc
RUN /opt/conda/envs/venv/bin/Rscript -e 'renv::activate("/opt/conda");renv::restore();renv::install("/protqc")'

# install dependencies before adding the rest of the source to maximize caching
# backend dependencies
ADD project.clj .
RUN lein deps

# add the rest of the source
ADD . .

# build the app
RUN lein uberjar

# ###################
# # STAGE 2: runner
# ###################

FROM adoptopenjdk/openjdk8:x86_64-debianslim-jre8u345-b01 as runner

LABEL org.opencontainers.image.source https://github.com/chinese-quartet/quartet-protqc-report.git

ENV PATH="$PATH:/opt/conda/bin"
ENV PYTHONDONTWRITEBYTECODE=1
ENV FC_LANG en-US
ENV LC_CTYPE en_US.UTF-8

RUN apt-get update && apt-get install -y coreutils bash git wget make gettext
RUN echo "**** Install dev packages ****" && \
    apt-get update && \
    apt-get install -y curl && \
    \
    echo "*** Install common development dependencies" && \
    apt-get install -y libmariadb-dev libxml2-dev libcurl4-openssl-dev libssl-dev && \
    \
    echo "**** Cleanup ****" && \
    apt-get clean


WORKDIR /data

COPY --from=builder /opt/conda /opt/conda
COPY --from=builder /app/source/target/uberjar/quartet-protqc-report*.jar /quartet-protqc-report.jar

# Copy all the installed R packages from the builder.
COPY --from=builder /root/.local/share/renv /root/.local/share/renv

# Run it
ENTRYPOINT ["quartet-protqc-report"]