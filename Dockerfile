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
RUN apt-get update && apt-get install -y coreutils bash git wget make gettext g++

# lein:    backend dependencies and building
ADD ./bin/lein /usr/local/bin/lein
RUN chmod 744 /usr/local/bin/lein

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_22.11.1-1-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p /opt/conda
RUN /opt/conda/bin/conda install -c conda-forge -c bioconda -y python=3.9 r-base=4.1.2 r-renv blas lapack
ADD ./resources/requirements.txt /data/requirements.txt
ADD ./bin/quartet-protqc-report /opt/conda/bin/quartet-protqc-report
RUN /opt/conda/bin/pip install -r /data/requirements.txt

ADD ./resources/bin/protqc.sh /opt/conda/bin/protqc.sh
ADD ./resources/renv /opt/conda/renv
ADD ./resources/renv.lock /opt/conda/renv.lock
ADD ./build/Rprofile /opt/conda/etc/Rprofile
RUN Rscript /opt/conda/etc/Rprofile

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

# Run it
ENTRYPOINT ["quartet-protqc-report"]