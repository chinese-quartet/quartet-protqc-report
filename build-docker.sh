#!/bin/bash

set -e

VERSION=$(git describe --tags `git rev-list --tags --max-count=1` --always)

# dynamically pull more interesting stuff from latest git commit
HASH=$(git show-ref --head --hash=8 head)  # first 8 letters of hash should be enough; that's what GitHub uses

# Change the version in the project.clj and resources/tservice-plugin.yaml
TRIMMED_VERSION=$(echo $VERSION | sed 's/^v//')
# If running on macOS, use sed -i '' instead of sed -i
if [[ "$OSTYPE" == "darwin"* ]]; then
  sed -i "" "s/(defproject quartet-protqc-report \".*\"/(defproject quartet-protqc-report \"${TRIMMED_VERSION}\"/g" project.clj
  sed -i "" "s/version: v.*$/version: v${TRIMMED_VERSION}-${HASH}/g" resources/tservice-plugin.yaml
else
  sed -i "s/(defproject quartet-protqc-report \".*\"/(defproject quartet-protqc-report \"${TRIMMED_VERSION}\"/g" project.clj
  sed -i "s/version: v.*$/version: v${TRIMMED_VERSION}-${HASH}/g" resources/tservice-plugin.yaml
fi

# Build base docker image
docker build -t quartet-protqc-report:${VERSION}-${HASH} .

if [ "$1" == "--push" ]; then
    docker tag quartet-protqc-report:${VERSION}-${HASH} ghcr.io/chinese-quartet/quartet-protqc-report:${VERSION}-${HASH} && \
    docker push ghcr.io/chinese-quartet/quartet-protqc-report:${VERSION}-${HASH}
fi
