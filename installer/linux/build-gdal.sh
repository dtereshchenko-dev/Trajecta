#!/usr/bin/env bash
set -euo pipefail

prefix=$(realpath -m "${1:?Usage: build-gdal.sh INSTALL_PREFIX}")
work=$(mktemp -d)
trap 'rm -rf -- "$work"' EXIT

version=3.12.2
archive="$work/gdal.tar.gz"
curl --fail --location --retry 3 -o "$archive" \
    "https://github.com/OSGeo/gdal/releases/download/v${version}/gdal-${version}.tar.gz"
printf '%s  %s\n' 458a899feea38000258144517fedc6662ebba255971669d2901ba77e9e8fbf79 "$archive" | sha256sum -c -
tar -xzf "$archive" -C "$work"
cmake -S "$work/gdal-$version" -B "$work/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" \
    -DBUILD_APPS=OFF -DBUILD_TESTING=OFF \
    -DBUILD_PYTHON_BINDINGS=OFF -DBUILD_JAVA_BINDINGS=OFF -DBUILD_CSHARP_BINDINGS=OFF \
    -DGDAL_BUILD_OPTIONAL_DRIVERS=OFF -DOGR_BUILD_OPTIONAL_DRIVERS=OFF \
    -DGDAL_USE_EXTERNAL_LIBS=OFF -DGDAL_USE_OPENJPEG=ON \
    -DGDAL_USE_PNG=ON -DGDAL_USE_LIBAEC=ON -DGDAL_USE_JPEG=ON \
    -DGDAL_ENABLE_DRIVER_GRIB=ON -DGDAL_ENABLE_DRIVER_JP2OPENJPEG=ON \
    -DOGR_ENABLE_DRIVER_CSV=ON
cmake --build "$work/build" --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-4}"
cmake --install "$work/build"
