#!/bin/sh
set -eu

app_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
export LD_LIBRARY_PATH="$app_dir/../lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export QT_PLUGIN_PATH="$app_dir/../plugins"
export QML2_IMPORT_PATH="$app_dir/../qml"
export GDAL_DATA="$app_dir/../share/gdal"
export PROJ_DATA="$app_dir/../share/proj"
case "$(basename "$0")" in
    trajecta-cli) executable=TrajectaCLI ;;
    strewn-cli) executable=StrewnCLI ;;
    *) executable=Trajecta ;;
esac
exec "$app_dir/$executable" "$@"
