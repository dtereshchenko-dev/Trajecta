#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/../.." && pwd)
build=$(realpath "${1:?Usage: deploy.sh BUILD IFW_DIR LINUXDEPLOY OUTPUT [VERSION] [GDAL_DATA_DIR]}")
ifw=$(realpath "${2:?Missing Qt Installer Framework directory}")
linuxdeploy=$(realpath "${3:?Missing linuxdeploy executable}")
output=$(realpath -m "${4:?Missing output directory}")
version=${5:-$(tr -d '\r\n' < "$build/TrajectaUI/trajecta-version.txt")}
gdal_data=${6:-/usr/share/gdal}
[[ $version =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || { echo 'Invalid version' >&2; exit 1; }

stage=$(mktemp -d)
trap 'rm -rf -- "$stage"' EXIT
app="$stage/packages/com.trajecta.app/data"
meta="$stage/packages/com.trajecta.app/meta"
mkdir -p "$app/usr/bin" "$meta" "$stage/config" "$output"
sed "s/@TRAJECTA_VERSION@/$version/g" "$root/installer/linux/config.xml" > "$stage/config/config.xml"
sed "s/@TRAJECTA_VERSION@/$version/g" "$root/installer/linux/package.xml" > "$meta/package.xml"
cat "$root/installer/linux/controller.qs" "$root/installer/common/lifecycle.qs" > "$stage/config/controller.qs"
cp "$root/installer/linux/installscript.qs" "$meta/installscript.qs"
cp "$root/installer/config/icon.png" "$root/installer/config/banner.png" "$stage/config/"
cp "$root/installer/packages/com.trajecta.app/meta/license.txt" "$meta/license.txt"

export APPIMAGE_EXTRACT_AND_RUN=1
export QMAKE=/usr/bin/qmake6
export QML_SOURCES_PATHS="$root/TrajectaUI/qml"
export EXTRA_PLATFORM_PLUGINS='libqwayland.so'
wayland_egl="$(qmake6 -query QT_INSTALL_PLUGINS)/wayland-graphics-integration-client/libqt-plugin-wayland-egl.so"
[[ -f $wayland_egl ]] || { echo 'Qt Wayland EGL plugin not found' >&2; exit 1; }
echo 'Deploying application and dependencies'
if ! "$linuxdeploy" --verbosity 2 --appdir "$app" --plugin qt \
    --executable "$build/TrajectaUI/Trajecta" --executable "$build/TrajectaCLI/TrajectaCLI" \
    --executable "$build/StrewnField/StrewnCLI" \
    --library "$build/OrbitalCore/libOrbitalCore.so" --library "$build/StrewnField/libStrewnCore.so" \
    --library "$wayland_egl" \
    --desktop-file "$root/installer/linux/trajecta.desktop" --icon-filename trajecta \
    --icon-file "$root/TrajectaUI/resources/icons/workspace/app-logo.svg" > "$stage/linuxdeploy.log" 2>&1; then
    tail -n 100 "$stage/linuxdeploy.log" >&2
    exit 1
fi

mkdir -p "$app/usr/plugins/wayland-graphics-integration-client"
cp "$wayland_egl" "$app/usr/plugins/wayland-graphics-integration-client/"
rm -f -- "$app/usr/lib/libqt-plugin-wayland-egl.so"
for style in Fusion Imagine Material Universal; do
    rm -rf -- "$app/usr/qml/QtQuick/Controls/$style"
    rm -f -- "$app/usr/lib/libQt6QuickControls2${style}"*.so*
done
rm -rf -- "$app/usr/qml/QtQuick/Controls/designer" "$app/usr/plugins/qmltooling"

cp "$root/installer/linux/launch.sh" "$app/usr/bin/trajecta"
cp "$root/installer/linux/launch.sh" "$app/usr/bin/trajecta-cli"
cp "$root/installer/linux/launch.sh" "$app/usr/bin/strewn-cli"
cp "$root/installer/linux/install-data.sh" "$app/usr/bin/install-data.sh"
chmod 755 "$app/usr/bin/trajecta" "$app/usr/bin/trajecta-cli" "$app/usr/bin/strewn-cli"
mkdir -p "$app/usr/share/trajecta-seed"
cp "$root/data/linux_p1550p2650.440" "$root/data/OBSCODE.dat" "$root/data/RADCODE.dat" "$app/usr/share/trajecta-seed/"
mkdir -p "$app/usr/bin/debias" "$app/usr/bin/weights" "$app/usr/bin/earth_orientation" "$app/usr/share"
cp "$root/data/debias/efcc18.dat" "$root/data/debias/fcct14.dat" "$app/usr/bin/debias/"
cp "$root/data/weights/vfcc17.rules" "$app/usr/bin/weights/"
cp "$root/data/earth_orientation/eopc04.62-now" "$root/data/earth_orientation/eopc04.cnc" \
    "$root/data/earth_orientation/ET-UT.dat" "$root/data/earth_orientation/finals2000A.all" "$app/usr/bin/earth_orientation/"
cp -a "$gdal_data" "$app/usr/share/gdal"
cp -a /usr/share/proj "$app/usr/share/"
cp -a "$root/data/examples" "$app/examples"
find "$app/examples" -type d -exec chmod 755 {} +
find "$app/examples" -type f -exec chmod 644 {} +
mkdir -p "$app/licenses"
cp "$root/LICENSE" "$root/licenses/LGPL-3.0.txt" "$app/licenses/"
cp "$root/TrajectaUI/THIRD_PARTY.md" "$app/"
cp "$root/TrajectaUI/resources/fonts/LICENSE-SourceSans3.txt" "$app/licenses/"
cp "$root/TrajectaUI/resources/fonts/LICENSE-CascadiaMono.txt" "$app/licenses/"
cp "$root/OrbitalCore/third_party/erfa/LICENSE" "$app/licenses/ERFA.txt"
cp "$root/TrajectaUI/resources/icons/workspace/LICENSE" "$app/licenses/workspace-icons.txt"
for name in leaflet plotly lucide; do
    cp "$root/StrewnField/output/web/vendor/$name/LICENSE" "$app/licenses/$name.txt"
done

"$ifw/bin/repogen" --compression 1 -p "$stage/packages" "$output/repository"
"$ifw/bin/binarycreator" --hybrid --compression 1 -c "$stage/config/config.xml" \
    -p "$stage/packages" "$output/TrajectaSetup-x86_64.run"
