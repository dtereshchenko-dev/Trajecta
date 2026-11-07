#!/usr/bin/env bash
set -euo pipefail

setup=$(realpath "${1:?Usage: install-smoke.sh SETUP_FILE [STREWN_TEST_BINARY]}")
strewn_tests=${2:-}
if [[ -n $strewn_tests ]]; then
    strewn_tests=$(realpath "$strewn_tests")
fi
tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
export HOME="$tmp/home" XDG_DATA_HOME="$tmp/data" XDG_CONFIG_HOME="$tmp/config" QT_QPA_PLATFORM=offscreen
mkdir -p "$HOME/.local/share/applications"
target="$tmp/installed"
"$setup" --root "$target" --default-answer --confirm-command --accept-licenses install
test -s "$XDG_DATA_HOME/Trajecta/Trajecta/ephemerides/linux_p1550p2650.440"
test ! -e "$target/usr/share/trajecta-seed"
test ! -e "$target/usr/bin/linux_p1550p2650.440"
test -f "$target/usr/plugins/wayland-graphics-integration-client/libqt-plugin-wayland-egl.so"
if find "$target/examples" -type d -perm -0002 -print -quit | grep -q .; then
    echo 'Example directory is writable by other users.' >&2
    exit 1
fi
for style in Fusion Imagine Material Universal; do
    test ! -e "$target/usr/qml/QtQuick/Controls/$style"
done
"$target/usr/bin/trajecta-cli" "$target/examples/NEO/impacts/2008_TC3.txt" --quiet > "$tmp/fit.log"
cat "$tmp/fit.log"
grep -q 'converged=true' "$tmp/fit.log"
if [[ -n $strewn_tests ]]; then
    env LD_LIBRARY_PATH="$target/usr/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
        GDAL_DATA="$target/usr/share/gdal" PROJ_DATA="$target/usr/share/proj" \
        "$strewn_tests" --gtest_filter='WeatherTest.RequiredRasterDriversAreAvailable:WeatherTest.RuntimeRasterResourcesAreAvailable:WeatherTest.RasterSamplingPreservesAxesAndPixelCenterConvention:JobTest.MapsAreLibraryOwnedAndReplayable'
fi
test -d "$XDG_DATA_HOME/Trajecta/work"
test ! -e "$target/usr/bin/work"
test -s "$XDG_DATA_HOME/Trajecta/Trajecta/ephemerides/linux_p1550p2650.440"
test -s "$XDG_DATA_HOME/Trajecta/Trajecta/observatories/OBSCODE.dat"
set +e
timeout 20s xvfb-run -a env QT_QPA_PLATFORM=xcb QT_QUICK_BACKEND=software "$target/usr/bin/trajecta" > "$tmp/ui.log" 2>&1
status=$?
set -e
cat "$tmp/ui.log"
if [[ $status -ne 124 ]] || grep -Eq 'Error decoding:|QQmlApplicationEngine failed to load|module .* is not installed' "$tmp/ui.log"; then
    echo "UI smoke test failed (exit: $status)" >&2
    exit 1
fi
"$target/TrajectaMaintenance" --default-answer --confirm-command purge
test -s "$XDG_DATA_HOME/Trajecta/Trajecta/ephemerides/linux_p1550p2650.440"
test -d "$XDG_DATA_HOME/Trajecta/work"
echo 'Linux installer, CLI fit, UI startup and data preservation passed.'
