#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/../.." && pwd)
ifw=$(realpath "${1:?Usage: update-smoke.sh IFW_DIR}")
tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
export HOME="$tmp/home" XDG_DATA_HOME="$tmp/data" XDG_CONFIG_HOME="$tmp/config" QT_QPA_PLATFORM=offscreen
mkdir -p "$HOME/.local/share/applications" "$XDG_DATA_HOME/Trajecta/work" "$XDG_CONFIG_HOME/Trajecta" \
    "$tmp/packages/com.trajecta.app/data" "$tmp/packages/com.trajecta.app/meta"
meta="$tmp/packages/com.trajecta.app/meta"
data="$tmp/packages/com.trajecta.app/data"
target="$tmp/installed with spaces"
cp "$root/installer/linux/installscript.qs" "$meta/installscript.qs"
cat "$root/installer/linux/controller.qs" "$root/installer/common/lifecycle.qs" > "$tmp/controller.qs"
cp "$root/installer/packages/com.trajecta.app/meta/license.txt" "$meta/license.txt"
cat > "$tmp/config.xml" <<EOF
<Installer><Name>Trajecta smoke</Name><Version>1.0.1</Version><Title>Trajecta smoke</Title>
<Publisher>Test</Publisher><MaintenanceToolName>TrajectaMaintenance</MaintenanceToolName>
<ControlScript>controller.qs</ControlScript></Installer>
EOF
write_version() {
    sed -e "s/@TRAJECTA_VERSION@/$1/g" -e 's/@TRAJECTA_RELEASE_DATE@/2026-01-01/g' \
        "$root/installer/linux/package.xml" > "$meta/package.xml"
    if [[ $1 == 1.0.1 ]]; then
        sed -i 's/<Essential>false<\/Essential>/<Essential>true<\/Essential>/' "$meta/package.xml"
    fi
    printf '%s' "$1" > "$data/version.txt"
    mkdir -p "$data/usr/bin" "$data/usr/share/trajecta-seed"
    printf '%s\n' '#include <unistd.h>' 'int main(void) { sleep(120); return 0; }' > "$tmp/process.c"
    cc "$tmp/process.c" -o "$data/usr/bin/Trajecta"
    cp "$root/installer/linux/install-data.sh" "$data/usr/bin/install-data.sh"
    for name in OBSCODE.dat RADCODE.dat linux_p1550p2650.440; do
        printf '%s' "packaged $name $1" > "$data/usr/share/trajecta-seed/$name"
    done
}
write_version 1.0.1
"$ifw/bin/binarycreator" --hybrid -c "$tmp/config.xml" -p "$tmp/packages" "$tmp/setup"
"$tmp/setup" --root "$target" --default-answer --confirm-command --accept-licenses install
test "$(cat "$target/version.txt")" = 1.0.1
test "$(cat "$XDG_DATA_HOME/Trajecta/Trajecta/ephemerides/linux_p1550p2650.440")" = 'packaged linux_p1550p2650.440 1.0.1'
test ! -e "$target/usr/share/trajecta-seed"
grep -Fq "Exec=\"$target/usr/bin/trajecta\"" "$HOME/.local/share/applications/trajecta.desktop"
printf '%s' 'user work' > "$XDG_DATA_HOME/Trajecta/work/keep.txt"
printf '%s' 'user settings' > "$XDG_CONFIG_HOME/Trajecta/Trajecta.conf"
printf '%s' 'updated stations' > "$XDG_DATA_HOME/Trajecta/Trajecta/observatories/OBSCODE.dat"
set +e
"$tmp/setup" --root "$target" --auto-answer OverwriteTargetDirectory=Yes --default-answer \
    --confirm-command --accept-licenses install > "$tmp/reinstall.log" 2>&1
status=$?
set -e
if [[ $status -eq 0 ]] || ! grep -Eq 'Trajecta is already installed|TargetDirectoryInUse' "$tmp/reinstall.log"; then
    cat "$tmp/reinstall.log" >&2
    echo 'Existing installation was not protected.' >&2
    exit 1
fi
test "$(cat "$target/version.txt")" = 1.0.1
test ! -e "$target/TrajectaMaintenance.new"
write_version 1.0.2
"$ifw/bin/repogen" -p "$tmp/packages" "$tmp/repository"
"$target/TrajectaMaintenance" --set-temp-repository "file://$tmp/repository" check-updates > "$tmp/updates.xml"
grep -q 'version="1.0.2"' "$tmp/updates.xml"
"$target/usr/bin/Trajecta" 120 &
app_pid=$!
trap 'kill "$app_pid" 2>/dev/null || true; rm -rf -- "$tmp"' EXIT
sleep 0.1
kill -0 "$app_pid"
set +e
"$target/TrajectaMaintenance" --set-temp-repository "file://$tmp/repository" --default-answer --confirm-command update > "$tmp/running-update.log" 2>&1
status=$?
set -e
if [[ $status -eq 0 || $status -eq 6 ]] || ! kill -0 "$app_pid" 2>/dev/null; then
    cat "$tmp/running-update.log" >&2
    echo 'Update did not protect the running application.' >&2
    exit 1
fi
test "$(cat "$target/version.txt")" = 1.0.1
test "$(cat "$XDG_DATA_HOME/Trajecta/work/keep.txt")" = 'user work'
set +e
"$target/TrajectaMaintenance" --default-answer --confirm-command purge RemoveUserData=1 > "$tmp/running-purge.log" 2>&1
status=$?
set -e
if [[ $status -ne 0 ]] || ! kill -0 "$app_pid" 2>/dev/null || ! grep -q 'Close Trajecta before removing it' "$tmp/running-purge.log"; then
    cat "$tmp/running-purge.log" >&2
    echo 'Unattended removal did not warn about the running application.' >&2
    exit 1
fi
test "$(cat "$XDG_DATA_HOME/Trajecta/work/keep.txt")" = 'user work'
kill "$app_pid"
wait "$app_pid" 2>/dev/null || true
trap 'rm -rf -- "$tmp"' EXIT
"$tmp/setup" --root "$target" --default-answer --confirm-command --accept-licenses install
test "$(cat "$target/version.txt")" = 1.0.1
test "$(cat "$XDG_DATA_HOME/Trajecta/work/keep.txt")" = 'user work'
set +e
"$target/TrajectaMaintenance" --set-temp-repository "file://$tmp/repository" --default-answer --confirm-command update
status=$?
set -e
[[ $status -eq 0 ]]
test "$(cat "$target/version.txt")" = 1.0.2
test "$(cat "$XDG_DATA_HOME/Trajecta/work/keep.txt")" = 'user work'
test "$(cat "$XDG_DATA_HOME/Trajecta/Trajecta/observatories/OBSCODE.dat")" = 'updated stations'
test "$(cat "$XDG_DATA_HOME/Trajecta/Trajecta/ephemerides/linux_p1550p2650.440")" = 'packaged linux_p1550p2650.440 1.0.1'
test ! -e "$target/usr/share/trajecta-seed"
"$target/TrajectaMaintenance" --default-answer --confirm-command purge
test -f "$XDG_DATA_HOME/Trajecta/work/keep.txt"
test -f "$XDG_CONFIG_HOME/Trajecta/Trajecta.conf"
"$tmp/setup" --root "$target" --default-answer --confirm-command --accept-licenses install
"$target/TrajectaMaintenance" --default-answer --confirm-command purge RemoveUserData=1
test ! -e "$XDG_DATA_HOME/Trajecta"
test ! -e "$XDG_CONFIG_HOME/Trajecta/Trajecta.conf"
echo 'Linux IFW update, uninstall and user-data choice passed.'
