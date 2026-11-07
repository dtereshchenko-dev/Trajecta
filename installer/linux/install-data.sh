#!/bin/sh
set -eu

seed_dir=${1:?Missing packaged data directory}
data_root=${XDG_DATA_HOME:-"$HOME/.local/share"}/Trajecta/Trajecta
mkdir -p "$data_root/observatories" "$data_root/ephemerides"

install_if_missing() {
    source=$1
    destination=$2
    [ -s "$source" ] || { echo "Missing packaged data: $source" >&2; return 1; }
    if [ -e "$destination" ]; then return; fi
    temporary=$(mktemp "${destination}.XXXXXX")
    if ! cp "$source" "$temporary"; then
        rm -f "$temporary"
        return 1
    fi
    if ! ln "$temporary" "$destination" 2>/dev/null && [ ! -e "$destination" ]; then
        rm -f "$temporary"
        return 1
    fi
    rm -f "$temporary"
}

for name in OBSCODE.dat RADCODE.dat; do
    install_if_missing "$seed_dir/$name" "$data_root/observatories/$name"
done
install_if_missing "$seed_dir/linux_p1550p2650.440" "$data_root/ephemerides/linux_p1550p2650.440"
rm -f "$seed_dir/OBSCODE.dat" "$seed_dir/RADCODE.dat" "$seed_dir/linux_p1550p2650.440"
rmdir "$seed_dir"
