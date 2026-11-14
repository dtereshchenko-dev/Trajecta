#!/bin/sh
set -eu

seed_dir=${1:?Missing packaged data directory}
data_root=${XDG_DATA_HOME:-"$HOME/.local/share"}/Trajecta/Trajecta
mkdir -p "$data_root/observatories" "$data_root/ephemerides" "$data_root/earth_orientation"

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
for name in earth_historic.bpc earth_latest.bpc earth_predict.bpc eopc04.62-now eopc04.cnc ET-UT.dat; do
    install_if_missing "$seed_dir/$name" "$data_root/earth_orientation/$name"
done
rm -f "$seed_dir/OBSCODE.dat" "$seed_dir/RADCODE.dat" "$seed_dir/linux_p1550p2650.440" \
    "$seed_dir/earth_historic.bpc" "$seed_dir/earth_latest.bpc" "$seed_dir/earth_predict.bpc" \
    "$seed_dir/eopc04.62-now" "$seed_dir/eopc04.cnc" "$seed_dir/ET-UT.dat"
rmdir "$seed_dir"
