#!/usr/bin/env bash
set -euo pipefail

# Runs 15 simulations for angle bins: 0-1, 1-2, ..., 14-15 degrees.
# For each bin:
#   - creates a temporary 2-line differential cross section file (0.5 / 0.5 weights),
#   - runs offset_runner.sh with zero cell offset and 10M events,
#   - deletes the temporary cross section file.


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
RUNNER="${BUILD_DIR}/offset_runner.sh"

cd "$BUILD_DIR"
cmake ..
make -j2

if [[ ! -x "$RUNNER" ]]; then
    echo "Error: $RUNNER not found or not executable."
    exit 1
fi

tmp_files=()
cleanup() {
    for f in "${tmp_files[@]:-}"; do
        [[ -f "$f" ]] && rm -f "$f"
    done
}
trap cleanup EXIT

for deg in $(seq 0 14); do
    next_deg=$((deg + 1))
    xsec_file="testdhkm_${deg}-${next_deg}_tmp.txt"
    xsec_path="${BUILD_DIR}/${xsec_file}"
    tmp_files+=("$xsec_path")

    # Write two-line cross section in radians with equal probabilities (0.5, 0.5)
    awk -v d0="$deg" -v d1="$next_deg" 'BEGIN {
        pi = atan2(0, -1)
        printf("%.18f 0.5\n", d0*pi/180.0)
        printf("%.18f 0.5\n", d1*pi/180.0)
    }' > "$xsec_path"

    echo "Running bin ${deg}-${next_deg} deg using $(basename "$xsec_file")"

    # start_x end_x step_x start_y end_y step_y z events xsection_file
    # Zero offset -> generated macro will use: /Norma/DetectorConstruction/Cell/setOffset 0 0 0 um
    # offset_runner.sh currently uses CELLSIZE=10.0 (10 um) internally.
    "$RUNNER" 0 0 1 0 0 1 0 10000000 "$xsec_file"

    rm -f "$xsec_path"
done

echo "All 15 degree-by-degree simulations completed."

# Run txtToHist on the generated output files
OUTPUT_DIR="${BUILD_DIR}/output"
mkdir -p "$OUTPUT_DIR"

# Norma names outputs from NUMERIC_MIE_FPATH directly:
# output<NUMERIC_MIE_FPATH><radius>_<g>_<p>.txt
# We therefore locate files by prefix rather than assuming output_<deg>-<deg>.txt.
for deg in $(seq 0 14); do
    next_deg=$((deg + 1))
    xsec_file="testdhkm_${deg}-${next_deg}_tmp.txt"
    output_prefix="output${xsec_file}"
    latest_match=""

    while IFS= read -r match; do
        latest_match="$match"
        break
    done < <(find "$BUILD_DIR" -maxdepth 1 -type f -name "${output_prefix}*.txt" -printf "%T@ %p\n" | sort -nr | awk '{print $2}')

    if [[ -n "$latest_match" && -f "$latest_match" ]]; then
        root_prefix="${OUTPUT_DIR}/output_${deg}-${next_deg}"
        echo "Processing ${latest_match} -> ${root_prefix}_0_alpha_output.root"
        root -x -l -b -q "$SCRIPT_DIR/macros/txtToHist.C(\"${latest_match}\", \"${root_prefix}\")"
    else
        echo "Warning: No Norma output matched ${output_prefix}*.txt in ${BUILD_DIR}. Skipping ${deg}-${next_deg}."
    fi
done

cd "$SCRIPT_DIR"