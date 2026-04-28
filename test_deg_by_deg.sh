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
MACRO_DIR="${SCRIPT_DIR}/macros"

DO_SIM=true
DO_HIST=true

usage() {
    echo "Usage: $0 [--hist-only|--sim-only]"
    echo "  --hist-only  Only run txtToHist post-processing"
    echo "  --sim-only   Only run Geant simulation loop"
}

if [[ $# -gt 1 ]]; then
    usage
    exit 1
fi

if [[ $# -eq 1 ]]; then
    case "$1" in
        --hist-only)
            DO_SIM=false
            DO_HIST=true
            ;;
        --sim-only)
            DO_SIM=true
            DO_HIST=false
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage
            exit 1
            ;;
    esac
fi

cd "$BUILD_DIR"

if [[ "$DO_SIM" == true ]]; then
    cmake ..
    make -j2
fi

if [[ "$DO_SIM" == true && ! -x "$RUNNER" ]]; then
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

if [[ "$DO_SIM" == true ]]; then
    for deg in $(seq 0 14); do
        next_deg=$((deg + 1))
        xsec_file="testdhkm_${deg}-${next_deg}_tmp.txt"
        xsec_path="${BUILD_DIR}/${xsec_file}"
        tmp_files+=("$xsec_path")

        # Write two-line cross section in radians with equal probabilities (0.5, 0.5)
        LC_ALL=C awk -v d0="$deg" -v d1="$next_deg" 'BEGIN {
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
fi

# Run txtToHist on the generated output files

# Norma names outputs from NUMERIC_MIE_FPATH directly:
# output<NUMERIC_MIE_FPATH><radius>_<g>_<p>.txt
# We therefore locate files by prefix rather than assuming output_<deg>-<deg>.txt.
if [[ "$DO_HIST" == true ]]; then
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
            input_basename="$(basename "$latest_match")"
            plot_prefix="output_${deg}-${next_deg}"
            input_from_macros="../build/${input_basename}"

            echo "Processing ${input_from_macros} with prefix ${plot_prefix}"
            (
                cd "$MACRO_DIR"
                root -x -l -b -q "txtToHist.C(\"${input_from_macros}\", \"${plot_prefix}\")"
            )
        else
            echo "Warning: No Norma output matched ${output_prefix}*.txt in ${BUILD_DIR}. Skipping ${deg}-${next_deg}."
        fi
    done
fi

# Plot beamprofiles
if [[ "$DO_HIST" == true ]]; then
    cd "$MACRO_DIR"
    root -b -q plot_det3_beamprofiler.C
    echo "Beam profile plots generated."
fi


cd "$SCRIPT_DIR"