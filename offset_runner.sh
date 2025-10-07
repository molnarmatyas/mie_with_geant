#!/bin/bash

# Function to format numbers
format_number() 
{
    local num=$1
    if [[ "$num" =~ [eE] ]]; then
        num=$(printf "%.10f" "$num")
    fi
    # Remove trailing zeros and possible . from printf
    # num=$(echo "$num" | sed -e 's/\.0*$//' -e 's/^\./0./' -e 's/0*$//')
    echo "$num"
}

# Function to display usage information
usage() 
{
    echo "Usage: $0 <start_x> <end_x> <step_x> <start_y> <end_y> <step_y> <z_value> <num_events> <crossection_inputfile>"
    echo "  All coordinate values (start/end/z) must be floats between -100 and 100"
    echo "  Step sizes must be positive floats"
    echo "  Number of events must be positive integer"
    exit 1
}

# Validate argument count
if [ $# -ne 9 ]; then
    echo "Error: Incorrect number of arguments"
    usage
fi

# Assign arguments
start_x=$1
end_x=$2
step_x=$3
start_y=$4
end_y=$5
step_y=$6
z_value=$7
num_events=$8
inputxsection="$9"

# Function to validate floats with range check
validate_float() {
    local num=$1
    local name=$2
    local min=$3
    local max=$4
    
    # Check if it's a valid float (supports scientific notation too)
    if ! [[ "$num" =~ ^-?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ || "$num" =~ ^-?\.[0-9]+([eE][-+]?[0-9]+)?$ ]]; then
        echo "Error: $name must be a number (can be float)"
        usage
    fi
    
    # Check range
    if (( $(echo "$num < $min || $num > $max" | bc -l) )); then
        echo "Error: $name must be between $min and $max"
        usage
    fi
}

# Function to validate positive floats
validate_positive_float() {
    local num=$1
    local name=$2
    
    validate_float "$num" "$name" 0.0000001 100  # Smallest positive value
    if (( $(echo "$num <= 0" | bc -l) )); then
        echo "Error: $name must be greater than 0"
        usage
    fi
}

# Function to validate positive integers
validate_int() {
    local num=$1
    local name=$2
    
    if ! [[ "$num" =~ ^[0-9]+$ ]]; then
        echo "Error: $name must be a positive integer"
        usage
    fi
    
    if [ "$num" -le 0 ]; then
        echo "Error: $name must be greater than 0"
        usage
    fi
}

# Validate all inputs
validate_float "$start_x" "Start X" -100 100
validate_float "$end_x" "End X" -100 100
validate_positive_float "$step_x" "X step size"
validate_float "$start_y" "Start Y" -100 100
validate_float "$end_y" "End Y" -100 100
validate_positive_float "$step_y" "Y step size"
validate_float "$z_value" "Z value" -100 100
validate_int "$num_events" "Number of events"

# Check directionality
if (( $(echo "($end_x - $start_x) * $step_x < 0" | bc -l) )); then
    echo "Error: X step direction contradicts start/end points"
    exit 1
fi

if (( $(echo "($end_y - $start_y) * $step_y < 0" | bc -l) )); then
    echo "Error: Y step direction contradicts start/end points"
    exit 1
fi

# Calculate number of steps in each direction
x_steps=$(echo "scale=10; ($end_x - $start_x)/$step_x" | bc)
y_steps=$(echo "scale=10; ($end_y - $start_y)/$step_y" | bc)

# Round up to nearest integer
x_steps=$(echo "scale=0; ($x_steps + 0.999)/1" | bc)
y_steps=$(echo "scale=0; ($y_steps + 0.999)/1" | bc)

# Original macro file template
template_macro="run2_template.mac"
output_macro="generated_macro.mac"

if [ ! -f "$template_macro" ]; then
    echo "Error: Template macro file '$template_macro' not found"
    exit 1
fi

total_iterations=$((x_steps * y_steps))
current_iteration=0

echo "X range: $start_x to $end_x in steps of $step_x ($x_steps steps)"
echo "Y range: $start_y to $end_y in steps of $step_y ($y_steps steps)"
echo "Z value: $z_value"
echo "Total iterations: $total_iterations"

# This is where the fun begins!!!
for ((i=0; i<=x_steps; i++)); do
    current_x=$(echo "scale=10; $start_x + $i * $step_x" | bc)
    current_x=$(format_number "$current_x")
    
    # Ensure we don't exceed end_x due to floating point rounding
    if (( $(echo "$current_x > $end_x" | bc -l) )); then
        current_x=$end_x
    fi
    
    for ((j=0; j<=y_steps; j++)); do
        current_y=$(echo "scale=10; $start_y + $j * $step_y" | bc)
        current_y=$(format_number "$current_y")
        
        # Ensure we don't exceed end_y
        if (( $(echo "$current_y > $end_y" | bc -l) )); then
            current_y=$end_y
        fi
        
        current_iteration=$((current_iteration + 1))
        
        # Create the modified macro file
        {
            grep -v "/Norma/DetectorConstruction/Cell/setOffset" "$template_macro" | \
            grep -v "/run/beamOn"
            echo "/Norma/DetectorConstruction/Cell/setOffset $current_x $current_y $z_value um"
            echo "/run/beamOn $num_events"
        } > "$output_macro"
        
        printf "Running %d/%d: X=%-10s Y=%-10s Z=%-10s\n" \
               "$current_iteration" "$total_iterations" "$current_x" "$current_y" "$z_value"
        
        NUMERIC_MIE_FPATH=${inputxsection} CELL_RADIUS_UM=7.5 ./Norma -b 6.5 -m "$output_macro" > /dev/null # FIXME use correct cell size!!!
        
        #if [ $? -ne 0 ]; then
        #    echo "Error: Command failed for X=$current_x, Y=$current_y"
        #    exit 1
        #fi
    done
done

echo "Successfully completed $total_iterations iterations"
