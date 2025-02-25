#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --job-name=error_summary
#SBATCH --output=/scratch/lorthiois/logs/%A.o
#SBATCH --error=/scratch/lorthiois/logs/%A.e

# Path configuration
LOG_DIR="/scratch/lorthiois/logs"
OUTPUT_FILE="/work/gr-fe/lorthiois/project2/error_summary.html"
JOB_MAPPING_FILE="/work/gr-fe/lorthiois/project2/logs/model_job_mapping.txt"

# Create a temporary working directory
TEMP_DIR=$(mktemp -d)
ERROR_PATTERNS_FILE="$TEMP_DIR/error_patterns.txt"

# Common error patterns to search for {test}
cat > "$ERROR_PATTERNS_FILE" << EOF
Error:
ERROR:
Exception:
Failed:
fatal:
Fatal:
FATAL:
segmentation fault
memory allocation
not found
cannot open
undefined reference
invalid argument
time limit
OOM
panic:
Panic:
CANCELLED
EOF

# Check if mapping file exists
if [ ! -f "$JOB_MAPPING_FILE" ]; then
    echo "Job mapping file not found: $JOB_MAPPING_FILE"
    exit 1
fi

# Generate HTML header
cat > "$OUTPUT_FILE" << EOF
<!DOCTYPE html>
<html>
<head>
    <title>Deconvolution Methods Error Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #333; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { padding: 8px; text-align: left; border: 1px solid #ddd; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .error { color: #e74c3c; }
        .success { color: #2ecc71; }
        .warning { color: #f39c12; }
    </style>
</head>
<body>
    <h1>Deconvolution Methods Error Summary</h1>
    <p>Generated on $(date)</p>
    <table>
        <tr>
            <th>Method</th>
            <th>Job ID</th>
            <th>Status</th>
            <th>Error Type</th>
            <th>Error Details</th>
        </tr>
EOF

# Process each method and its job ID
while IFS=: read -r METHOD JOB_ID; do
    echo "Processing $METHOD (Job $JOB_ID)..."
    
    # Check if error file exists
    ERROR_FILE="$LOG_DIR/${JOB_ID}.e"
    if [ ! -f "$ERROR_FILE" ]; then
        # Write row with missing error file
        cat >> "$OUTPUT_FILE" << EOF
        <tr>
            <td>${METHOD}</td>
            <td>${JOB_ID}</td>
            <td class="warning">Unknown</td>
            <td>N/A</td>
            <td>Error log file not found</td>
        </tr>
EOF
        continue
    fi
    
    # Check file size - if zero, likely no errors
    if [ ! -s "$ERROR_FILE" ]; then
        cat >> "$OUTPUT_FILE" << EOF
        <tr>
            <td>${METHOD}</td>
            <td>${JOB_ID}</td>
            <td class="success">Success</td>
            <td>None</td>
            <td>No errors detected</td>
        </tr>
EOF
        continue
    fi
    
    # Search for error patterns
    ERROR_FOUND=false
    ERROR_TYPE=""
    ERROR_DETAILS=""
    
    while IFS= read -r PATTERN; do
        if grep -i "$PATTERN" "$ERROR_FILE" > /dev/null; then
            ERROR_FOUND=true
            # Capture error type from pattern
            ERROR_TYPE="$PATTERN"
            # Extract a concise error message (first occurrence)
            ERROR_LINE=$(grep -i -m 1 "$PATTERN" "$ERROR_FILE" | tr -d '\r' | sed 's/</\&lt;/g; s/>/\&gt;/g')
            # Get 1 line before and after for context
            LINE_NUM=$(grep -i -n -m 1 "$PATTERN" "$ERROR_FILE" | cut -d: -f1)
            if [ ! -z "$LINE_NUM" ]; then
                START_LINE=$((LINE_NUM > 1 ? LINE_NUM - 1 : 1))
                END_LINE=$((LINE_NUM + 1))
                ERROR_CONTEXT=$(sed -n "${START_LINE},${END_LINE}p" "$ERROR_FILE" | tr -d '\r' | sed 's/</\&lt;/g; s/>/\&gt;/g')
                ERROR_DETAILS=$(echo "$ERROR_CONTEXT" | sed 's/$/<br>/')
            else
                ERROR_DETAILS="$ERROR_LINE"
            fi
            break
        fi
    done < "$ERROR_PATTERNS_FILE"
    
    # If no specific error found but file exists, grab the last few lines
    if [ "$ERROR_FOUND" = false ] && [ -s "$ERROR_FILE" ]; then
        ERROR_TYPE="Unknown"
        ERROR_DETAILS=$(tail -5 "$ERROR_FILE" | sed 's/</\&lt;/g; s/>/\&gt;/g; s/$/<br>/')
    fi
    
    # Write row to HTML
    STATUS_CLASS="error"
    STATUS_TEXT="Failed"
    
    # Check if there's a success message
    if grep -i "success\|completed successfully" "$ERROR_FILE" > /dev/null; then
        STATUS_CLASS="success"
        STATUS_TEXT="Success with warnings"
    fi
    
    cat >> "$OUTPUT_FILE" << EOF
        <tr>
            <td>${METHOD}</td>
            <td>${JOB_ID}</td>
            <td class="${STATUS_CLASS}">${STATUS_TEXT}</td>
            <td>${ERROR_TYPE}</td>
            <td>${ERROR_DETAILS}</td>
        </tr>
EOF

done < "$JOB_MAPPING_FILE"

# Close HTML
cat >> "$OUTPUT_FILE" << EOF
    </table>
    
    <h2>Error Pattern Summary</h2>
    <table>
        <tr>
            <th>Error Pattern</th>
            <th>Count</th>
            <th>Affected Methods</th>
        </tr>
EOF

# Add error pattern summary
while IFS= read -r PATTERN; do
    # Skip empty lines
    [ -z "$PATTERN" ] && continue
    
    # Count occurrences of this pattern
    METHODS_WITH_ERROR=""
    ERROR_COUNT=0
    
    while IFS=: read -r METHOD JOB_ID; do
        ERROR_FILE="$LOG_DIR/${JOB_ID}.e"
        if [ -f "$ERROR_FILE" ] && grep -i "$PATTERN" "$ERROR_FILE" > /dev/null; then
            METHODS_WITH_ERROR="$METHODS_WITH_ERROR$METHOD, "
            ERROR_COUNT=$((ERROR_COUNT + 1))
        fi
    done < "$JOB_MAPPING_FILE"
    
    # Only add row if errors found
    if [ $ERROR_COUNT -gt 0 ]; then
        # Remove trailing comma and space
        METHODS_WITH_ERROR=${METHODS_WITH_ERROR%, }
        
        cat >> "$OUTPUT_FILE" << EOF
        <tr>
            <td>${PATTERN}</td>
            <td>${ERROR_COUNT}</td>
            <td>${METHODS_WITH_ERROR}</td>
        </tr>
EOF
    fi
done < "$ERROR_PATTERNS_FILE"

# Finish HTML
cat >> "$OUTPUT_FILE" << EOF
    </table>
</body>
</html>
EOF

echo "Error summary generated at: $OUTPUT_FILE"

# Clean up
rm -rf "$TEMP_DIR"