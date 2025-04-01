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
DEFAULT_DATASET_PREFIX="TB_D100-bulk_null"
DEFAULT_SAMPLE_FILTER="AB"  # Default sample filter

# Parse command line arguments
DATASET_PREFIX="${1:-$DEFAULT_DATASET_PREFIX}"  # Use default prefix if not provided
SAMPLE_FILTER="${2:-$DEFAULT_SAMPLE_FILTER}"    # Use default filter if not provided

# New format for log directory includes both prefix and filter
GLOBAL_LOG_DIR="/work/gr-fe/lorthiois/DeconBenchmark/logs/${DATASET_PREFIX}_${SAMPLE_FILTER}"
OUTPUT_FILE="${GLOBAL_LOG_DIR}/error_summary_${DATASET_PREFIX}_${SAMPLE_FILTER}.html"

# Use the new path format for mapping files
JOB_MAPPING_FILE="${GLOBAL_LOG_DIR}/model_job_mapping.txt"
PAIRED_MAPPING_FILE="${GLOBAL_LOG_DIR}/paired_plot_job_mapping.txt"
STATS_MAPPING_FILE="${GLOBAL_LOG_DIR}/stats_job_mapping.txt"

# Create directory if it doesn't exist
mkdir -p "$GLOBAL_LOG_DIR"

# If stats mapping file doesn't exist, create it by looking for benchmark stats jobs
if [ ! -f "$STATS_MAPPING_FILE" ] && [ -d "$GLOBAL_LOG_DIR" ]; then
    # Look for stats job logs that match pattern with both dataset prefix and sample filter
    stats_pattern="benchmarking_${DATASET_PREFIX}_${SAMPLE_FILTER}|benchmark_stats.*${DATASET_PREFIX}.*${SAMPLE_FILTER}|paired_benchmark_stats.sh"
    if grep -l -E "$stats_pattern" "${LOG_DIR}/"*.o 2>/dev/null | head -1 > /dev/null; then
        echo "# Auto-generated stats job mapping for ${DATASET_PREFIX}_${SAMPLE_FILTER}" > "$STATS_MAPPING_FILE"
        for log_file in $(grep -l -E "$stats_pattern" "${LOG_DIR}/"*.o 2>/dev/null); do
            job_id=$(basename "$log_file" .o)
            echo "stats_${DATASET_PREFIX}_${SAMPLE_FILTER}:${job_id}" >> "$STATS_MAPPING_FILE"
        done
        echo "Auto-generated stats mapping file: $STATS_MAPPING_FILE"
    fi
fi

DATA_NAME="${DATASET_PREFIX}_${SAMPLE_FILTER}"

# Create a temporary working directory
TEMP_DIR=$(mktemp -d)
ERROR_PATTERNS_FILE="$TEMP_DIR/error_patterns.txt"

# Common error patterns to search for
cat > "$ERROR_PATTERNS_FILE" << EOF
Error:
ERROR:
Exception:
Failed:
fatal:
Fatal:
FATAL:
Error
ERROR
Exception
Failed
fatal
Fatal
FATAL
segmentation fault
memory allocation
not found
cannot open
undefined reference
invalid argument
time limit
OOM
Panic:
Panic
CANCELLED
EOF

# Check if mapping files exist
if [ ! -f "$JOB_MAPPING_FILE" ]; then
    echo "Job mapping file not found: $JOB_MAPPING_FILE"
    # Try to look for a job mapping file without the sample filter (for backward compatibility)
    LEGACY_JOB_MAPPING_FILE="/work/gr-fe/lorthiois/DeconBenchmark/logs/${DATASET_PREFIX}/model_job_mapping.txt"
    if [ -f "$LEGACY_JOB_MAPPING_FILE" ]; then
        echo "Found legacy job mapping file: $LEGACY_JOB_MAPPING_FILE"
        JOB_MAPPING_FILE="$LEGACY_JOB_MAPPING_FILE"
    fi
    # Continue but note the missing file in the HTML output
fi

# Generate HTML header
cat > "$OUTPUT_FILE" << EOF
<!DOCTYPE html>
<html>
<head>
    <title>Deconvolution, Paired Plot, and Stats Error Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2 { color: #333; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { padding: 8px; text-align: left; border: 1px solid #ddd; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .error { color: #e74c3c; }
        .success { color: #2ecc71; }
        .warning { color: #f39c12; }
        .missing { color: #95a5a6; }
        .dashboard { display: flex; justify-content: space-between; margin-bottom: 20px; }
        .metric { border: 1px solid #ddd; padding: 15px; text-align: center; flex: 1; margin: 0 5px; }
        .metric-title { font-weight: bold; }
        .metric-value { font-size: 24px; margin: 10px 0; }
        .tabs { margin-top: 20px; }
        .tab { display: inline-block; padding: 10px 20px; cursor: pointer; background: #f2f2f2; border: 1px solid #ddd; }
        .tab.active { background: white; border-bottom: 1px solid white; }
        .tab-content { display: none; border: 1px solid #ddd; padding: 15px; margin-top: -1px; }
        .tab-content.active { display: block; }
    </style>
    <script>
        function openTab(evt, tabName) {
            var i, tabcontent, tablinks;
            tabcontent = document.getElementsByClassName("tab-content");
            for (i = 0; i < tabcontent.length; i++) {
                tabcontent[i].style.display = "none";
            }
            tablinks = document.getElementsByClassName("tab");
            for (i = 0; i < tablinks.length; i++) {
                tablinks[i].className = tablinks[i].className.replace(" active", "");
            }
            document.getElementById(tabName).style.display = "block";
            evt.currentTarget.className += " active";
        }
    </script>
</head>
<body>
    <h1>Deconvolution, Paired Plot, and Stats Error Summary</h1>
    <p>Generated on $(date) for dataset: ${DATA_NAME}</p>
    
    <div class="tabs">
        <button class="tab active" onclick="openTab(event, 'DeconvolutionTab')">Deconvolution Errors</button>
        <button class="tab" onclick="openTab(event, 'PairedPlotTab')">Paired Plot Errors</button>
        <button class="tab" onclick="openTab(event, 'StatsTab')">Stats Errors</button>
        <button class="tab" onclick="openTab(event, 'SummaryTab')">Error Summary</button>
    </div>
    
    <div id="DeconvolutionTab" class="tab-content active">
        <h2>Deconvolution Run Errors</h2>
EOF

if [ ! -f "$JOB_MAPPING_FILE" ]; then
    cat >> "$OUTPUT_FILE" << EOF
        <p class="warning">Warning: Job mapping file not found at: ${JOB_MAPPING_FILE}</p>
EOF
else
    cat >> "$OUTPUT_FILE" << EOF
        <table>
            <tr>
                <th>Method</th>
                <th>Job ID</th>
                <th>Status</th>
                <th>Error Type</th>
                <th>Error Details</th>
            </tr>
EOF
fi

# Initialize counters for dashboard
TOTAL_METHODS=0
SUCCESS_METHODS=0
FAILED_METHODS=0
WARNING_METHODS=0

# Process each method and its job ID for deconvolution errors
if [ -f "$JOB_MAPPING_FILE" ]; then
    while IFS=: read -r METHOD JOB_ID; do
        # Skip comment lines
        [[ "$METHOD" =~ ^#.*$ ]] && continue
        
        TOTAL_METHODS=$((TOTAL_METHODS + 1))
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
            WARNING_METHODS=$((WARNING_METHODS + 1))
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
            SUCCESS_METHODS=$((SUCCESS_METHODS + 1))
            continue
        fi

        # Check for success markers in output without significant errors
        if grep -i "runtime:\|completed\|saved to\|check: deconvolution completed" "$ERROR_FILE" > /dev/null && ! grep -i "error\|failed\|fatal\|panic\|exception" "$ERROR_FILE" > /dev/null; then
            # File contains normal output without errors - consider as success
            cat >> "$OUTPUT_FILE" << EOF
            <tr>
                <td>${METHOD}</td>
                <td>${JOB_ID}</td>
                <td class="success">Success</td>
                <td>None</td>
                <td>Standard execution log (no errors)</td>
            </tr>
EOF
            SUCCESS_METHODS=$((SUCCESS_METHODS + 1))
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
        STATUS_CLASS="warning"
        STATUS_TEXT="Completed with warnings"
        
        # Check if there's a success message
        if grep -i "success\|completed successfully" "$ERROR_FILE" > /dev/null; then
            STATUS_CLASS="success"
            STATUS_TEXT="Success with warnings"
            SUCCESS_METHODS=$((SUCCESS_METHODS + 1))
        elif grep -i "error\|failed\|fatal\|panic" "$ERROR_FILE" > /dev/null; then
            STATUS_CLASS="error"
            STATUS_TEXT="Failed"
            FAILED_METHODS=$((FAILED_METHODS + 1))
        else
            WARNING_METHODS=$((WARNING_METHODS + 1))
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
fi

# Close deconvolution section and start paired plot section
cat >> "$OUTPUT_FILE" << EOF
        </table>
    </div>
    
    <div id="PairedPlotTab" class="tab-content">
        <h2>Paired Plot Errors</h2>
EOF

if [ ! -f "$PAIRED_MAPPING_FILE" ]; then
    cat >> "$OUTPUT_FILE" << EOF
        <p class="warning">Warning: Paired plot job mapping file not found at: ${PAIRED_MAPPING_FILE}</p>
EOF
else
    cat >> "$OUTPUT_FILE" << EOF
        <table>
            <tr>
                <th>Method</th>
                <th>Plot Job ID</th>
                <th>Status</th>
                <th>Error Type</th>
                <th>Error Details</th>
            </tr>
EOF
fi

# Initialize counters for paired plot jobs
PLOT_TOTAL=0
PLOT_SUCCESS=0
PLOT_FAILED=0
PLOT_WARNING=0

# Add code to track paired plot job errors
if [ -f "$PAIRED_MAPPING_FILE" ]; then
    while IFS=: read -r METHOD PLOT_JOB_ID; do
        # Skip comment lines
        [[ "$METHOD" =~ ^#.*$ ]] && continue
        
        PLOT_TOTAL=$((PLOT_TOTAL + 1))
        
        # Check plot job error file
        PLOT_ERROR_FILE="$LOG_DIR/${PLOT_JOB_ID}.e"
        if [ ! -f "$PLOT_ERROR_FILE" ]; then
            cat >> "$OUTPUT_FILE" << EOF
        <tr>
            <td>${METHOD}</td>
            <td>${PLOT_JOB_ID}</td>
            <td class="warning">Unknown</td>
            <td>N/A</td>
            <td>Error log file not found</td>
        </tr>
EOF
            PLOT_WARNING=$((PLOT_WARNING + 1))
            continue
        fi
        
    # Check for successful completion
    if [ ! -s "$PLOT_ERROR_FILE" ]; then
        # Empty error file means success with no errors
        cat >> "$OUTPUT_FILE" << EOF
    <tr>
        <td>${METHOD}</td>
        <td>${PLOT_JOB_ID}</td>
        <td class="success">Success</td>
        <td>None</td>
        <td>No errors detected</td>
    </tr>
EOF
        PLOT_SUCCESS=$((PLOT_SUCCESS + 1))
        continue
    fi

    # Check for success markers in output
    if grep -i "saved to\|completed\|runtime:\|plot saved" "$PLOT_ERROR_FILE" > /dev/null && ! grep -i "error\|failed\|fatal\|panic\|exception" "$PLOT_ERROR_FILE" > /dev/null; then
        # File contains normal output without errors - consider as success
        cat >> "$OUTPUT_FILE" << EOF
    <tr>
        <td>${METHOD}</td>
        <td>${PLOT_JOB_ID}</td>
        <td class="success">Success</td>
        <td>None</td>
        <td>Standard execution log (no errors)</td>
    </tr>
EOF
        PLOT_SUCCESS=$((PLOT_SUCCESS + 1))
        continue
    fi
        
        # Search for error patterns
        ERROR_FOUND=false
        ERROR_TYPE=""
        ERROR_DETAILS=""
        
        while IFS= read -r PATTERN; do
            if grep -i "$PATTERN" "$PLOT_ERROR_FILE" > /dev/null; then
                ERROR_FOUND=true
                ERROR_TYPE="$PATTERN"
                ERROR_LINE=$(grep -i -m 1 "$PATTERN" "$PLOT_ERROR_FILE" | tr -d '\r' | sed 's/</\&lt;/g; s/>/\&gt;/g')
                LINE_NUM=$(grep -i -n -m 1 "$PATTERN" "$PLOT_ERROR_FILE" | cut -d: -f1)
                if [ ! -z "$LINE_NUM" ]; then
                    START_LINE=$((LINE_NUM > 1 ? LINE_NUM - 1 : 1))
                    END_LINE=$((LINE_NUM + 1))
                    ERROR_CONTEXT=$(sed -n "${START_LINE},${END_LINE}p" "$PLOT_ERROR_FILE" | tr -d '\r' | sed 's/</\&lt;/g; s/>/\&gt;/g')
                    ERROR_DETAILS=$(echo "$ERROR_CONTEXT" | sed 's/$/<br>/')
                else
                    ERROR_DETAILS="$ERROR_LINE"
                fi
                break
            fi
        done < "$ERROR_PATTERNS_FILE"
        
        # If no specific error found but file exists
        if [ "$ERROR_FOUND" = false ] && [ -s "$PLOT_ERROR_FILE" ]; then
            ERROR_TYPE="Unknown"
            ERROR_DETAILS=$(tail -5 "$PLOT_ERROR_FILE" | sed 's/</\&lt;/g; s/>/\&gt;/g; s/$/<br>/')
        fi
        
        # Determine status
        STATUS_CLASS="warning"
        STATUS_TEXT="Completed with warnings"
        
        if grep -i "success\|completed successfully" "$PLOT_ERROR_FILE" > /dev/null; then
            STATUS_CLASS="success"
            STATUS_TEXT="Success with warnings"
            PLOT_SUCCESS=$((PLOT_SUCCESS + 1))
        elif grep -i "error\|failed\|fatal\|panic" "$PLOT_ERROR_FILE" > /dev/null; then
            STATUS_CLASS="error"
            STATUS_TEXT="Failed"
            PLOT_FAILED=$((PLOT_FAILED + 1))
        else
            PLOT_WARNING=$((PLOT_WARNING + 1))
        fi
        
        cat >> "$OUTPUT_FILE" << EOF
        <tr>
            <td>${METHOD}</td>
            <td>${PLOT_JOB_ID}</td>
            <td class="${STATUS_CLASS}">${STATUS_TEXT}</td>
            <td>${ERROR_TYPE}</td>
            <td>${ERROR_DETAILS}</td>
        </tr>
EOF

    done < "$PAIRED_MAPPING_FILE"
else
    # Fallback message if mapping file doesn't exist
    cat >> "$OUTPUT_FILE" << EOF
        <p class="warning">No paired plot job mapping file found at ${PAIRED_MAPPING_FILE}</p>
EOF
fi

# Close paired plot section and start stats section
cat >> "$OUTPUT_FILE" << EOF
        </table>
    </div>
    
    <div id="StatsTab" class="tab-content">
        <h2>Statistics Calculation Errors</h2>
EOF

if [ ! -f "$STATS_MAPPING_FILE" ]; then
    cat >> "$OUTPUT_FILE" << EOF
        <p class="warning">Warning: Stats job mapping file not found at: ${STATS_MAPPING_FILE}</p>
EOF
else
    cat >> "$OUTPUT_FILE" << EOF
        <table>
            <tr>
                <th>Job Name</th>
                <th>Stats Job ID</th>
                <th>Status</th>
                <th>Error Type</th>
                <th>Error Details</th>
            </tr>
EOF
fi

# Initialize counters for stats jobs
STATS_TOTAL=0
STATS_SUCCESS=0
STATS_FAILED=0
STATS_WARNING=0

# Add code to track stats job errors
if [ -f "$STATS_MAPPING_FILE" ]; then
    while IFS=: read -r JOB_NAME STATS_JOB_ID; do
        # Skip comment lines
        [[ "$JOB_NAME" =~ ^#.*$ ]] && continue
        
        STATS_TOTAL=$((STATS_TOTAL + 1))
        
        # Check stats job error file
        STATS_ERROR_FILE="$LOG_DIR/${STATS_JOB_ID}.e"
        if [ ! -f "$STATS_ERROR_FILE" ]; then
            cat >> "$OUTPUT_FILE" << EOF
        <tr>
            <td>${JOB_NAME}</td>
            <td>${STATS_JOB_ID}</td>
            <td class="warning">Unknown</td>
            <td>N/A</td>
            <td>Error log file not found</td>
        </tr>
EOF
            STATS_WARNING=$((STATS_WARNING + 1))
            continue
        fi
        
        # Check for successful completion
        if [ ! -s "$STATS_ERROR_FILE" ]; then
            # Empty error file means success with no errors
            cat >> "$OUTPUT_FILE" << EOF
        <tr>
            <td>${JOB_NAME}</td>
            <td>${STATS_JOB_ID}</td>
            <td class="success">Success</td>
            <td>None</td>
            <td>No errors detected</td>
        </tr>
EOF
            STATS_SUCCESS=$((STATS_SUCCESS + 1))
            continue
        fi

        # Check for success markers in output
        if grep -i "saved to\|completed\|runtime:\|benchmark" "$STATS_ERROR_FILE" > /dev/null && ! grep -i "error\|failed\|fatal\|panic\|exception" "$STATS_ERROR_FILE" > /dev/null; then
            # File contains normal output without errors - consider as success
            cat >> "$OUTPUT_FILE" << EOF
        <tr>
            <td>${JOB_NAME}</td>
            <td>${STATS_JOB_ID}</td>
            <td class="success">Success</td>
            <td>None</td>
            <td>Standard execution log (no errors)</td>
        </tr>
EOF
            STATS_SUCCESS=$((STATS_SUCCESS + 1))
            continue
        fi
            
        # Search for error patterns
        ERROR_FOUND=false
        ERROR_TYPE=""
        ERROR_DETAILS=""
        
        while IFS= read -r PATTERN; do
            if grep -i "$PATTERN" "$STATS_ERROR_FILE" > /dev/null; then
                ERROR_FOUND=true
                ERROR_TYPE="$PATTERN"
                ERROR_LINE=$(grep -i -m 1 "$PATTERN" "$STATS_ERROR_FILE" | tr -d '\r' | sed 's/</\&lt;/g; s/>/\&gt;/g')
                LINE_NUM=$(grep -i -n -m 1 "$PATTERN" "$STATS_ERROR_FILE" | cut -d: -f1)
                if [ ! -z "$LINE_NUM" ]; then
                    START_LINE=$((LINE_NUM > 1 ? LINE_NUM - 1 : 1))
                    END_LINE=$((LINE_NUM + 1))
                    ERROR_CONTEXT=$(sed -n "${START_LINE},${END_LINE}p" "$STATS_ERROR_FILE" | tr -d '\r' | sed 's/</\&lt;/g; s/>/\&gt;/g')
                    ERROR_DETAILS=$(echo "$ERROR_CONTEXT" | sed 's/$/<br>/')
                else
                    ERROR_DETAILS="$ERROR_LINE"
                fi
                break
            fi
        done < "$ERROR_PATTERNS_FILE"
        
        # If no specific error found but file exists
        if [ "$ERROR_FOUND" = false ] && [ -s "$STATS_ERROR_FILE" ]; then
            ERROR_TYPE="Unknown"
            ERROR_DETAILS=$(tail -5 "$STATS_ERROR_FILE" | sed 's/</\&lt;/g; s/>/\&gt;/g; s/$/<br>/')
        fi
        
        # Determine status
        STATUS_CLASS="warning"
        STATUS_TEXT="Completed with warnings"
        
        if grep -i "success\|completed successfully" "$STATS_ERROR_FILE" > /dev/null; then
            STATUS_CLASS="success"
            STATUS_TEXT="Success with warnings"
            STATS_SUCCESS=$((STATS_SUCCESS + 1))
        elif grep -i "error\|failed\|fatal\|panic" "$STATS_ERROR_FILE" > /dev/null; then
            STATUS_CLASS="error"
            STATUS_TEXT="Failed"
            STATS_FAILED=$((STATS_FAILED + 1))
        else
            STATS_WARNING=$((STATS_WARNING + 1))
        fi
        
        cat >> "$OUTPUT_FILE" << EOF
        <tr>
            <td>${JOB_NAME}</td>
            <td>${STATS_JOB_ID}</td>
            <td class="${STATUS_CLASS}">${STATUS_TEXT}</td>
            <td>${ERROR_TYPE}</td>
            <td>${ERROR_DETAILS}</td>
        </tr>
EOF

    done < "$STATS_MAPPING_FILE"
else
    # Fallback message if mapping file doesn't exist
    cat >> "$OUTPUT_FILE" << EOF
        <p class="warning">No stats job mapping file found at ${STATS_MAPPING_FILE}</p>
EOF
fi

# Close stats section and start summary section
cat >> "$OUTPUT_FILE" << EOF
        </table>
    </div>
    
    <div id="SummaryTab" class="tab-content">
        <h2>Error Pattern Summary</h2>
        
        <div class="dashboard">
            <div class="metric">
                <div class="metric-title">Deconv Methods</div>
                <div class="metric-value">${TOTAL_METHODS}</div>
            </div>
            <div class="metric">
                <div class="metric-title">Success</div>
                <div class="metric-value" style="color: #2ecc71">${SUCCESS_METHODS}</div>
                <div>$(printf "%.1f%%" $((SUCCESS_METHODS * 100 / (TOTAL_METHODS > 0 ? TOTAL_METHODS : 1))))</div>
            </div>
            <div class="metric">
                <div class="metric-title">Warning</div>
                <div class="metric-value" style="color: #f39c12">${WARNING_METHODS}</div>
                <div>$(printf "%.1f%%" $((WARNING_METHODS * 100 / (TOTAL_METHODS > 0 ? TOTAL_METHODS : 1))))</div>
            </div>
            <div class="metric">
                <div class="metric-title">Failed</div>
                <div class="metric-value" style="color: #e74c3c">${FAILED_METHODS}</div>
                <div>$(printf "%.1f%%" $((FAILED_METHODS * 100 / (TOTAL_METHODS > 0 ? TOTAL_METHODS : 1))))</div>
            </div>
        </div>
        
        <div class="dashboard">
            <div class="metric">
                <div class="metric-title">Plots Total</div>
                <div class="metric-value">${PLOT_TOTAL}</div>
            </div>
            <div class="metric">
                <div class="metric-title">Plots Success</div>
                <div class="metric-value" style="color: #2ecc71">${PLOT_SUCCESS}</div>
                <div>$(printf "%.1f%%" $((PLOT_SUCCESS * 100 / (PLOT_TOTAL > 0 ? PLOT_TOTAL : 1))))</div>
            </div>
            <div class="metric">
                <div class="metric-title">Plots Warning</div>
                <div class="metric-value" style="color: #f39c12">${PLOT_WARNING}</div>
                <div>$(printf "%.1f%%" $((PLOT_WARNING * 100 / (PLOT_TOTAL > 0 ? PLOT_TOTAL : 1))))</div>
            </div>
            <div class="metric">
                <div class="metric-title">Plots Failed</div>
                <div class="metric-value" style="color: #e74c3c">${PLOT_FAILED}</div>
                <div>$(printf "%.1f%%" $((PLOT_FAILED * 100 / (PLOT_TOTAL > 0 ? PLOT_TOTAL : 1))))</div>
            </div>
        </div>
        
        <div class="dashboard">
            <div class="metric">
                <div class="metric-title">Stats Total</div>
                <div class="metric-value">${STATS_TOTAL}</div>
            </div>
            <div class="metric">
                <div class="metric-title">Stats Success</div>
                <div class="metric-value" style="color: #2ecc71">${STATS_SUCCESS}</div>
                <div>$(printf "%.1f%%" $((STATS_SUCCESS * 100 / (STATS_TOTAL > 0 ? STATS_TOTAL : 1))))</div>
            </div>
            <div class="metric">
                <div class="metric-title">Stats Warning</div>
                <div class="metric-value" style="color: #f39c12">${STATS_WARNING}</div>
                <div>$(printf "%.1f%%" $((STATS_WARNING * 100 / (STATS_TOTAL > 0 ? STATS_TOTAL : 1))))</div>
            </div>
            <div class="metric">
                <div class="metric-title">Stats Failed</div>
                <div class="metric-value" style="color: #e74c3c">${STATS_FAILED}</div>
                <div>$(printf "%.1f%%" $((STATS_FAILED * 100 / (STATS_TOTAL > 0 ? STATS_TOTAL : 1))))</div>
            </div>
        </div>
        
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
    
    # Count occurrences of this pattern across all logs
    METHODS_WITH_ERROR=""
    ERROR_COUNT=0
    
    # Check deconvolution jobs
    if [ -f "$JOB_MAPPING_FILE" ]; then
        while IFS=: read -r METHOD JOB_ID; do
            # Skip comment lines
            [[ "$METHOD" =~ ^#.*$ ]] && continue
            
            ERROR_FILE="$LOG_DIR/${JOB_ID}.e"
            if [ -f "$ERROR_FILE" ] && grep -i "$PATTERN" "$ERROR_FILE" > /dev/null; then
                METHODS_WITH_ERROR="$METHODS_WITH_ERROR$METHOD (deconv), "
                ERROR_COUNT=$((ERROR_COUNT + 1))
            fi
        done < "$JOB_MAPPING_FILE"
    fi
        
    # Check paired plot jobs
    if [ -f "$PAIRED_MAPPING_FILE" ]; then
        while IFS=: read -r METHOD PLOT_JOB_ID; do
            # Skip comment lines
            [[ "$METHOD" =~ ^#.*$ ]] && continue
            
            if [ ! -z "$PLOT_JOB_ID" ]; then
                PLOT_ERROR_FILE="$LOG_DIR/${PLOT_JOB_ID}.e"
                if [ -f "$PLOT_ERROR_FILE" ] && grep -i "$PATTERN" "$PLOT_ERROR_FILE" > /dev/null; then
                    METHODS_WITH_ERROR="$METHODS_WITH_ERROR$METHOD (plot), "
                    ERROR_COUNT=$((ERROR_COUNT + 1))
                fi
            fi
        done < "$PAIRED_MAPPING_FILE"
    fi
    
    # Check stats jobs
    if [ -f "$STATS_MAPPING_FILE" ]; then
        while IFS=: read -r JOB_NAME STATS_JOB_ID; do
            # Skip comment lines
            [[ "$JOB_NAME" =~ ^#.*$ ]] && continue
            
            if [ ! -z "$STATS_JOB_ID" ]; then
                STATS_ERROR_FILE="$LOG_DIR/${STATS_JOB_ID}.e"
                if [ -f "$STATS_ERROR_FILE" ] && grep -i "$PATTERN" "$STATS_ERROR_FILE" > /dev/null; then
                    METHODS_WITH_ERROR="$METHODS_WITH_ERROR$JOB_NAME (stats), "
                    ERROR_COUNT=$((ERROR_COUNT + 1))
                fi
            fi
        done < "$STATS_MAPPING_FILE"
    fi
    
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
    </div>
    
    <script>
        // Open first tab by default
        document.getElementsByClassName("tab")[0].click();
    </script>
</body>
</html>
EOF

echo "Error summary generated at: $OUTPUT_FILE"

# Clean up
rm -rf "$TEMP_DIR"