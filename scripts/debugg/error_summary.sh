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
OUTPUT_FILE="/work/gr-fe/lorthiois/DeconBenchmark/error_summary.html"
JOB_MAPPING_FILE="/work/gr-fe/lorthiois/DeconBenchmark/logs/model_job_mapping.txt"
PAIRED_MAPPING_FILE="/work/gr-fe/lorthiois/DeconBenchmark/logs/paired_plot_job_mapping.txt"
DATA_NAME=$(basename $(ls /work/gr-fe/lorthiois/DeconBenchmark/deconv_results/results_*_*.rda 2>/dev/null | head -1) | sed -E 's/results_[^_]+_(.+)\.rda/\1/' || echo "unknown")

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
    exit 1
fi

# Generate HTML header
cat > "$OUTPUT_FILE" << EOF
<!DOCTYPE html>
<html>
<head>
    <title>Deconvolution and Paired Plot Error Summary</title>
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
    <h1>Deconvolution and Paired Plot Error Summary</h1>
    <p>Generated on $(date) for dataset: ${DATA_NAME}</p>
    
    <div class="tabs">
        <button class="tab active" onclick="openTab(event, 'DeconvolutionTab')">Deconvolution Errors</button>
        <button class="tab" onclick="openTab(event, 'PairedPlotTab')">Paired Plot Errors</button>
        <button class="tab" onclick="openTab(event, 'SummaryTab')">Error Summary</button>
    </div>
    
    <div id="DeconvolutionTab" class="tab-content active">
        <h2>Deconvolution Run Errors</h2>
        <table>
            <tr>
                <th>Method</th>
                <th>Job ID</th>
                <th>Status</th>
                <th>Error Type</th>
                <th>Error Details</th>
            </tr>
EOF

# Initialize counters for dashboard
TOTAL_METHODS=0
SUCCESS_METHODS=0
FAILED_METHODS=0
WARNING_METHODS=0

# Process each method and its job ID for deconvolution errors
while IFS=: read -r METHOD JOB_ID; do
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

# Close deconvolution section
cat >> "$OUTPUT_FILE" << EOF
        </table>
    </div>
    
    <div id="PairedPlotTab" class="tab-content">
        <h2>Paired Plot Errors</h2>
        <table>
            <tr>
                <th>Method</th>
                <th>Plot Job ID</th>
                <th>Status</th>
                <th>Error Type</th>
                <th>Error Details</th>
            </tr>
EOF

# Initialize counters for paired plot jobs
PLOT_TOTAL=0
PLOT_SUCCESS=0
PLOT_FAILED=0
PLOT_WARNING=0

# Now add code to track paired plot job errors
if [ -f "$PAIRED_MAPPING_FILE" ]; then
    while IFS=: read -r METHOD PLOT_JOB_ID; do
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
        <tr>
            <td colspan="5" class="warning">No paired plot job mapping file found at ${PAIRED_MAPPING_FILE}</td>
        </tr>
EOF
fi

# Close paired plot section
cat >> "$OUTPUT_FILE" << EOF
        </table>
    </div>
    
    <div id="SummaryTab" class="tab-content">
        <h2>Error Pattern Summary</h2>
        
        <div class="dashboard">
            <div class="metric">
                <div class="metric-title">Total Methods</div>
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
    while IFS=: read -r METHOD JOB_ID; do
        ERROR_FILE="$LOG_DIR/${JOB_ID}.e"
        if [ -f "$ERROR_FILE" ] && grep -i "$PATTERN" "$ERROR_FILE" > /dev/null; then
            METHODS_WITH_ERROR="$METHODS_WITH_ERROR$METHOD (deconv), "
            ERROR_COUNT=$((ERROR_COUNT + 1))
        fi
    done < "$JOB_MAPPING_FILE"
        
    # Check paired plot jobs
    if [ -f "$PAIRED_MAPPING_FILE" ]; then
        while IFS=: read -r METHOD PLOT_JOB_ID; do
            if [ ! -z "$PLOT_JOB_ID" ]; then
                PLOT_ERROR_FILE="$LOG_DIR/${PLOT_JOB_ID}.e"
                if [ -f "$PLOT_ERROR_FILE" ] && grep -i "$PATTERN" "$PLOT_ERROR_FILE" > /dev/null; then
                    METHODS_WITH_ERROR="$METHODS_WITH_ERROR$METHOD (plot), "
                    ERROR_COUNT=$((ERROR_COUNT + 1))
                fi
            fi
        done < "$PAIRED_MAPPING_FILE"
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