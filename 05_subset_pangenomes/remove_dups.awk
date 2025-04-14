#!/usr/bin/awk

BEGIN {
    FS = "\t";    # Set input field separator to tab
    OFS = "\t";   # Set output field separator to tab
}

{
    counts[$column]++;               # Count occurrences of each value in the specified column
    rows[NR] = $0;                  # Store the row for potential output
    column_values[NR] = $column;    # Store the column value for each row
}

END {
    for (i = 1; i <= NR; i++) {
        if (counts[column_values[i]] == 1) {
            print rows[i];          # Print rows where the value in the column is unique
        }
    }
}