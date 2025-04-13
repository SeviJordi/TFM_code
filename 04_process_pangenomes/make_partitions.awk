#!/usr/bin/awk
 
/^FT   feature/ { 
        start = $2; 
        end = $3; 
        gsub(/\.\./, "-", end)  # Convert ".." to "-"
    }
/^FT   misc_feature/ { 
        start = $2; 
        end = $3; 
        gsub(/\.\./, "-", end)  # Convert ".." to "-"
    }   
/\/label=/ { 
        gsub("/label=", "", $0);  
        gsub(/\.aln/, "", $0);  
        print $2 " = " end
    } 