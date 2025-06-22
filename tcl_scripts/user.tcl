proc factor { fname } {
  set all [atomselect top all]
  set frame 0
  set infile [open $fname r] 
  while {[gets $infile line] >=0} {
    switch -- [string range $line 0 3] {
    END {
     $all frame $frame
     # reads 1st column 
     set value [lindex $line 0] 
     # put $value in user field 
     $all set user $value
     incr frame
        }
     }
  }
}
