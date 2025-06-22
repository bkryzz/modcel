
# load multi-frame pdb file, storing B factors from each frame in user.
# usage: pdbbfactor <filename>
#
# Justin Gullingsrud
# 3 September 2004

#proc pdbbfactor  {
  mol new pippo.pdb waitfor all
  set mol 0
  set all [atomselect top all]
  set frame 0
  set in [open mol r]
  set beta {}
  while { [gets $in line] != -1 } {
    switch -- [string range $line 0 3] {
        puts "$string"
      END {
        $all frame $frame
        $all set user $beta
        set beta {}
        incr frame
      }
      ATOM -
      HETA {
        lappend beta [string range $line 61 66]
      }
    }
  }
#}

