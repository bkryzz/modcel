# load multi-frame pdb file, storing occupancy factors 'my scripts store centro-symmetry parameter in occupancy slot' from each frame in user.
# usage: pdbbfactor <filename>
# Justin Gullingsrud & Modified by cjo to take input from occupancy field
# 3 September 2004 & 1/2007
proc factor { fname } {
   mol new $fname waitfor all
   set all [atomselect top all]
   set frame 0
   set in [open $fname r]
   set occupancy {}
   while { [gets $in line] != -1 } {
     switch -- [string range $line 0 3] {
       END {
         $all frame $frame
         $all set user $beta
         set beta {}
         incr frame
       }
       ATOM {
         lappend beta [string range $line 61 66]
       }
         puts "$user"
     }
   }
 }
