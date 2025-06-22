set mol 1
set all [atomselect top all]

set nframes [molinfo $mol get numframes]

for {set i 0} {$i < $nframes} {incr i} {
	$all frame $i
	set dlist ""
        foreach ATOM [$all string range 61 66] {
		lappend dlist [$x $y $z]
		puts "$dlist"
	}
	$all set user $dlist
}
