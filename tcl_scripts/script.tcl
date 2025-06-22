  set filelist [glob maps*.pdb] 
  mol default color beta
  mol default representation vdw
  foreach file $filelist { 
    mol new $file waitfor all 
  }
  for {set n 0} {$n < 40} {incr n} { 
    mol scaleminmax $n 0 0.0 4.0
  }
