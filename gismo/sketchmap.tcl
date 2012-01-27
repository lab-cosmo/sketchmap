#
# An interface to michelino's sketch map routines
#

package provide sketchmap 1.0

namespace eval ::sketchmap {
   namespace export retrieveParams       ; # Return all the parameters required for out of sample embedding on a list 
   namespace export setParameters        ; # Set parameters for out of sample embedding
   namespace export setPFileParameters   ; # Set parameters for creating projection files
   namespace export MDSWindow            ; # Create the MDS window
   namespace export sketchMapWindow      ; # Create the sketch map window
   namespace export outOfSampleWindow    ; # Create the outOfSample window
   namespace export projectionFileWindow ; # Create the window for producing projection files
   namespace export printProjectionData  ; # Print data to projection file (bespoke.defs)
   namespace export runSketchMap         ; # Actually runs sketch map
   namespace export runOutOfSample       ; # Actually run out of sample embedding
   namespace export setProgressBar       ; # Stores the location of the progress bar 
# Must export these variables so they can be used in windows
   namespace export dimredConfig
   namespace export osampleConfig
   namespace export pfileConfig
# The actual variables we use to make things go
   variable  dimredConfig                ; # Configuration for dimensionality reduction 
   variable  osampleConfig               ; # Configuration for out of sample embedding
   variable  pfileConfig                 ; # Configuration for output of projection file
   variable  progressBar                 ; # The progress bar for the currently active process   
   variable  tmpdir                      ; # The tempory directory for out of sample calculations
}

proc sketchmap::retrieveParams {args} {
   variable dimredConfig

   set type [getargs $args "-type" "No -type flag in call to createDimRedGoButtons"]

   set paramlist {}
   lappend paramlist "TYPE" 
   lappend paramlist $type
   if { $type=="smap" } {
      lappend paramlist "SIGMA" ; lappend paramlist $dimredConfig(sigma)
      lappend paramlist "aD"    ; lappend paramlist $dimredConfig(aD)
      lappend paramlist "bD"    ; lappend paramlist $dimredConfig(bD)
      lappend paramlist "ad"    ; lappend paramlist $dimredConfig(ad)
      lappend paramlist "bd"    ; lappend paramlist $dimredConfig(bd)
   }
   return $paramlist
}

proc sketchmap::setParameters {args} {
   variable osampleConfig

   set parameters [getargs $args "-parameters" "No -parameters flag in call to setParameters"]

   set osampleConfig(rtype) [retrieveParam $parameters -name "TYPE"]
   if { $osampleConfig(rtype)=="smap" } {
      set osampleConfig(sigma) [retrieveParam $parameters -name "SIGMA"]
      set osampleConfig(aD) [retrieveParam $parameters -name "aD"]
      set osampleConfig(bD) [retrieveParam $parameters -name "bD"]
      set osampleConfig(ad) [retrieveParam $parameters -name "ad"]
      set osampleConfig(bd) [retrieveParam $parameters -name "bd"] 
   }
}

proc sketchmap::setPfileParameters {args} {
   variable pfileConfig

   set parameters [getargs $args "-parameters" "No -parameters flag in call to setParameters"]
   
   set pfileConfig(rtype) [retrieveParam $parameters -name "TYPE"]
   if { $pfileConfig(rtype)=="smap" } {
      set pfileConfig(sigma) [retrieveParam $parameters -name "SIGMA"]
      set pfileConfig(aD) [retrieveParam $parameters -name "aD"]
      set pfileConfig(bD) [retrieveParam $parameters -name "bD"]
      set pfileConfig(ad) [retrieveParam $parameters -name "ad"]
      set pfileConfig(bd) [retrieveParam $parameters -name "bd"] 
   }
}

proc sketchmap::retrieveParam {params args} {
   set name [getargs $args "-name" "No -name flag in call to retrieveParam"]
   set pos [ lsearch $params $name ]
   if { $pos<0 } { puts "Something terrible has gone wrong - can't find $name parameter" }
   return [ lindex $params [expr $pos+1] ]
}

proc sketchmap::runOutOfSample {args} {
   variable osampleConfig
   global env

   set ncv [getargs $args "-ncv" "No -ncv flag in call to runOutOfSample"]
   set nlandmarks [getargs $args "-nlandmarks" "No -nlandmarks flag in all to runOutOfSample"]
   set npoints [getargs $args "-npoints" "No -npoints flag in call to runOutOfSample"]
   set ifile [getargs $args "-ifile" "No -ifile flag in call to runOutOfSample"]
   set lhfile [getargs $args "-lhfile" "No -lhfile flag in call to runOutOfSample"]
   set llfile [getargs $args "-llfile" "No -llfile flag in call to runOutOfSample"]

   if { $osampleConfig(useweights)==1 } {
     set wflag "-w"
   } else {
     set wflag ""
   }

   if { $osampleConfig(isperiodic)==1 } {
      set pflag "-pi $osampleConfig(period)"
   } else {
      set pflag ""
   }

   # Make sure we have valid input for all required stuff
   if { ![ string is integer -strict $osampleConfig(nconjgrad) ] } {
      tk_messageBox -icon error -type ok -title Message -message "Invalid value for number of steps of conjugate gradient"
      return 0
   }
   if { ![ string is integer -strict $osampleConfig(ngridcoarse) ] } {
      tk_messageBox -icon error -type ok -title Message -message "Invalid value for number of points in coarse grained grid"
      return 0
   }
   if { ![ string is integer -strict $osampleConfig(ngridfine) ] } {
      tk_messageBox -icon error -type ok -title Message -message "Invalid value for number of points in fine grained grid"
      return 0
   }
  
   # Check we have all skech-map parameters set
   if { $osampleConfig(rtype)=="smap" } {
      if { ![string is double -strict $osampleConfig(sigma) ] } {
         tk_messageBox -icon error -type ok -title Message -message "Invalid value for sigma"
         return 0
      }
      if { ![ string is integer -strict $osampleConfig(aD) ] } {
         tk_messageBox -icon error -type ok -title Message -message "Invalid value for aD"
         return 0
      }  
      if { ![ string is integer -strict $osampleConfig(ad) ] } {
         tk_messageBox -icon error -type ok -title Message -message "Invalid value for ad"
         return 0
      }  
      if { ![ string is integer -strict $osampleConfig(bD) ] } {
         tk_messageBox -icon error -type ok -title Message -message "Invalid value for bD"
         return 0
      }  
      if { ![ string is integer -strict $osampleConfig(bd) ] } {
         tk_messageBox -icon error -type ok -title Message -message "Invalid value for bd"
         return 0
      }  
   }

   # Retrive the width of the grid based on the spread of landmark points
   if { [readMikFile -ifile $llfile -errors "Something has gone badly wrong for out of sample" -npoints $nlandmarks -projname TMP]==0 } { return 0 }
   set gwidth [ expr  1.1 * [ retrieveDataWidth ] ]

   # Run out of sample projection algorithm 
   if { $osampleConfig(rtype)=="mds" } { 
      puts "Executing: dimproj -D $ncv -d 2 $pflag $wflag -P $lhfile -p $llfile -g1 $osampleConfig(ngridcoarse) -g2 $osampleConfig(ngridfine) -gw $gwidth -fun identity -cgmin $osampleConfig(nconjgrad) < $ifile > OSAMPLES"
      catch { exec $env(smapdir)/bin/dimproj -D $ncv -d 2 $pflag $wflag -P $lhfile -p $llfile -g1 $osampleConfig(ngridcoarse) -g2 $osampleConfig(ngridfine) -gw $gwidth -fun identity -cgmin $osampleConfig(nconjgrad) < $ifile > OSAMPLES } osample_stdout
      if { [readMikFile -ifile OSAMPLES -errors $osample_stdout -npoints $npoints -projname MDS]==0 } { return 0 }
   } elseif { $osampleConfig(rtype)=="smap" } { 
      puts "Executing: dimproj -D $ncv -d 2 $pflag $wflag -P $lhfile -p $llfile -g1 $osampleConfig(ngridcoarse) -g2 $osampleConfig(ngridfine) -gw $gwidth -fun xsigmoid -sigma $osampleConfig(sigma) -expa $osampleConfig(aD) -lexpa $osampleConfig(ad) -expb $osampleConfig(bD) -lexpb $osampleConfig(bd) -cgmin $osampleConfig(nconjgrad) < $ifile > OSAMPLES"
      catch { exec $env(smapdir)/bin/dimproj -D $ncv -d 2 $pflag $wflag -P $lhfile -p $llfile -g1 $osampleConfig(ngridcoarse) -g2 $osampleConfig(ngridfine) -gw $gwidth -fun xsigmoid -sigma $osampleConfig(sigma) -expa $osampleConfig(aD) -lexpa $osampleConfig(ad) -expb $osampleConfig(bD) -lexpb $osampleConfig(bd) -cgmin $osampleConfig(nconjgrad) < $ifile > OSAMPLES } osample_stdout
      if { [readMikFile -ifile OSAMPLES -errors $osample_stdout -npoints $npoints -projname SMAP]==0 } { return 0 }
   } else {
      puts "Invalid osampleConfig(rtype) ($osampleConfig(rtype)) in call to runOutOfSample"
      return 0
   }
   return 1 
}

proc sketchmap::runSketchMap {args} {
   variable dimredConfig
   variable progressBar
   global env
 
   set type [getargs $args "-type" "No -type flag in call to createDimRedGoButtons"]
   set ncv [getargs $args "-ncv" "No -ncv flag in call to runSketchMap"]
   set ifile [getargs $args "-ifile" "No -ifile flag in call to runSketchMap"]

  if { $type=="mds" && $dimredConfig(noiter)==1 } {
     ::progressbar::config $progressBar -ntasks 1
  } elseif { $type=="mds" } {
     ::progressbar::config $progressBar -ntasks 2
  } else {
     ::progressbar::config $progressBar -ntasks 7
  }

  if { $dimredConfig(useweights)==1 } {
     set wflag "-w"
  } else {
     set wflag ""
  }

  if { $dimredConfig(isperiodic)==1 } {
     set pflag "-pi $dimredConfig(period)"
  } else {
     set pflag ""
  }

  # Run MDS algorithm with eigenvector/eigenvalue solution
  puts "Executing: dimred -D $ncv -d 2 $pflag -mode MDS < $ifile > MDS_DATA"
  catch { exec $env(smapdir)/bin/dimred -D $ncv -d 2 $pflag $wflag -mode MDS < $ifile > MDS_DATA } dimr_stdout
  if { [readMikFile -ifile MDS_DATA -errors $dimr_stdout -npoints [molinfo top get numframes] -projname MDS]==0 } { return 0 }

  # Increment the progress bar
  ::progressbar::increment $progressBar

  if { $type=="mds" && $dimredConfig(noiter)==1 } { return 1 }

  # Make sure we have valid input for all required stuff
  if { ![ string is integer -strict $dimredConfig(nconjgrad) ] } {
     tk_messageBox -icon error -type ok -title Message -message "Invalid value for number of steps of conjugate gradient"
     return 0
  }
  if { ![ string is integer -strict $dimredConfig(ngridcoarse) ] } {
     tk_messageBox -icon error -type ok -title Message -message "Invalid value for number of points in coarse grained grid"
     return 0
  }
  if { ![ string is integer -strict $dimredConfig(ngridfine) ] } {
     tk_messageBox -icon error -type ok -title Message -message "Invalid value for number of points in fine grained grid"
     return 0
  }

  # Retrive the width of the grid
  set gwidth [ expr  1.1 * [ retrieveDataWidth ] ]

  # Run iterative MDS algorithm
  puts "Executing: dimred -D $ncv -d 2 $pflag -mode GMDS -fun identity -steps $dimredConfig(nconjgrad) -imode conjgrad $wflag -g1 $dimredConfig(ngridcoarse) -g2 $dimredConfig(ngridfine) -gw $gwidth < $ifile > IMDS_DATA"
  catch { exec $env(smapdir)/bin/dimred -D $ncv -d 2 $pflag -mode GMDS -fun identity -steps $dimredConfig(nconjgrad) -imode conjgrad $wflag -g1 $dimredConfig(ngridcoarse) -g2 $dimredConfig(ngridfine) -gw $gwidth < $ifile > IMDS_DATA } dimr_stdout 
  
  # Read the file containing the output data
  if { [readMikFile -ifile IMDS_DATA -errors $dimr_stdout -npoints [molinfo top get numframes] -projname IMDS ]==0 } { return 0 }

  # Increment the progress bar
  ::progressbar::increment $progressBar

  if { $type=="mds" } { return 1 }

  # Check we have all skech-map parameters set
  if { ![string is double -strict $dimredConfig(sigma) ] } {
     tk_messageBox -icon error -type ok -title Message -message "Invalid value for sigma"
     return 0
  }
  if { ![ string is integer -strict $dimredConfig(aD) ] } {
     tk_messageBox -icon error -type ok -title Message -message "Invalid value for aD"
     return 0
  }
  if { ![ string is integer -strict $dimredConfig(ad) ] } {
     tk_messageBox -icon error -type ok -title Message -message "Invalid value for ad"
     return 0
  }
  if { ![ string is integer -strict $dimredConfig(bD) ] } {
     tk_messageBox -icon error -type ok -title Message -message "Invalid value for bD"
     return 0
  }
  if { ![ string is integer -strict $dimredConfig(bd) ] } {
     tk_messageBox -icon error -type ok -title Message -message "Invalid value for bd"
     return 0
  }

  # Now move the IMDS output file to the output for sketch-map 
  file rename -force IMDS_DATA SMAP_OUT

  # Run sketch-map proper
  for { set i 1 } { $i<=5 } { incr i } {
     # Delete output previous cycle
     catch { file delete SMAP_IN }
     # Move the output for the previous run of sketch-map to the input
     file rename -force SMAP_OUT SMAP_IN

     # Set the switching parameter
     set alpha [ expr 1.0 - $i*0.2 ]
     # Retrieve the width of the grid
     set gwidth [ expr  1.1 * [ retrieveDataWidth ] ]

     # Run sketch-map algorithm
     puts "Executing: dimred -D $ncv -d 2 $pflag -mode GMDS -fun xsigmoid -sigma $dimredConfig(sigma) -epxa $dimredConfig(aD) -lexpa $dimredConfig(ad) -expb $dimredConfig(bD) -lexpb $dimredConfig(bd) steps $dimredConfig(nconjgrad) -imode conjgrad $wflag -g1 $dimredConfig(ngridcoarse) -g2 $dimredConfig(ngridfine) -gw $gwidth -init SMAP_IN -imix $alpha < $ifile > SMAP_OUT"
     catch { exec $env(smapdir)/bin/dimred -D $ncv -d 2 $pflag -mode GMDS -fun xsigmoid -sigma $dimredConfig(sigma) -epxa $dimredConfig(aD) -lexpa $dimredConfig(ad) -expb $dimredConfig(bD) -lexpb $dimredConfig(bd) steps $dimredConfig(nconjgrad) -imode conjgrad $wflag -g1 $dimredConfig(ngridcoarse) -g2 $dimredConfig(ngridfine) -gw $gwidth -init SMAP_IN -imix $alpha < $ifile > SMAP_OUT } dimr_stdout

     # Read the file containing the output data
     if { [readMikFile -ifile SMAP_OUT -errors $dimr_stdout -npoints [molinfo top get numframes] -projname SMAP ]==0 } { return 0 }
    
     # Increment the progress bar
     ::progressbar::increment $progressBar
  }

  return 1
}

proc sketchmap::readMikFile {args} {
  variable dimredConfig

  set npoints [getargs $args "-npoints" "No number of points specified in call to readMikFile"]
  set ifile [getargs $args "-ifile" "No filename specified in call to readMikFile"]
  set errors [getargs $args "-errors" "No errors specified in call to readMikFile"]
  set dimredConfig(projname) [getargs $args "-projname" "No -projname specified in call to readMikFile"]

  # Check file has been generated
  if { ![file exists $ifile] } {
    puts $errors
    puts "-------------"
    tk_messageBox -icon error -type ok -title Message -message "Something went wrong when doing dimensionality reduction.  Please see messages in log"
    return 0
  }

  # Read in the file
  set dimredConfig(xx) {} ; set dimredConfig(yy) {}   ; set nproj 0
  set fd [open $ifile r]
  while {[gets $fd line] != -1} {
    set sline [regexp -inline -all -- {\S+} $line]
    if { [lindex $sline 0]=="#" } {
          continue
    } else {
        lappend dimredConfig(xx) [lindex $sline 0]
        lappend dimredConfig(yy) [lindex $sline 1]
        incr nproj
    }
  }
  close $fd

  # Check that the correct number of points is present
  if { $nproj!=$npoints } {
    tk_messageBox -icon error -type ok -title Message -message "Dimensionality reduction did not produce correct number of frames"
    return 0
  }

  return 1
}

proc sketchmap::retrieveDataWidth {args} {
  variable dimredConfig

  if { [llength $dimredConfig(xx)] != [ llength $dimredConfig(yy) ] } {
     tk_messageBox -icon error -type ok -title Message -message "Disaster: Number of x components does not equal number of y components"
     return 100.0
  }

  set min 0 ; set max 0
  for { set i 0 } { $i<[llength $dimredConfig(xx)] } { incr i } { 
      if { [lindex $dimredConfig(xx) $i]>$max } { set max [lindex $dimredConfig(xx) $i] } 
      if { [lindex $dimredConfig(yy) $i]>$max } { set max [lindex $dimredConfig(yy) $i] }
      if { [lindex $dimredConfig(xx) $i]<$min } { set min [lindex $dimredConfig(xx) $i] }
      if { [lindex $dimredConfig(yy) $i]<$min } { set min [lindex $dimredConfig(yy) $i] }
  }

  return [ expr ( $max - $min ) / 2 ]
}

proc sketchmap::setProgressBar {bar} {
   variable progressBar
   set progressBar $bar
}

proc sketchmap::projectionFileWindow {wind args} {
   variable pfileConfig

   frame $wind
   frame $wind.controls -padx 1m -pady 1m
   grid [ label $wind.controls.typelab -text "type" ] -row 1 -column 1 -columnspan 2
   grid [ entry $wind.controls.dtype -width 10 -textvariable sketchmap::pfileConfig(rtype) ] -row 1 -column 3 -columnspan 2
   grid [ label $wind.controls.sigmalab -text "sigma" ] -row 1 -column 5 -columnspan 2
   grid [ entry $wind.controls.sigma -width 10 -textvariable sketchmap::pfileConfig(sigma) ] -row 1 -column 7 -columnspan 2
   grid [ label $wind.controls.aDlab -text "aD" ] -row 2 -column 1
   grid [ entry $wind.controls.aD -width 5 -textvariable sketchmap::pfileConfig(aD) ] -row 2 -column 2
   grid [ label $wind.controls.bDlab -text "bD" ] -row 2 -column 3
   grid [ entry $wind.controls.bD -width 5 -textvariable sketchmap::pfileConfig(bD) ] -row 2 -column 4
   grid [ label $wind.controls.adlab -text "ad" ] -row 2 -column 5
   grid [ entry $wind.controls.ad -width 5 -textvariable sketchmap::pfileConfig(ad) ] -row 2 -column 6
   grid [ label $wind.controls.bdlab -text "bd" ] -row 2 -column 7
   grid [ entry $wind.controls.bd -width 5 -textvariable sketchmap::pfileConfig(bd) ] -row 2 -column 8
   pack $wind.controls -in $wind -fill both -expand 1
   return $wind
}

proc sketchmap::printProjectionData {args} {
   variable pfileConfig

   set od [getargs $args "-to" "No -to flag in call to printProjectionData"]

   if { $pfileConfig(rtype)=="mds" } {
      puts $od "LOW_D_FUNCTION TYPE distance"
      puts $od "HIGH_D_FUNCTION TYPE distance"
      return 1
   } elseif { $pfileConfig(rtype)=="smap" } {
      if { ![string is double -strict $pfileConfig(sigma) ] } {
         tk_messageBox -icon error -type ok -title Message -message "Invalid value for sigma"
         return 0
      }
      if { ![ string is integer -strict $pfileConfig(aD) ] } {
         tk_messageBox -icon error -type ok -title Message -message "Invalid value for aD"
         return 0
      }
      if { ![ string is integer -strict $pfileConfig(ad) ] } {
         tk_messageBox -icon error -type ok -title Message -message "Invalid value for ad"
         return 0
      }  
      if { ![ string is integer -strict $pfileConfig(bD) ] } {
         tk_messageBox -icon error -type ok -title Message -message "Invalid value for bD"
         return 0
      }
      if { ![ string is integer -strict $pfileConfig(bd) ] } {
         tk_messageBox -icon error -type ok -title Message -message "Invalid value for bd"
         return 0
      }
      puts $od "LOW_D_FUNCTION TYPE general SIGMA $pfileConfig(sigma) POWERS $pfileConfig(ad) $pfileConfig(bd)"
      puts $od "HIGH_D_FUNCTION TYPE general SIGMA $pfileConfig(sigma) POWERS $pfileConfig(aD) $pfileConfig(bD)"
      return 1
   } else {
      tk_messageBox -icon error -type ok -title Message -message "Invalid projection type"
      return 0
   }
}

proc sketchmap::outOfSampleWindow {wind args} {
   variable osampleConfig
   variable tmpdir
   global   env

   if { ![file exists "$env(smapdir)/bin/dimproj"] } {
     tk_messageBox -icon error -type ok -title Message -message "Could not find dimproj code.  You must define an smapdir inside your .bashrc file.  dimproj should then be in \$smapdir/bin/dimproj"
     return
   }
   set tmpdir [getargs $args "-tmpdir" "No tempory directory in call to outOfSampleWindow"]

  frame $wind
  frame $wind.controls -padx 1m -pady 1m
  grid [ button $wind.controls.but1 -text "read plumed input" -relief raised -command [namespace code sketchmap::readBespokeFile] ] -row 1 -column 1 -columnspan 4
  grid [ button $wind.controls.but2 -text "read dimred output" -relief raised -command [namespace code sketchmap::readDimredFile] ] -row 1 -column 5 -columnspan 4
  grid [ label $wind.controls.typelab -text "type" ] -row 2 -column 1 -columnspan 2
  grid [ entry $wind.controls.dtype -width 10 -textvariable sketchmap::osampleConfig(rtype) ] -row 2 -column 3 -columnspan 2
  grid [ label $wind.controls.sigmalab -text "sigma" ] -row 2 -column 5 -columnspan 2
  grid [ entry $wind.controls.sigma -width 10 -textvariable sketchmap::osampleConfig(sigma) ] -row 2 -column 7 -columnspan 2
  grid [ label $wind.controls.aDlab -text "aD" ] -row 3 -column 1
  grid [ entry $wind.controls.aD -width 5 -textvariable sketchmap::osampleConfig(aD) ] -row 3 -column 2
  grid [ label $wind.controls.bDlab -text "bD" ] -row 3 -column 3
  grid [ entry $wind.controls.bD -width 5 -textvariable sketchmap::osampleConfig(bD) ] -row 3 -column 4
  grid [ label $wind.controls.adlab -text "ad" ] -row 3 -column 5
  grid [ entry $wind.controls.ad -width 5 -textvariable sketchmap::osampleConfig(ad) ] -row 3 -column 6
  grid [ label $wind.controls.bdlab -text "bd" ] -row 3 -column 7
  grid [ entry $wind.controls.bd -width 5 -textvariable sketchmap::osampleConfig(bd) ] -row 3 -column 8
  pack $wind.controls -in $wind -fill both -expand 1
  # Optimization algorithm controls
  labelframe $wind.optc -relief ridge -bd 2 -text "optimization controls" -padx 2m -pady 2m
  set osampleConfig(period) 0
  grid [label $wind.optc.labe1 -text "periodic variables"] -row 1 -column 1
  grid [checkbutton $wind.optc.cb1 -variable sketchmap::osampleConfig(isperiodic) ] -row 1 -column 2
  grid [label $wind.optc.labe2 -text "period"] -row 2 -column 1
  grid [entry $wind.optc.entrr -width 5 -textvariable sketchmap::osampleConfig(period) ] -row 2 -column 2
  grid [label $wind.optc.lab0 -text "Use landmark weights"] -row 3 -column 1
  set osampleConfig(useweights) 1
  grid [checkbutton $wind.optc.b1 -variable sketchmap::osampleConfig(useweights) ] -row 3 -column 2
  grid [label $wind.optc.lab1 -text "(local optimizer) number of conjugate gradient steps"] -row 4 -column 1
  set osampleConfig(nconjgrad) 3
  grid [entry $wind.optc.ent1 -width 5 -textvariable sketchmap::osampleConfig(nconjgrad) ] -row 4 -column 2
  grid [label $wind.optc.lab2 -text "(global optimizer) number of points in the coarse grained grid" ] -row 5 -column 1
  set osampleConfig(ngridcoarse) 20
  grid [entry $wind.optc.ent2 -width 5 -textvariable sketchmap::osampleConfig(ngridcoarse) ] -row 5 -column 2
  grid [label $wind.optc.lab3 -text "(global optimizer) number of points in the fine grained grid" ] -row 6 -column 1
  set osampleConfig(ngridfine) 200
  grid [entry $wind.optc.ent3 -width 5 -textvariable sketchmap::osampleConfig(ngridfine) ] -row 6 -column 2
#  grid [label $wind.optc.lab4 -text "(global optimizer) amount of smearing to use in averaging" ] -row 7 -column 1
#  grid [entry $wind.optc.ent4 -width 5 -textvariable sketchmap::osampleConfig(smear) ] -row 7 -column 2
  pack $wind.optc -in $wind -fill both -expand 1

  return $wind
}

proc sketchmap::readDimredFiles {args} {
  variable tmpdir

  # Delete any old files containting the projection
  catch { eval file delete [glob -dir $tmpdir "LANDMARK_FILE"] }
  catch { eval file delete [glob -dir $tmpdir "PROJECTION_FILE"] }
  
  set osampleConfig(gotfiles) 1
  set highdimfile [ tk_getOpenFile -initialdir [pwd] -title "Select a file containing the high dimensionality landmarks" ]
  set lowdimfile [ tk_getOpenFile -initialdir [pwd] -title "Select a file containing the projections of the landmarks" ]
  file copy -force $highdimfile "$tmpdir/LANDMARK_FILE"
  file copy -force $lowdimfile "$tmpdir/PROJECTION_FILE"
}

proc sketchmap::readBespokeFile {args} {
  variable osampleConfig
  variable tmpdir
  
  set bfile [ tk_getOpenFile -initialdir [pwd] -title "Select a plumed input file for bespoke cvs" ]

  # Now read the bespoke file
  set fd [open $bfile r]   ;   set usingsmap 0
  set foundhF 0 ; set foundlF 0  ; set block "none"
  set lowdp {} ; set highdp {} ; set weights {}
  while { [ gets $fd line] !=-1 } {
     set sline [regexp -inline -all -- {\S+} $line]
     if { $block=="lowd" } {
         if { [lsearch $sline "LOW_D<"]>=0 } {
            set block "none"
         } else {
            lappend lowdp $line
         }
     } elseif { $block=="highd" } {
         if { [lsearch $sline "HIGH_D<"]>=0 } {
            set block "none"
         } else {
            lappend highdp $line
         }
     } elseif { $block=="weights" } {
         if { [lsearch $sline "WEIGHTS<"]>=0 } {
            set block "none"
         } else {
            lappend weights $line
         }
     } elseif { $block=="none" } {
        # Find start of low_D block
        if { [lsearch $sline "LOW_D>"]>=0 } {
             if { $block!="none" && [llength $lowdp]!=0 } {
                tk_messageBox -icon error -type ok -title Message -message "Bespoke file $bfile is nonsensical"
                return
             }
             set block "lowd"
        # Find start of high_D block
        } elseif { [lsearch $sline "HIGH_D>"]>=0 } {
            if { $block!="none" && [llength $highdp]!=0 } {
                tk_messageBox -icon error -type ok -title Message -message "Bespoke file $bfile is nonsensical"
                return
             }
             set block "highd"
        # Find start of weights block
        } elseif { [lsearch $sline "WEIGHTS>"]>=0 } {
             if { $block!="none" && [llength $weights]!=0 } {
                tk_messageBox -icon error -type ok -title Message -message "Bespoke file $bfile is nonsensical"
                return
             }
             set block "weights"
        # Find low dimensional function
        } elseif { [lsearch $sline "LOW_D_FUNCTION"]>=0 } {
           set foundlF 1
           set tpos [lsearch $sline "TYPE"]
           if { $tpos<0 } {
              tk_messageBox -icon error -type ok -title Message -message "Did not find TYPE for low dimensional function in bespoke file $bfile"
              return
           }
           set type [ lindex $sline [expr $tpos+1] ]
           if { $type=="distance" } {
               if { $usingsmap!=0 } {
                  tk_messageBox -icon error -type ok -title Message -message "Bespoke file $bfile is invalid. Cannot combine sigmoid function in one space with distances in other"
                  return
               }
               set osampleConfig(rtype) "mds"
           } elseif { $type=="compress" } {
               set usingsmap 1
               set osampleConfig(rtype) "smap"
               set osampleConfig(ad) 1
               set osampleConfig(bd) 1
           } elseif { $type=="sigmoid" } {
               set usingsmap 1
               set osampleConfig(rtype) "smap"
               set osampleConfig(ad) 2
               set osampleConfig(bd) 2
           } elseif { $type=="general" } {
               set usingsmap 1
               set osampleConfig(rtype) "smap"
               set spos [lsearch $sline "POWERS"]
               if { $tpos<0 } {
                  tk_messageBox -icon error -type ok -title Message -message "Did not find POWERS for low dimensional function in bespoke file $bfile"
                  return
               }
               set osampleConfig(ad) [ lindex $sline [expr $spos+1] ]
               set osampleConfig(bd) [ lindex $sline [expr $spos+2] ]
           } else {
              tk_messageBox -icon error -type ok -title Message -message "Found invalid function type for low dimensional function in bespoke file $bfile"
              return
           }
        # Search for the parameters of the high dimensional function
        } elseif { [lsearch $sline "HIGH_D_FUNCTION"]>=0 } {
           set foundhF 1
           set tpos [lsearch $sline "TYPE"]
           if { $tpos<0 } {
              tk_messageBox -icon error -type ok -title Message -message "Did not find TYPE for high dimensional function in bespoke file $bfile"
              return
           }
           set type [ lindex $sline [expr $tpos+1] ]
           if { $type=="distance" } {
               if { $usingsmap!=0 } {
                  tk_messageBox -icon error -type ok -title Message -message "Bespoke file $bfile is invalid. Cannot combine sigmoid function in one space with distances in other"
                  return
               }
               set osampleConfig(rtype) "mds"
           } elseif { $type=="compress" } {
               set usingsmap 1
               set osampleConfig(rtype) "smap"
               set osampleConfig(aD) 1
               set osampleConfig(bD) 1
           } elseif { $type=="sigmoid" } {
               set usingsmap 1
               set osampleConfig(rtype) "smap"
               set osampleConfig(aD) 2
               set osampleConfig(bD) 2
           } elseif { $type=="general" } {
               set usingsmap 1
               set osampleConfig(rtype) "smap"
               set spos [lsearch $sline "POWERS"]
               if { $tpos<0 } {
                  tk_messageBox -icon error -type ok -title Message -message "Did not find POWERS for high dimensional function in bespoke file $bfile"
                  return
               }
               set osampleConfig(aD) [ lindex $sline [expr $spos+1] ]
               set osampleConfig(bD) [ lindex $sline [expr $spos+2] ]
           } else {
              tk_messageBox -icon error -type ok -title Message -message "Found invalid function type for high dimensional function in bespoke file $bfile"
              return
           }
           # Find sigma
           if { $usingsmap==1 } {
                set sigpos [lsearch $sline "SIGMA"]
                if { $sigpos<0 } {
                     tk_messageBox -icon error -type ok -title Message -message "Did not find SIGMA parameter in bespoke file $bfile"
                     return
                }
                set osampleConfig(sigma) [ lindex $sline [expr $sigpos+1] ]
           }
        }
     }
  }
  close $fd

  # Checks on the input file
  if { [llength $lowdp]==0 } {
     tk_messageBox -icon error -type ok -title Message -message "Did not find any low dimensionality points in input file $bfile"
     return
  } elseif { [llength $highdp]==0 } {
     tk_messageBox -icon error -type ok -title Message -message "Did not find any high dimensionality points in input file $bfile"
     return
  } elseif { [llength $weights]==0 } {
     tk_messageBox -icon error -type ok -title Message -message "Did not find any weights in input file $bfile"
     return
  } elseif { [llength $lowdp]!=[llength $highdp] } {
     tk_messageBox -icon error -type ok -title Message -message "In file $bfile number of landmarks does not match number of projections"
     return
  } elseif { [llength $lowdp]!=[llength $weights] } {
     tk_messageBox -icon error -type ok -title Message -message "In file $bfile number of landmarks does not match number of weights"
     return
  } elseif { $foundlF!=1 } {
     tk_messageBox -icon error -type ok -title Message -message "Found no low dimensionality function defined in $bfile"
     return
  } elseif { $foundhF!=1 } {
     tk_messageBox -icon error -type ok -title Message -message "Found no high dimensionality function defined in $bfile"
     return
  }

  # Delete any old files containting the projection
  catch { eval file delete [glob -dir $tmpdir "LANDMARK_FILE"] }
  catch { eval file delete [glob -dir $tmpdir "PROJECTION_FILE"] }
  set ldf [open "$tmpdir/PROJECTION_FILE" w]  ;   set hdf [open "$tmpdir/LANDMARK_FILE" w]
  for { set i 0 } { $i<[llength $lowdp] } { incr i } {
     puts $ldf [lindex $lowdp $i]
     puts $hdf "[lindex $highdp $i] [lindex $weights $i]"
  } 
  close $ldf ; close $hdf
  # We always have weights when we use bespoke input
  set osampleConfig(useweights) 1
}

proc sketchmap::MDSWindow {wind args} {
   variable dimredConfig
   global   env

   if { ![file exists "$env(smapdir)/bin/dimred"] } { 
      tk_messageBox -icon error -type ok -title Message -message "Could not find dimred code.  You must define an smapdir inside your .bashrc file.  dimred should then be in \$smapdir/bin/dimred"
      return
   }

   frame $wind
   # check button to switch off all optimization
   labelframe $wind.cont -relief ridge -bd 2 -text "parameters" -padx 2m -pady 2m
   grid [label $wind.cont.lab -text "Skip iterative optimization"] -row 1 -column 1
   set dimredConfig(noiter) 0
   grid [ checkbutton $wind.cont.noit -variable sketchmap::dimredConfig(noiter) ] -row 1 -column 2
   pack $wind.cont -in $wind -fill both -expand 1
   # These are the controls for iterative optimization
   pack [ createOptimizationControl $wind.ocontrol ] -in $wind -fill both -expand 1
   return $wind
}

proc sketchmap::sketchMapWindow {wind args} {
   variable dimredConfig
   global   env
  
   if { ![file exists "$env(smapdir)/bin/dimred"] } { 
      tk_messageBox -icon error -type ok -title Message -message "Could not find dimred code.  You must define an smapdir inside your .bashrc file.  dimred should then be in \$smapdir/bin/dimred"
      return
   }

   frame $wind
   # The sketch map parameters        
   labelframe $wind.cont -relief ridge -bd 2 -text "parameters" -padx 2m -pady 2m
   grid [label $wind.cont.slab -text "sigma"] -row 1 -column 1
   grid [entry $wind.cont.sigma -width 5 -textvariable sketchmap::dimredConfig(sigma) ] -row 1 -column 2
   grid [label $wind.cont.aDlab -text "aD"] -row 1 -column 3
   grid [entry $wind.cont.aD -width 5 -textvariable sketchmap::dimredConfig(aD) ] -row 1 -column 4
   grid [label $wind.cont.bDlab -text "bD"] -row 1 -column 5
   grid [entry $wind.cont.bD -width 5 -textvariable sketchmap::dimredConfig(bD) ] -row 1 -column 6
   grid [label $wind.cont.rlab -text "switch rate"] -row 2 -column 1
   grid [entry $wind.cont.rate -width 5 -textvariable sketchmap::dimredConfig(srate) ] -row 2 -column 2
   grid [label $wind.cont.adlab -text "ad"] -row 2 -column 3
   grid [entry $wind.cont.ad -width 5 -textvariable sketchmap::dimredConfig(ad) ] -row 2 -column 4
   grid [label $wind.cont.bdlab -text "bd"] -row 2 -column 5
   grid [entry $wind.cont.bd -width 5 -textvariable sketchmap::dimredConfig(bd) ] -row 2 -column 6
   pack $wind.cont -in $wind -fill both -expand 1

   # These are the controls for iterative optimization
   pack [ createOptimizationControl $wind.ocontrol ] -in $wind -fill both -expand 1

   return $wind
}

proc sketchmap::createOptimizationControl {wind args} {
  variable dimredConfig

  labelframe $wind -relief ridge -bd 2 -text "optimization controls" -padx 2m -pady 2m
  set dimredConfig(period) 0
  grid [label $wind.labe1 -text "periodic variables"] -row 1 -column 1
  grid [checkbutton $wind.cb1 -variable sketchmap::dimredConfig(isperiodic) ] -row 1 -column 2
  grid [label $wind.labe2 -text "period"] -row 2 -column 1
  grid [entry $wind.entrr -width 5 -textvariable sketchmap::dimredConfig(period) ] -row 2 -column 2
  grid [label $wind.lab0 -text "Use landmark weights"] -row 3 -column 1
  set dimredConfig(useweights) 1
  grid [checkbutton $wind.b1 -variable sketchmap::dimredConfig(useweights) ] -row 3 -column 2
  grid [label $wind.lab1 -text "(local optimizer) number of conjugate gradient steps"] -row 4 -column 1
  set dimredConfig(nconjgrad) 10
  grid [entry $wind.ent1 -width 5 -textvariable sketchmap::dimredConfig(nconjgrad) ] -row 4 -column 2
  grid [label $wind.lab2 -text "(global optimizer) number of points in the coarse grained grid" ] -row 5 -column 1
  set dimredConfig(ngridcoarse) 20
  grid [entry $wind.ent2 -width 5 -textvariable sketchmap::dimredConfig(ngridcoarse) ] -row 5 -column 2
  grid [label $wind.lab3 -text "(global optimizer) number of points in the fine grained grid" ] -row 6 -column 1
  set dimredConfig(ngridfine) 200
  grid [entry $wind.ent3 -width 5 -textvariable sketchmap::dimredConfig(ngridfine) ] -row 6 -column 2
  return $wind
}

proc sketchmap::getargs {arglist tag {warning 0} } {
    set pos [lsearch $arglist $tag]
    if { $pos<0 && $warning!=0 } {
       tk_messageBox -icon error -type ok -title Message -message "$warning"
       return
    } elseif { $pos<0 } {
       set val 0
    } else {
       set val [ lindex $arglist [expr $pos+1 ] ]
    }
    return $val
}


