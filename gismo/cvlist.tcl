#
# Tools for manipulating lists of colvars 
#

source $env(VMDDIR)/scripts/init.d/progress_bar.tcl
source $env(VMDDIR)/scripts/init.d/sketchmap.tcl

package provide cvlist 1.0
package require progressbar 1.0
package require sketchmap 1.0

namespace eval ::cvlist {
   namespace export readColvar        ; # Read in a colvar file
   namespace export exists            ; # Checks to see if there is cvdata with a given tag
   namespace export updateMenu        ; # Update a menu that contains a list of all the colvar names
   namespace export getDataCol        ; # Returns a list containing a list of cvdata
   namespace export getDataRow        ; # Returns a list containing an instant of cvdata 
   namespace export addDataCol        ; # Add a new column of data to a dataset
   namespace export clone             ; # Clones a list of frames from a cvlist to a cvlist with a new tag
   namespace export resetMenus        ; # Reset the menus and destroy anything like fps
   namespace export deleteTopMolData  ; # Delete all cvdata for the top molecule
# Stuff for interface to michelino's dimlandmarks code
   namespace export landmarkConfig   
   namespace export colvarSelection
   namespace export landmarkWindow     ; # Creates a control window for running landmark code
   namespace export createLandmarks    ; # Actually run dimlandmarks on the data 
   namespace export getLandmarkWeights ; # Return the list of landmark weights
# Stuff for dimensionality reduction
   namespace export dimredFinished
   namespace export dimensionalityReductionWindow
# Stuff for out of sample embedding
   namespace export osampleFinished
   namespace export useDataFromMol
   namespace export outOfSampleWindow
   namespace export destroyOutOfSampleWindow
   namespace export outOfSampleEmbedding
# Stuff for outputting projection files
   namespace export projectionFileFinished 
   namespace export projectionFileWindow
   namespace export destroyProjectionFileWindow
   namespace export printProjectionFile
# The variables   
   variable osample_blocksize 100   ; # The number of out of sample points to do at a time
   variable data                    ; # this will hold all the cvdata
   variable landmarkConfig          ; # The configuration we are using to create landmark data points
   variable dimredFinished          ; # Tells plumed_gui that dimensionality reduction is completed
   variable osampleFinished         ; # Tells plumed_gui that out of sample embedding is completed
   variable useDataFromMol          ; # Used for out of sample embedding
   variable osampleMenu             ; # The dynamic menu for out of sample embedding 
   variable osampleProjNumber       ; # The particular projection of the data we are using
   variable colvarSelection         ; # A selection of colvars to use for something or other
   variable dataOperationTag        ; # A tag that remembers what data we need to operate something on
   variable tmpdir                  ; # A tempory directory
   variable osampleTabs             ; # Holds all the tabs for the osampleWindow
   variable osampleProgress         ; # Holds the location of the osample progress bar
   variable projectionFileFinished  ; # Tells plumed_gui that projection file printing is finished
   variable pfileConfig             ; # Holds everything for creating a projection file
   variable pfileTabs               ; # Holds all the tabs for the printProjectionFileWindow
   variable pfileProjNumber         ; # The particular projection of the data we want to print out
   variable pfileMol                ; # The molecule for which we are creating a projection file
}

proc cvlist::projectionFileWindow {wind args} {
  variable dataOperationTag
  variable pfileConfig
  variable pfileTabs
  variable pfileProjNumber
  variable pfileMol
  variable data

  set dtag [getargs $args "-datatag" "No -datatag flag in call to projectionFileWindow"]
  set pfileMol $dtag:[findMolecule]

  frame $wind -padx 1m -pady 1m 

  # Create a list of dimensionality reductions that have been performed on this molecule
  frame $wind.dimred -padx 2m -pady 2m
  grid [label $wind.dimred.l -text "Dimensionality reduction data available for mol: [findMolecule]" ] -row 1 -column 1
  grid [menubutton $wind.dimred.mbut -relief raised -bd 2 -direction flush -width 10 -menu $wind.dimred.mbut.menu ] -row 1 -column 2
  menu $wind.dimred.mbut.menu
  $wind.dimred.mbut.menu delete 0 end
  # Find the parameter sets for the particular molecule and set up commands so that we can retrieve data
  set mtag "$dtag:[findMolecule]" 
  foreach { p1 dat } [array get data] {
     if { [lindex $p1 0]==$mtag && [lindex $p1 3]=="name" } {
        set number [lindex $p1 2]
        $wind.dimred.mbut.menu add radiobutton -label "$dat" -value "$number" -variable cvlist::pfileProjNumber
     }
  }
  pack $wind.dimred -in $wind -fill both -expand 1
 
  # Coordinates in the high and low dimensional spaces
  pack [ createColvarSelector $wind.highdc -datatag $dtag -selectiontag "highd" -text "Coordinates for D dimensional space" -ncols 6 ] -in $wind -side top -fill both
  pack [ createColvarSelector $wind.lowdc -datatag $dtag -selectiontag "lowd" -text "Coordinates for d dimensional space" -ncols 6 ] -in $wind -side top -fill both

  # Something for including weights
  frame $wind.wfram -padx 3m -pady 3m
  grid [label $wind.wfram.l -text "Weights are uniform" ] -row 1 -column 1
  set pfileConfig(uniformWeights) 0
  grid [checkbutton $wind.wfram.cb -variable cvlist::pfileConfig(uniformWeights)] -row 1 -column 2
  pack $wind.wfram -in $wind -fill both -expand 1 

  # Method specific stuff
  pack [ttk::notebook $wind.nb] -fill both -expand 1
  set pfileTabs(notebook) $wind.nb   
  # Stuff for sketch-map tab
  $wind.nb add [ ::sketchmap::projectionFileWindow $wind.nb.smap ] -text "sketch-map"
  set pfileTabs(sketch-map) $wind.nb.smap   

  # And finish the notebook
  ttk::notebook::enableTraversal $wind.nb
  $wind.nb select $wind.nb.smap

  trace add variable cvlist::pfileProjNumber write [namespace code cvlist::updatePfileParameters]
  return $wind
}

proc cvlist::destroyProjectionFileWindow {args} {
  trace remove variable cvlist::pfileProjNumber write [namespace code cvlist::updatePfileParameters]
}

proc cvlist::printProjectionFile {args} {
  variable projectionFileFinished
  variable pfileConfig
  variable pfileTabs
  variable data

  # Retrieve the data we are printing
  set dtag [getargs $args "-datatag" "No -datatag flag in call to printProjectionFile"]
  # Get the name of the file we are printing to
  set filen [tk_getSaveFile -initialfile "bespoke.defs" -title "Print projection data to file" -parent . \
    -filetypes [list {{Postscript files} {.defs}} {{All files} {*}}]]
  if { $filen == "" } { return }

  # Open a file for output
  set od [open $filen w]

  # Print number of landmarks
  set molno [findMolecule]
  puts $od "NLANDMARKS [molinfo $molno get numframes]" 
  puts $od " "

  # Find which tab is curretly on top
  set tabname [$pfileTabs(notebook) tab current -text]

  # Print method specific data
  if { $tabname=="sketch-map" } {
     if { [ ::sketchmap::printProjectionData -to $od ]==0 } { close $od ; file delete -force $filen ; return }
  } else {
     puts "This method ($tabname) has no method for producing projection files"
  }
  puts $od " "

  # Print high dimensional collective coordinates
  puts $od "HIGH_D>"
  if { [printSelectedColvars -to $od -datatag $dtag -all -selectiontag "highd"]==0 } { close $od ; file delete -force $filen ; return }
  puts $od "HIGH_D<"
  puts $od " "
  # Print low dimensional collective coordinates
  puts $od "LOW_D>"
  if { [printSelectedColvars -to $od -datatag $dtag -all -selectiontag "lowd"]==0 } { close $od ; file delete -force $filen ; return }
  puts $od "LOW_D<"
  puts $od " "
  # Check if we have weights
  set areweights 0
  set pgtag "$dtag:[findMolecule]"
  foreach name $data([list $pgtag cv_names]) {
     if { $name=="weights" } { set areweights 1 }
  }
  if { $areweights==0 } { set pfileConfig(uniformWeights) 1 }

  # Print weights
  puts $od "WEIGHTS>"
  if { $pfileConfig(uniformWeights)==1 } { 
     for { set i 0 } { $i<[molinfo $molno get numframes] } { incr i } { puts $od "1.0" }
  } else {
     selectCollectiveVariables -moltag $dtag:[findMolecule] -variables [list "weights"] -selectiontag "weights"
     if { [printSelectedColvars -to $od -datatag $dtag -all -selectiontag "weights"]!=1 } { 
        tk_messageBox -icon error -type ok -title Message -message "Found more than 1 column of weights in data"
        close $od ; file delete -force $filen  
        return
     }
  }
  puts $od "WEIGHTS<"

  # Close the file
  close $od

  # This destroys the projection file window
  set projectionFileFinished 0
}

proc cvlist::updatePfileParameters {args} {
  variable pfileProjNumber
  variable pfileMol
  variable pfileTabs
  variable pfileConfig
  variable colvarSelection
  variable data

  if { ![string is integer -strict $pfileProjNumber] } { return } 

  # Now select the cvs
  selectCollectiveVariables -moltag $pfileMol -variables $data([list $pfileMol dimred $pfileProjNumber cvlist]) -selectiontag "highd"

  # Check if weights is selected, unselect it and turn off uniform weights
  set pfileConfig(uniformWeights) 1
  foreach name $data([list $pfileMol cv_names]) {
     if { $name=="weights" } { set pfileConfig(uniformWeights) 0 }
  }
  if { $pfileConfig(uniformWeights)==0 } {
     if { $colvarSelection([list $pfileMol weights highd])==1 } { 
        set colvarSelection([list $pfileMol weights highd]) 0
        set pfileConfig(uniformWeights) 0
     } else {
        set pfileConfig(uniformWeights) 1
     }
  }

  # Select the projection coordinates
  set projselect {}
  set projname $data([list $pfileMol dimred $pfileProjNumber name])
  foreach cvthere $data([list $pfileMol cv_names]) {
     if { [string first $projname $cvthere] != -1 } { lappend projselect $cvthere }
  }
  # Print the projection coordinates
  selectCollectiveVariables -moltag $pfileMol -variables $projselect -selectiontag "lowd"

  # Get the type and retrieve the parameters
  set tno [ lsearch $data([list $pfileMol dimred $pfileProjNumber parameters]) "TYPE" ]
  set type [ lindex $data([list $pfileMol dimred $pfileProjNumber parameters]) [expr $tno + 1] ]
  if { $type=="smap" || $type=="mds" } {
      ::sketchmap::setPfileParameters -parameters $data([list $pfileMol dimred $pfileProjNumber parameters])
  } else {
     puts "There is no way of setting the projection file parameters for dimensionality reduction method: $type"
     return
  }

  # Switch the notebook to the appropriate tab
  if { $type=="mds" || $type=="smap" } {
     $pfileTabs(notebook) select $pfileTabs(sketch-map)
  } else {
     $pfileTabs(notebook) select $pfileTabs($type)
  }
}

proc cvlist::createLandmarks {args} {
  variable data
  variable landmarkConfig
  global   env

  set dtag [getargs $args "-datatag" "No -datatag flag in call to createLandmarks"]
  set tag $dtag:[findMolecule]

  if { ![info exists landmarkConfig([list $tag tmpdir])] } {
     tk_messageBox -icon error -type ok -title Message -message "Something has gone wrong in creating landmarks for this data set"
     return
  }
  # Check that user has specified a number of landmarks
  if { ![string is integer -strict $landmarkConfig([list $tag npoints]) ] } {
     tk_messageBox -icon error -type ok -title Message -message "Please specify the number of landmarks that you require"
     return
  }

  # Change to the working directory
  set odir [pwd] ; cd $landmarkConfig([list $tag tmpdir])
  puts "tempory directory is $landmarkConfig([list $tag tmpdir])"
  set od [open "$landmarkConfig([list $tag tmpdir])/CV_DATA" w]
  set ncv [printSelectedColvars -to $od -datatag $dtag -all -selectiontag "print"]
  if { $ncv==0 } { close $od ; cd $odir ; return }
  close $od  

  # And the periodicity stuff
  if { $landmarkConfig([list $tag isperiodic])==0 } { 
      set pflag {}
  } else {
      set pflag "-pi $landmarkConfig([list $tag period])"
  } 

  puts "Execute: dimlandmark -D $ncv -n $landmarkConfig([list $tag npoints]) $pflag -mode $landmarkConfig([list $tag mode]) -i -w -lowmem < CV_DATA > landmarks"
  catch { exec $env(smapdir)/bin/dimlandmark -D $ncv -n $landmarkConfig([list $tag npoints]) $pflag -mode $landmarkConfig([list $tag mode]) -i -w -lowmem < CV_DATA > landmarks } diml_stdout

  if { [file exists landmarks] } {
     set landmarks {}  ; set fd [open landmarks r] 
     set landmarkConfig([list $tag weights]) {}
     while {[gets $fd line] != -1} {
       set sline [regexp -inline -all -- {\S+} $line]
       if { [lindex $sline 0]=="#" } { 
            continue 
       } elseif { [llength $sline]!=[expr $ncv+2] } {
            puts $diml_stdout
            puts "-------------" 
            tk_messageBox -icon error -type ok -title Message -message "Found dubious line in landmarks file.  Please see messages in log"  
            return
       } else {
            lappend landmarks [lindex $sline 0]
            lappend landmarkConfig([list $tag weights]) [lindex $sline $ncv+1]
       }
     }
     if { [llength $landmarks]!=$landmarkConfig([list $tag npoints]) } {
        puts $diml_stdout
        puts "-------------" 
        tk_messageBox -icon error -type ok -title Message -message "Wrong number of landmarks generated.  Please see messages in log"
        return
     }
  } else {
     puts $diml_stdout
     puts "-------------"
     tk_messageBox -icon error -type ok -title Message -message "Something went wrong when running dimlandmarks.  Please see messages in log"
  }

  cd $odir
  return $landmarks
}

proc cvlist::getLandmarkWeights {args} {
   variable landmarkConfig
   set dtag [getargs $args "-datatag" "No -datatag flag in call to getLandmarkWeights"]
   set tag $dtag:[findMolecule]
   return $landmarkConfig([list $tag weights])
}

proc cvlist::outOfSampleEmbedding {args} {
   variable tmpdir
   variable data
   variable osampleProjNumber
   variable useDataFromMol
   variable osampleFinished
   variable osampleProgress
   variable osample_blocksize
   variable osampleTabs

   if { ![info exists tmpdir] } {
     tk_messageBox -icon error -type ok -title Message -message "Something has gone wrong during out of sample"
     return
   }

   # We must configure the progress bar
   set nblocks [ expr floor( [molinfo top get numframes ] / $osample_blocksize ) ]
   ::progressbar::config $osampleProgress -ntasks [ expr $nblocks + 1 ]
   
   set dtag [getargs $args "-datatag" "No -datatag flag in call to outOfSampleEmbedding"]

   # Get the number of colvars
   set tag $dtag:[findMolecule]
   # set ncv [ llength $data([list $tag sub_names]) ] 

   # Change to the tempory directory
   set odir [pwd] ; cd $tmpdir

   # Get the number of landmarks
   set fd [open "LANDMARK_FILE" r]
   set nlandmarks 0
   while { [gets $fd line] !=-1 } { incr nlandmarks }
   close $fd

   # Now do each block of out of sample data
   set projection(0) {}   ;  set projection(1) {}
   for { set iblock 0 } { $iblock<=$nblocks } { incr iblock } {
      if { $iblock==$nblocks } { 
         set blockend  [molinfo top get numframes ] 
         set blocksize [ expr [molinfo top get numframes ] - $nblocks*$osample_blocksize ]
      } else {
         set blockend [ expr ($iblock+1)*$osample_blocksize ]
         set blocksize $osample_blocksize
      }

      # Print out the colvars
      set od [open "$tmpdir/CV_DATA" w]
      set ncv [printSelectedColvars -to $od -datatag $dtag -start [expr $iblock*$osample_blocksize] -end $blockend -selectiontag "print"]
      if { $ncv==0 } { return }
      close $od

      set tabname [$osampleTabs(notebook) tab current -text]
      if { $tabname=="sketch-map" } {
         if { [::sketchmap::runOutOfSample -ncv $ncv -nlandmarks $nlandmarks -npoints $blocksize -ifile "CV_DATA" -lhfile "LANDMARK_FILE" -llfile "PROJECTION_FILE" ]==0 } { 
              cd $odir ; return 
         }
         set nproj 2
         if { [llength $::sketchmap::dimredConfig(xx)] != [llength $::sketchmap::dimredConfig(yy)] } { 
              puts "!ERROR: something has gone wrong in retrieving out of sample data"
         }

         for { set j 0 } { $j<[llength $::sketchmap::dimredConfig(xx)] } { incr j } {
             lappend projection(0) [lindex $::sketchmap::dimredConfig(xx) $j]
             lappend projection(1) [lindex $::sketchmap::dimredConfig(yy) $j]
         }
      } else {
         puts "This method ($tabname) has no out of sample projection method"
      }
      ::progressbar::increment $osampleProgress
   }
   
   cd $odir   ; unset tmpdir

   # Get a name for this projection
   if { $osampleProjNumber!="" } {
      set projname $data([list $useDataFromMol dimred $osampleProjNumber name])
   } else {
      set projname "user"
   }

   # Add output from dimensionality reduction to data for this molecule
   for { set i 0  } { $i<$nproj } { incr i } {
     addDataCol -datatag $dtag -name "proj-$projname-$i" -data $projection($i)
   }
   # Setting this variable destroys the window
   set osampleFinished 0
}

proc cvlist::outOfSampleWindow {wind args} {
  variable dataOperationTag
  variable osampleTabs
  variable osampleMenu
  variable osampleProgress
  variable osampleProjNumber
  variable tmpdir

  set dtag [getargs $args "-datatag" "No -datatag flag in call to outOfSampleWindow"]
  set dataOperationTag $dtag:[findMolecule]
  set tmpdir [getargs $args "-tmpd" "No -tmpd flag in call to outOfSampleWindow"]
#  set cvnames [getargs $args "-name" "No -name flag in call to outOfSampleWindow"]

  frame $wind -padx 1m -pady 1m
  pack [ createColvarSelector $wind.sel -datatag $dtag -selectiontag "print" -text "select the cvs to use" -ncols 6 ] -in $wind -side top -fill both

  frame $wind.controls -padx 2m -pady 2m
  grid [ label $wind.controls.lmol -text "Use molecule" ] -row 1 -column 1 
  grid [ menubutton $wind.controls.mmol -relief raised -bd 2 -direction flush -width 10 -menu $wind.controls.mmol.menu ] -row 1 -column 2
  menu $wind.controls.mmol.menu
  foreach m [molinfo list] {
     $wind.controls.mmol.menu add radiobutton -label "molecule $m" -value "$dtag:$m" -variable cvlist::useDataFromMol
  }
  grid [ label $wind.controls.lred -text "Data" ] -row 1 -column 3
  grid [ menubutton $wind.controls.mred -relief raised -bd 2 -direction flush -width 10 -menu $wind.controls.mred.menu ] -row 1 -column 4     
  menu $wind.controls.mred.menu
  set osampleMenu $wind.controls.mred.menu
  pack $wind.controls -in $wind -fill both -expand 1

  pack [ttk::notebook $wind.nb] -fill both -expand 1
  set osampleTabs(notebook) $wind.nb
  # Stuff for sketch-map tab
  $wind.nb add [ ::sketchmap::outOfSampleWindow $wind.nb.smap -tmpdir $tmpdir ] -text "sketch-map"
  set osampleTabs(sketch-map) $wind.nb.smap

  # And finish the notebook
  ttk::notebook::enableTraversal $wind.nb
  $wind.nb select $wind.nb.smap

  # This is a progress bar
  frame $wind.bar -padx 3m -pady 3m
  pack [ label $wind.bar.lab -text "progress" ] -in $wind.bar -side left
  pack [ ::progressbar::create $wind.bar.progress -width 360 -height 10 ] -in $wind.bar -side right
  set osampleProgress $wind.bar.progress
  pack $wind.bar -in $wind -side top

  trace add variable cvlist::useDataFromMol write [namespace code cvlist::updateOSampleMenu]
  trace add variable cvlist::osampleProjNumber write [namespace code cvlist::updateParameters]

  return $wind
}

proc cvlist::destroyOutOfSampleWindow {args} {
  trace remove variable cvlist::useDataFromMol write [namespace code cvlist::updateOSampleMenu]
  trace remove variable cvlist::osampleProjNumber write [namespace code cvlist::updateParameters]
}

proc cvlist::updateOSampleMenu {args} {
  variable data
  variable useDataFromMol
  variable osampleProjNumber
  variable osampleMenu

  set mtag $useDataFromMol   ; # mtag is the molecule we have previously embedded
 
  $osampleMenu delete 0 end
  # Find the parameter sets for the particular molecule and set up commands so that we can retrieve data
  foreach { p1 dat } [array get data] {
     if { [lindex $p1 0]==$mtag && [lindex $p1 3]=="name" } {
        set number [lindex $p1 2]
        $osampleMenu add radiobutton -label "$dat" -value "$number" -variable cvlist::osampleProjNumber
     } 
  }
}

proc cvlist::updateParameters {args} {
# This routine sets the parameters based on an old parameter set that dimensionality reduction was performed on
  variable osampleProjNumber
  variable useDataFromMol
  variable osampleTabs
  variable colvarSelection
  variable data
  variable tmpdir

  if { ![string is integer -strict $osampleProjNumber ] } { return }
#  puts "In updateParameters molecule is $useDataFromMol projection data set is $osampleProjNumber"
#  puts "The tempory director is $tmpdir" 

  # Delete any old files containting the projection
  catch { eval file delete [glob -dir $tmpdir "LANDMARK_FILE"] }
  catch { eval file delete [glob -dir $tmpdir "PROJECTION_FILE"] }

  # Select the projection coordinates 
  set projselect {}
  set projname $data([list $useDataFromMol dimred $osampleProjNumber name])
  foreach cvthere $data([list $useDataFromMol cv_names]) {
     if { [string first $projname $cvthere] != -1 } { lappend projselect $cvthere }
  }
  # Print the projection coordinates
  selectCollectiveVariables -moltag $useDataFromMol -variables $projselect -selectiontag "print"
  set od1 [open "$tmpdir/PROJECTION_FILE" w]
  if { [printSelectedColvars -to $od1 -datatag [lindex [split $useDataFromMol :] 0] -molecule [lindex [split $useDataFromMol :] 1] -all -selectiontag "print"]==0} { return }
  close $od1 

  # Now select the cvs
  selectCollectiveVariables -moltag $useDataFromMol -variables $data([list $useDataFromMol dimred $osampleProjNumber cvlist]) -selectiontag "print" 
  set od2 [open "$tmpdir/LANDMARK_FILE" w]
  if { [printSelectedColvars -to $od2 -datatag [lindex [split $useDataFromMol :] 0] -molecule [lindex [split $useDataFromMol :] 1] -all -selectiontag "print"]==0 } { return }
  close $od2

  # Now select the cvs in the top molecule -- these are the ones we will print at the start of outOfSampleEmbedding
  foreach cvselect $data([list $useDataFromMol dimred $osampleProjNumber cvlist]) {
     if { $cvselect!="weights" } { lappend govariables $cvselect ; }
  }
  set topmoltag [lindex [split $useDataFromMol :] 0]:[findMolecule]
  selectCollectiveVariables -moltag $topmoltag -variables $govariables -selectiontag "print"

  # Find out whether or not we are using weights
  set useweights 0
  foreach cvthere $data([list $useDataFromMol cv_names]) { 
     if { $cvthere=="weights" && $colvarSelection([list $useDataFromMol $cvthere print])==1 } { set useweights 1 }
  }

  # Get the type and retrieve the parameters
  set tno [ lsearch $data([list $useDataFromMol dimred $osampleProjNumber parameters]) "TYPE" ] 
  set type [ lindex $data([list $useDataFromMol dimred $osampleProjNumber parameters]) [expr $tno + 1] ] 
  if { $type=="smap" || $type=="mds" } {
      ::sketchmap::setParameters -parameters $data([list $useDataFromMol dimred $osampleProjNumber parameters])
      set ::sketchmap::osampleConfig(useweights) $useweights
  } else {
     puts "This method ($type) is not implemented inside the dimensionality reduction guff"
     return
  }
  # Switch the notebook to the appropriate tab
  if { $type=="mds" || $type=="smap" } { 
     $osampleTabs(notebook) select $osampleTabs(sketch-map)
  } else {
     $osampleTabs(notebook) select $osampleTabs($type)
  }
}

proc cvlist::dimensionalityReductionWindow {wind args} {
   variable dataOperationTag
   variable tmpdir

   set dtag [getargs $args "-datatag" "No -datatag flag in call to landmarkWindow"]
   set dataOperationTag $dtag:[findMolecule]
   set tmpdir [getargs $args "-tmpd" "No -tmpd flag in call to landmarkWindow"]

   frame $wind -padx 1m -pady 1m
   pack [ createColvarSelector $wind.sel -datatag $dtag -selectiontag "print" -text "select the cvs to use" -ncols 6 ] -in $wind -side top -fill both

   pack [ttk::notebook $wind.nb] -fill both -expand 1
   $wind.nb add [ ::sketchmap::MDSWindow $wind.nb.mds ] -text "MDS"
   # Buttons to make it go
   frame $wind.nb.mds.buttons -padx 1m -pady 1m
   pack [ button $wind.nb.mds.buttons.cancel -text "cancel" -relief raised -command {set cvlist::dimredFinished 0} ] -in $wind.nb.mds.buttons -side left
   pack [ button $wind.nb.mds.buttons.ok -text "run" -relief raised -command [namespace code {cvlist::dimensionalityReduction -type mds}] ] -in $wind.nb.mds.buttons -side right
   pack $wind.nb.mds.buttons -in $wind.nb.mds -fill both -expand 1
   $wind.nb add [ ::sketchmap::sketchMapWindow $wind.nb.smap ] -text "sketch-map"
   # Buttons to make it go
   frame $wind.nb.smap.buttons -padx 1m -pady 1m
   pack [ button $wind.nb.smap.buttons.cancel -text "cancel" -relief raised -command {set cvlist::dimredFinished 0} ] -in $wind.nb.smap.buttons -side left
   pack [ button $wind.nb.smap.buttons.ok -text "run" -relief raised -command [namespace code {cvlist::dimensionalityReduction -type smap}] ] -in $wind.nb.smap.buttons -side right
   pack $wind.nb.smap.buttons -in $wind.nb.smap -fill both -expand 1
   $wind.nb select $wind.nb.smap
   ttk::notebook::enableTraversal $wind.nb

   # We need our status bar here
   frame $wind.bar -padx 3m -pady 3m
   pack [ label $wind.bar.lab -text "progress" ] -in $wind.bar -side left
   pack [ ::progressbar::create $wind.bar.progress -width 360 -height 10 ] -in $wind.bar -side right
   ::sketchmap::setProgressBar $wind.bar.progress
   pack $wind.bar -in $wind -side top

   return $wind
}

proc cvlist::dimensionalityReduction {args} {
  variable dimredFinished
  variable dataOperationTag
  variable tmpdir
  variable colvarSelection
  variable data

  set tag $dataOperationTag
  set type [getargs $args "-type" "No -type flag in call to createDimRedGoButtons"]

  if { ![info exists tmpdir] } {
     tk_messageBox -icon error -type ok -title Message -message "Something has gone wrong during dimensionality reduction, namely tmpdir"
     return
  }

  # Check that some cvs have been selected
  set nselected 0
  foreach name $data([list $tag cv_names]) {
    if { $colvarSelection([list $tag $name print])==1 } { incr nselected }
  }
  if { $nselected==0 } {
     tk_messageBox -icon error -type ok -title Message -message "You have not made a collective variable selection"
     return
  }

  if { [molinfo top get numframes]>1000 } {
     set ans [tk_messageBox -icon question -type yesno -title Message \
              -message "There are more than 1000 points in the trajectory.  Do you really want to do dimensionality reduction - it will take ages?"]
     if { $ans=="no" } { return }
  }

  if { $type=="smap" || $type=="mds" } {
     # Find out if we need weights data
     set useweights  $::sketchmap::dimredConfig(useweights) 
     # Find the number of coordinates we are reducing to
     set nproj 2
  } else {
     puts "This method ($type) is not implemented inside the dimensionality reduction guff"
     return
  }
   
  set foundweights 0
  if { $useweights==1 } {
     foreach name $data([list $tag cv_names ]) {
       if { $name=="weights" } { set foundweights 1; set colvarSelection([list $tag $name print]) 1 }
     }
  }

  # Store the cvs used for this particular molecule
  set cvselection {}
  foreach name $data([list $tag cv_names]) {
    if { $colvarSelection([list $tag $name print])==1 } { lappend cvselection $name }
  }

  # Output the collective variables we need ready for smap projection
  set odir [pwd] ; cd $tmpdir
  puts "tempory directory is $tmpdir"
  set dtag [lindex [split $tag :] 0]
  set od [open "$tmpdir/CV_DATA" w]
  set ncv [expr [printSelectedColvars -to $od -datatag $dtag -all -selectiontag "print"] - $foundweights]
  if { $ncv==0 } { cd $odir; close $od ; return }
  close $od

  # Get the number of colvars
  #set ncv [expr [ llength $data([list $tag sub_names]) ] - $foundweights]

  puts "Printed out $ncv collective variables for dimensionality reduction and foundweights is $foundweights" 

  # Run the dimensionalilty reduction algorithm
  if { $type=="smap" || $type=="mds" } { 
     if { [::sketchmap::runSketchMap -type $type -ncv $ncv -ifile "CV_DATA"]==0 } { cd $odir; return }
     set projection(0) $::sketchmap::dimredConfig(xx)
     set projection(1) $::sketchmap::dimredConfig(yy)
  } else {
     puts "This method ($type) is not implemented inside the dimensionality reduction guff"
     return
  } 
  cd $odir        ;   unset tmpdir

  # Come up with a number for this particular dimensionality reduction
  if { ![info exists data([list $tag ndimred])] } {
      set data([list $tag ndimred]) 0
  } else {
      incr data([list $tag ndimred])
  }
  set ndimred $data([list $tag ndimred])
  # Store the name and type of the dimensionality reduction
  set data([list $tag dimred $ndimred name]) $type:$ndimred

  # Add output from dimensionality reduction to data for this molecule
  for { set i 0  } { $i<$nproj } { incr i } {
     addDataCol -datatag $dtag -name "$type:$data([list $tag ndimred])-$i" -data $projection($i)
  }

  # Store the parameters of the dimensionality reduction
  if { $type=="smap" || $type=="mds" } {
     set data([list $tag dimred $ndimred parameters]) [ ::sketchmap::retrieveParams -type $type ]
  } else {
    puts "There is no method implemented for storing parameteters of dimensionality reduction type: $type" 
  }

  # Store the cvs used for this particular molecule
  set data([list $tag dimred $ndimred cvlist]) $cvselection

  # Setting this should destroy the dimensionality reduction window
  set dimredFinished 0
}

proc cvlist::landmarkWindow {wind args} {
   variable landmarkConfig
   global   env

   set dtag [getargs $args "-datatag" "No -datatag flag in call to landmarkWindow"]
   set tag $dtag:[findMolecule]
   set landmarkConfig([list $tag tmpdir]) [getargs $args "-tmpd" "No -tmpd flag in call to landmarkWindow"]
#   set cvnames [getargs $args "-name" "No -name flag in call to landmarkWindow"] 
   set landmarkConfig([list $tag mode]) "minmax" 
   set landmarkConfig([list $tag isperidic]) 0
   set landmarkConfig([list $tag period]) 0

   if { ![file exists "$env(smapdir)/bin/dimlandmark"] } {
      tk_messageBox -icon error -type ok -title Message -message "Could not find dimlandmark code.  You must define an smapdir inside your .bashrc file.  dimlandmark should then be in \$smapdir/bin/dimlandmark"
      return
   }

   frame $wind -padx 1m -pady 1m
   pack [ createColvarSelector $wind.sel -datatag $dtag -selectiontag "print" -text "select the cvs to use" -ncols 6 ] -in $wind -side top -fill both
   # Needs an entry for the number of points and a menu for the mode we are using
   frame $wind.bfr -padx 1m -pady 1m
   label $wind.bfr.plab1 -text "periodic variables"
   checkbutton $wind.bfr.isp -variable cvlist::landmarkConfig([list $tag isperiodic])  
   label $wind.bfr.plab2 -text "period equals"
   entry $wind.bfr.pval -width 5 -textvariable cvlist::landmarkConfig([list $tag period])
   label $wind.bfr.npoints -text "Number of landmarks"
   entry $wind.bfr.np -width 5 -textvariable cvlist::landmarkConfig([list $tag npoints])
   label $wind.bfr.lopt -text "landmark selection mode"
   menubutton $wind.bfr.opt -relief raised -bd 2 -direction flush -width 5 \
         -textvariable cvlist::landmarkConfig([list $tag mode]) -menu $wind.bfr.opt.menu
   menu $wind.bfr.opt.menu
   $wind.bfr.opt configure -state disabled
   $wind.bfr.opt.menu delete 0 end
   $wind.bfr.opt configure -state normal
   $wind.bfr.opt.menu add radiobutton -label stride -value stride -variable cvlist::landmarkConfig([list $tag mode])
   $wind.bfr.opt configure -state normal
   $wind.bfr.opt.menu add radiobutton -label random -value rnd -variable cvlist::landmarkConfig([list $tag mode])
   $wind.bfr.opt configure -state normal
   $wind.bfr.opt.menu add radiobutton -label fps -value minmax -variable cvlist::landmarkConfig([list $tag mode])

   grid $wind.bfr.plab1 -row 1 -column 1 -sticky we 
   grid $wind.bfr.isp -row 1 -column 2 -sticky we
   grid $wind.bfr.plab2 -row 1 -column 3 -sticky we
   grid $wind.bfr.pval -row 1 -column 4 -sticky we
   grid $wind.bfr.npoints -row 2 -column 1 -sticky we
   grid $wind.bfr.np -row 2 -column 2 -sticky we
   grid $wind.bfr.lopt -row 2 -column 3 -sticky we
   grid $wind.bfr.opt -row 2 -column 4 -sticky we
   pack $wind.bfr

   return $wind
}

proc cvlist::createColvarSelector {wind args} {
  variable data
  variable colvarSelection

  set dtag [getargs $args "-datatag" "No -datatag flag in call to createColvarSelector"]
  set tag $dtag:[findMolecule]
#  set sname [getargs $args "-name" "No -name flag in call to createColvarSelector"]
  set ncols [getargs $args "-ncols" "No -ncols flag in call to createColvarSelector"]
  set seltag [getargs $args "-selectiontag" "No -selectiontag in call to createColvarSelector"]
  set textl [getargs $args "-text" "No -text in call to createColvarSelector"]

#  if { $sname=="all" } {
   set ncolvar 0
   foreach name $data([list $tag cv_names ]) {
      set colvarSelection([list $tag $name $seltag]) 0
      if { $name!="time" } { 
        set cv_names($ncolvar) $name
        incr ncolvar
      }
   }
#   } else {
#      set ncolvar 0
#      foreach name [lsearch -all -inline $data([ list $tag cv_names ]) $sname ] {
#        set cv_names($ncolvar) $name
#        incr ncolvar
#      }
#      foreach name $data([ list $tag cv_names ]) { set colvarSelection([list $tag $name $seltag]) 0 }
#   }
   # Calculate the number of rows
   set nrows [expr int( $ncolvar / $ncols ) + 1]

   labelframe $wind -relief ridge -bd 2 -text $textl -padx 2m -pady 2m
   # This creates all our little cv checkboxes - moment of gloating here this is so cool
   set n 0
   frame $wind.fcv -padx 3m -pady 3m
   for { set i 0 } { $i<$nrows } { incr i } {
      for { set j 0 } { $j<$ncols } { incr j } {
          set n [expr $i*$ncols + $j]
          if { $n>=$ncolvar } { break }
          label $wind.fcv.l$cv_names($n) -text "$cv_names($n)"
          checkbutton $wind.fcv.b$cv_names($n) -variable cvlist::colvarSelection([list $tag $cv_names($n) $seltag])
          grid $wind.fcv.l$cv_names($n) -row [expr $i+1] -column [expr 4*$j+1]
          grid $wind.fcv.b$cv_names($n) -row [expr $i+1] -column [expr 4*$j+2]
      }
      if { $n>=$ncolvar } { break }
   }
   pack $wind.fcv
   return $wind
}

proc cvlist::selectCollectiveVariables {args} {
   variable data
   variable colvarSelection

   set tag [getargs $args "-moltag" "no moltag specified in call to selectCollectiveVariables"]
   set variables [getargs $args "-variables" "no variables specified in call to selectCollectiveVariables"]
   set seltag [getargs $args "-selectiontag" "no -selectiontag specified in call to selectCollectiveVariables"] 

   # First unselct all the cvs
   foreach cvthere $data([list $tag cv_names]) { set colvarSelection([list $tag $cvthere $seltag]) 0 }

   # Now select the cvs that are actually required
   foreach cvreq $variables {
      set found 0
      foreach cvthere $data([list $tag cv_names]) {
        if { $cvreq==$cvthere } { set colvarSelection([list $tag $cvthere $seltag]) 1 ; set found 1 ; break }   
      }
      if { $found==0 } {
        tk_messageBox -icon error -type ok -title Message -message "You are trying to select a colvar that isn't present"
        return
      }
   }
}

proc cvlist::printSelectedColvars {args} {
   variable data
   variable colvarSelection 
  
   set dtag [getargs $args "-datatag" "No -datatag flag in call to printSelectedColvars"]
   set pos [lsearch $args "-molecule"]
   if { $pos<0 } {
      set printtag $dtag:[findMolecule]
   } else {
      set printtag $dtag:[ lindex $args [expr $pos+1] ]
   } 
   #set filen [getargs $args "-to" "No -to flag in call to printSelectedColvars"]
   set od [getargs $args "-to" "No -to flag in call to printSelectedColvars"]
   set seltag [getargs $args "-selectiontag" "No -selectiontag in call to printSelectedColvars"]

   # Based on the input setup how much of the data we are printing out
   set start -1  ; set end -1
   set pos [lsearch $args "-all"]
   set molno [lindex [split $printtag :] 1]
   if { $pos>0 } { 
      set start 0
      set end [molinfo $molno get numframes]
   } else {
      set start [getargs $args "-start" "No -start flag in call to printSelectedColvars use -start and -end or -all"]
      set end [getargs $args "-end" "No -end flag in call to printSelectedColvars use -start and -end or -all"]
   }

   # Check for errors in start and end points
   if { $start>=$end } { 
      puts "Error: start and end points for colvar writing are invalid"
      return 0
   }
   if { $start>[molinfo $molno get numframes] || $end>[molinfo $molno get numframes] } {
      puts "Error: start or end point for colvar writing is greater than total number of points"
      return 0
   }

   # Count the number of cvs and store the names of we have data for
   set sub_names {}
   foreach name $data([list $printtag cv_names]) {
       if { $colvarSelection([list $printtag $name $seltag])==1 } { lappend sub_names $name }
   }
   
   # Check that cvs have been selected
   if { [llength $sub_names]==0 } { 
      tk_messageBox -icon error -type ok -title Message -message "No cvs have been selected"
      return 0
   }

   # Output the data
   for { set i $start } { $i<$end } { incr i } {
       set tmplist {}
       foreach name $sub_names {
           lappend tmplist [lindex $data([list $printtag $name]) $i]  
       }
       puts $od "$tmplist"
   }
   return [llength $sub_names]
}

proc cvlist::clone {args} {
   variable data

   set tmpd [getargs $args "-tmpdir" "No tempory directory given in call to clone"]
   set ntag [getargs $args "-gtag" "Call to clone without -gtag flag"]
   set selection [getargs $args "-selection" "Call to clone without selection flag"]

   puts "In clone selection is $selection"
   set owd [pwd] ; cd $tmpd 
   set nout 0
   foreach fp $selection {
       [ atomselect top all frame $fp ] writepdb select_tmp_$nout.pdb
       incr nout
   }
   
   # Find the current top molecule so we can copy its represetations
   set molid [findMolecule]
   if { $nout>0 } {
     # tell vmd to read in the data and load it as a new molecule (N.B. the selection will be the new top molecule)
     mol new select_tmp_0.pdb type pdb
     for { set i 1 } { $i<$nout } { incr i } { mol addfile select_tmp_$i.pdb type pdb }

     mol delrep 0 top     ; # Delete the default representation of the new molecule (we copy those of the old)

     # Get the list of representations for the old molecule 
     set nreps [molinfo $molid get numreps]
     for { set i 0 } { $i<$nreps } { incr i } {
         set toget "{ rep $i} {selection $i} {color $i}"
         set rep [ molinfo $molid get $toget ]
         mol addrep top
         mol modstyle $i top [lindex $rep 0]
         mol modselect $i top [lindex $rep 1]
         mol modcolor $i top [lindex $rep 2]
     }
     # Clean up all the pdb files we have created
     eval file delete [glob -dir [pwd] select_tmp_*]
   }
   cd $owd

   # Find the current top molecule (this is the new one)
   set newmol [findMolecule]
   
   # And copy the new molecules cv data
   set oldtag $ntag:$molid ; set newtag $ntag:$newmol
   set data([list $newtag cv_names]) $data([list $oldtag cv_names])
   foreach frame $selection { 
      foreach cv_name $data([ list $oldtag cv_names ]) {
         lappend data([ list $newtag $cv_name ]) [lindex $data([list $oldtag $cv_name]) $frame] 
      } 
   }
   set data($newtag) 1
   mol top $molid
   return $newmol
}

proc cvlist::getDataRow {args} {
   variable data
   set dtag [getargs $args "-datatag" "Call to exists without -name flag"]
   set tag $dtag:[findMolecule]
   set row [getargs $args "-row" "Call to getData without -row flag"]

   foreach name $data([ list $tag cv_names ]) {
       if { $name!="time" } {
          lappend outd [lindex $data([list $tag $name]) $row]
       }
   } 
    
   return $outd
}

proc cvlist::getDataCol {args} {
   variable data 
   set dtag [getargs $args "-datatag" "Call to exists without -name flag"]
   set tag $dtag:[findMolecule]
   set name [getargs $args "-name" "Call to getData without -name flag"]

   if { ![ info exists data([list $tag $name]) ] } { return "error" }
   if { [ llength $data([list $tag $name]) ]!=[molinfo top get numframes] } { return "error" }   

   return $data([list $tag $name])
}

proc cvlist::exists {args} {
   variable data
   set tag [getargs $args "-datatag" "Call to exists without -datatag flag"]
   if { [info exists data($tag:[findMolecule])] } { 
       return 1 
   } else {
       return 0
   }
}

proc cvlist::updateMenu {mnu args} {
   variable data

   set dtag [getargs $args "-datatag" "Call to updateMenu without -datatag flag"]
   set tag $dtag:[findMolecule]
   set var [getargs $args "-variable" "Call to updateMenu with variable specification"]

   $mnu.menu delete 0 end
   $mnu configure -state disabled

   # Create none switches
   $mnu configure -state normal
   $mnu.menu add radiobutton -label "none" -value "none" -variable $var

   if { ![info exists data($tag)] } { set $var "none"; return }

   foreach cvname $data([list $tag cv_names]) {
      $mnu configure -state normal
      $mnu.menu add radiobutton -value $cvname -label "$cvname" -variable $var
   }
   set $var "none"
}

proc cvlist::addDataCol {args} {
   variable data

   set dtag [getargs $args "-datatag" "Call to cvlist::addDataCol without -datatag flag"]
   set pos [lsearch $args "-molecule"]
   if { $pos<0 } {
      set tag $dtag:[findMolecule]
   } else {
      set tag $dtag:[ lindex $args [expr $pos+1] ]
   }
   set name [getargs $args "-name" "Call to cvlist::addDataCol without -name flag"]
   set newdata [getargs $args "-data" "Call to cvlist::addDataCol without -data flag"]

   puts "In add data col new data is named $name and is to be added to datatag $tag"

   if { [info exists data([list $tag cv_names])] } {
     foreach oname $data([list $tag cv_names]) {
        if { $oname==$name } { 
           tk_messageBox -icon error -type ok -title Message -message "Name of new data column has already been used.  Reset molecule and reload"
           return
        }
     }
   } else {
     set data($tag) 1
   }

   # Make sure we have the number for the molecule we are loading data into
   set molno [lindex [split $tag :] 1] 
   if { [llength $newdata]!=[molinfo $molno get numframes] } {
        tk_messageBox -icon error -type ok -title Message -message "Size of new data column does not match number of frames [llength $newdata] [molinfo top get numframes]"
        return 
   }

   # Add the data to the array
   lappend data([list $tag cv_names]) $name
   for { set i 0 } { $i<[molinfo top get numframes] } { incr i } {
       lappend data([list $tag $name]) [lindex $newdata $i]
   }
}

proc cvlist::deleteTopMolData {args} {
  variable data

  set dtag [getargs $args "-datatag" "Call to deleteTopMolData without -datatag flag"]
  set tag $dtag:[findMolecule]

  foreach { l1 val } [array get data] {
     if { [lindex $l1 0]==$tag } { unset data($l1) }
  }
}

proc cvlist::readColvar {args} {
   set filename [getargs $args "-filename" "Call to read colvar but filename not specified"]
   set dtag [getargs $args "-datatag" "Call to readColvar without -datatag flag"]

   set fd [open $filename r]                 ; # Open COLVAR FILE
   set ncv [ expr 0 ]                        ; # Counter over Column number
   gets $fd header                           ; # Get the first line of the COLVAR file
   set sheader [regexp -inline -all -- {\S+} $header]
 
   set ncv [llength $sheader]                ; # Get the number of fields in the COLVAR file

   if { [lindex $sheader 1] != "FIELDS" } {
      tk_messageBox -icon error -type ok -title Message -message "Input file is not a COLVAR file - cannot read"
      return 0
   }

   # Get the names of the cvs
   set cvnames {}
   for { set id 2 } { $id < $ncv } { incr id } { lappend cvnames [lindex $sheader $id] }

   # Read the colvar file
   set ncv [expr $ncv - 2]    ; # Don't count the first two fields in the colvar header
   set ncolvar 0              ; # Counter of the number of colvars read in
   while { [gets $fd line] !=-1 } {
      set sline [regexp -inline -all -- {\S+} $line]
      for { set id 0 } { $id < $ncv } { incr id } { lappend tmpdata([lindex $cvnames $id]) [lindex $sline $id] }
      incr ncolvar
   }
   close $fd

   # Now add the data to the data array
   for { set id 0 } { $id < $ncv } { incr id } {
       set cvname [lindex $cvnames $id]
       addDataCol -datatag $dtag -name $cvname -data $tmpdata($cvname)
   }
   return $ncolvar
}

# This gets arguments or returns a warning   
proc cvlist::getargs {arglist tag {warning 0} } {
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

proc cvlist::findMolecule {args} {
   # Get the new molid
   foreach m [molinfo list] {
      if { [ molinfo $m get top ] } { set topmol $m }
   }
   return $topmol
}  

