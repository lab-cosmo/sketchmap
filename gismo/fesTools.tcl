# 
# Utilities for drawing free energy surfaces
#
# $Id: fesTools.tcl, v 1.0 2011/04/24
# 

source $env(VMDDIR)/scripts/init.d/gtplot.tcl

package provide fesTools 1.0
package require gtPlot 2.0

namespace eval ::fesTools {
# Procedures to export
    namespace export plotFes
    namespace export unpotFes
    namespace export appearance
    namespace export integrate
    namespace export clearSelection
    namespace export plotFesDiff
    namespace export drawConvergencePlot
    namespace export updateMenu
    namespace export createSpinbox
    namespace export configureSpinbox
    namespace export deleteTopMolData
# Variables to export
    namespace export kT
    variable fescol "jet"      ; # The color scale for free energy surfaces
# The variables for our free energy surfaces
    variable xcoord            ; # Stores the current x coordinate
    variable ycoord            ; # Stores the current y coordinate 
    variable feslist           ; # This is the list of all the free energy surfaces we have 
    variable fesdata           ; # We use this to store data on the free energy surface
    variable sumhills          ; # We use this to store the data on how to sum hills
    variable kT                ; # This stores the value of kT
}

proc fesTools::createSpinbox {sbox args} {
   variable fesdata
   lappend fesdata(frames) 1
   set var [getargs $args "-textvariable" "call to createSpinbox without -textvariable flag"]
   eval spinbox $sbox -textvariable $var -values $fesdata(frames) -width 5
}

proc fesTools::configureSpinbox {sbox args} {
   variable fesdata

   set dir [getargs $args "-directory" "Call to countFes without -directory flag"]
   set var [getargs $args "-textvariable" "call to createSpinbox without -textvariable flag"]
   set nam [ getFesFilename -directory $dir ]
 
   set fesdata(frames) [ countFes -directory $dir -name $nam ]
   $sbox configure -values $fesdata(frames)
   set $var [lindex $fesdata(frames) [expr [ llength $fesdata(frames) ] - 1] ]
}

proc fesTools::fesCreator {args} {
   variable feslist
   variable xcoord
   variable ycoord

   set name [getargs $args "-style" "Call to fesCreator without -style flag"]
   set filename [getargs $args "-filen" "Call to fesCreator without -filename flag"]
   set tmpdir [getargs $args "-tmpd" "Call to fesCreator without -tmpd flag"]

   if { ![info exists feslist(totalfes) ] } {
      set feslist(totalfes) 1
   } else {
      incr feslist(totalfes)
   }

   if { $name=="readfes" } {
     set ans [tk_messageBox -icon question -type yesno -title Message -message "This is a free energy surface for cv $xcoord vs cv $ycoord?"]
     if { $ans=="no" } { return "error" }
   } elseif { $name=="sumhills" } {
     set ans [tk_messageBox -icon question -type yesno -title Message -message "You will sum hills to get free energy as a function of cvs $xcoord and $ycoord?"]
     if { $ans=="no" } { return "error" }
   }

   set molid [findMolecule]
   if { ![info exists feslist([list $molid $xcoord $ycoord name])] } {
      set feslist([list $molid $xcoord $ycoord name]) {}
      set feslist([list $molid $xcoord $ycoord dname]) {}
      set feslist([list $molid $xcoord $ycoord fname]) {}
   }
   # Store the name of this free energy surface 
   lappend feslist([list $molid $xcoord $ycoord name]) $name
   # Create a directory to store this fes
   set tmpd "$tmpdir/Fes_$feslist(totalfes)"
   file mkdir $tmpd
   lappend feslist([list $molid $xcoord $ycoord dname]) $tmpd

   if { $name=="readfes" } {
      # Make sure we know what this fes data is called
      set namelist [ lassign [split $filename /] ]
      set pos [ expr [llength $namelist] - 1 ]
      set nam [lindex $namelist $pos]
      # Get the directory name containing the free energy surfaces
      for {set i 0} {$i<$pos} {incr i} { lappend tpath [lindex $namelist $i] }
      set path [join $tpath /]
      # Copy all files to the tempory directory    
      eval file copy -force [glob -dir $path $nam*] $tmpd
      lappend feslist([list $molid $xcoord $ycoord fname]) $nam
      return $tmpd
   } elseif { $name=="sumhills" } {
      # Store the name of the free energy surface
      lappend feslist([list $molid $xcoord $ycoord fname]) "fes.dat"
      # Open a window to run sum hills
      openSumHills -filename $filename -tmpdir $tmpd
   }
}

proc fesTools::openSumHills {args} {
    variable sumhills
    variable kT
    global env

    # Check sum hills is present     
    if { ![file exists "$env(plumedir)/utilities/sum_hills/sum_hills"] } { 
       tk_messageBox -icon error -type ok -title Message -message "Could not find sum_hills should be in $env(plumedir)/utilities/sum_hills/sum_hills"
       return
    }
  
    # Get everything from args       
    set filename [getargs $args "-filename" "Call to openSumHills without directive -filename"] 
    set sumhills(filename) $filename
    set sumhills(tmpd) [getargs $args "-tmpdir" "Call to openSumHills without directive -tmpdir"]

    set fd [open $filename r]                        ; # Open the HILLS file
    gets $fd line                                    ; # Read the first line
    close $fd                                        ; # Close the HILLS file
    set sline [regexp -inline -all -- {\S+} $line]   ; # Convert the line to a list
    set nline [llength $sline]                       ; # Count the number of elements
    set ncv [expr ( $nline - 3 ) / 2 ]               ; # Calculate the number of cvs  

    # Check for integral ncv and that there enough cvs in the input COLVAR
    if { [expr $ncv - int($ncv)] > 0 } {
       tk_messageBox -icon error -type ok -title Message -message "Found a non-integral number of cvs in HILLS file abandoning"
       return
    }
    set sumhills(ncv) $ncv          ; # Store number of cvs
    set sumhills(stride) 1000       ; # A tentative stride
    set kT(T) 300                   ; # A tentative temperature 
    set kT(units) 0.00831447        ; # Default value for kT is kJ mol-1

    set hills [toplevel ".hills"]
    wm title $hills "Sum Hills"
    wm resizable $hills 0 0

    frame $hills.top -padx 1m -pady 1m
    frame $hills.fmess -padx 5m -pady 5m
    pack [label $hills.mess -text "Please specify the all information required to sum hills" -anchor s] -in $hills.fmess -side top
    pack $hills.fmess -in $hills.top -side top -fill both

    # Stride and kT
    frame $hills.mid -padx 1m -pady 1m
    labelframe $hills.dat -relief ridge -bd 2 -text "General" -padx 2m -pady 2m
    label $hills.dat.lstride -text "stride" -anchor e
    entry $hills.dat.stride -width 5 -textvariable fesTools::sumhills(stride)
    kTMenu $hills.dat.kT
    grid $hills.dat.lstride -row 1 -column 1 -sticky e
    grid $hills.dat.stride  -row 1 -column 2 -sticky e
    grid $hills.dat.kT      -row 2 -column 1 -columnspan 2 -sticky e
    pack $hills.dat -in $hills.mid -side left -fill both

    # CV Stuff
    labelframe $hills.cvs -relief ridge -bd 2 -text "Select CVs" -padx 2m -pady 2m
 
    label $hills.cvs.title -text "cv number" -anchor e
    label $hills.cvs.activ -text "active" -anchor e
    label $hills.cvs.period -text "period" -anchor e
    grid $hills.cvs.title -row 1 -column 1 -sticky e
    grid $hills.cvs.activ -row 1 -column 2 -sticky e
    grid $hills.cvs.period -row 1 -column 3 -sticky e
    for { set i 1 } { $i<=$ncv } { incr i } {
       label $hills.cvs.t$i -text "cv$i" -anchor e
       set fesTools::sumhills([list active $i]) 0
       checkbutton $hills.cvs.a$i -variable fesTools::sumhills([list active $i])
       menubutton $hills.cvs.p$i -relief raised -bd 2 -direction flush -width 5 \
              -textvariable fesTools::sumhils([list period $i]) -menu $hills.cvs.p$i.menu
       menu $hills.cvs.p$i.menu
       set sumhills([list period $i]) "none"
       $hills.cvs.p$i.menu delete 0 end
       $hills.cvs.p$i configure -state disabled
       $hills.cvs.p$i configure -state normal
       $hills.cvs.p$i.menu add radiobutton -label "none" -value "none" -variable fesTools::sumhills([list period $i])
       $hills.cvs.p$i configure -state normal
       $hills.cvs.p$i.menu add radiobutton -label "pi" -value pi -variable fesTools::sumhills([list period $i])
       $hills.cvs.p$i configure -state normal
       $hills.cvs.p$i.menu add radiobutton -label "2pi" -value 2pi -variable fesTools::sumhills([list period $i])
       grid $hills.cvs.t$i -row [expr $i+1] -column 1 -sticky e
       grid $hills.cvs.a$i -row [expr $i+1] -column 2 -sticky e
       grid $hills.cvs.p$i -row [expr $i+1] -column 3 -sticky e
    }
    pack $hills.cvs -in $hills.mid -side right -fill both
    pack $hills.mid -in $hills.top -side top -fill both

    # Setup the buttons
    frame $hills.but -padx 3m -pady 3m
    pack [button $hills.but.cancel -text "Cancel" -relief raised -command {destroy .hills ; return} ] -in $hills.but -side right
    pack [button $hills.but.ok -text "Sum Hills" -relief raised -command {[namespace code ::fesTools::sumHills] ; destroy .hills } ] -in $hills.but -side right
    pack $hills.but -in $hills.top -side top -fill both
    pack $hills.top -fill both
}

proc fesTools::sumHills {args} {
    variable sumhills
    variable kT
    global env

    # Calculate the value of kT
    set kT(kT) [expr $kT(T)*$kT(units)]

    # Setup array containing all period data
    set period {} 
    for { set i 1 } { $i<=$sumhills(ncv) } { incr i } {
       if { $sumhills([list active $i])==1 } { lappend ndw $i }
       if { $sumhills([list period $i])=="2pi" } {
          lappend period "-2pi $i"
       } elseif {$sumhills([list period $i])=="pi" } {
          lappend period "-pi $i"
       }
    }

    # These are all our sanity checks
    if { $kT(kT)==0 } {
       tk_messageBox -icon error -type ok -title Message -message "kT was not set correctly - check units"
       return
    } elseif { ![string is integer -strict $sumhills(stride)] } {
       tk_messageBox -icon error -type ok -title Message -message "stride is not an integer - please remedy"
       return
    } elseif { [llength $ndw] != 2 && [llength $ndw] != 1 } {
       tk_messageBox -icon error -type ok -title Message -message "Please select one or two cvs only"
       return
    }

    # Copy HILLS file
    file copy -force $sumhills(filename) $sumhills(tmpd)
    set owd [pwd]   ;    cd $sumhills(tmpd)

    # Get the name of the hills file (i.e. without the path) 
    set namelist [ lassign [split $sumhills(filename) /] ]
    set pos [ expr [llength $namelist] - 1 ]
    set fname [lindex $namelist $pos]

    puts "Executing: $env(plumedir)/utilities/sum_hills/sum_hills -ndim $sumhills(ncv) -ndw [join $ndw] [join $period] -stride $sumhills(stride) -kt $kT(kT) -file $fname -out fes.dat"
    catch {eval exec $env(plumedir)/utilities/sum_hills/sum_hills -ndim $sumhills(ncv) -ndw [join $ndw] [join $period] -stride $sumhills(stride) -kt $kT(kT) -file $fname -out fes.dat} sumhills_stdout

    # Check it worked
    if { ![file exists $sumhills(tmpd)/fes.dat] } {
       puts $sumhills_stdout
       puts "-------------"
       tk_messageBox -icon error -type ok -title Message -message "Something went wrong when running sum_hills.  Please see messages in log"
       return
    }
    cd $owd
}

proc fesTools::getFesFilename {args} {
   variable feslist 
   variable xcoord
   variable ycoord

   set dir [getargs $args "-directory" "Call to getFesFilename without -directory flag"]

   # Find the name of the molecule   
   set molid [findMolecule]
   for { set i 0 } { $i<[llength $feslist([list $molid $xcoord $ycoord dname])] } { incr i } {
       if { [lindex $feslist([list $molid $xcoord $ycoord dname]) $i]==$dir } {
            set ncoord [lindex $feslist([list $molid $xcoord $ycoord fname]) $i]
       }
   }
   return $ncoord 
}

proc fesTools::drawConvergencePlot {args} {
    variable kT

    set dir [getargs $args "-directory" "Call to countFes without -directory flag"]
    #set nam [getargs $args "-name" "Call to countFes without -name flag"]
    set nam [ getFesFilename -directory $dir ] 

    # Loop over all the free energy surfaces
    set n 1
    set fesl [countFes -directory $dir -name $nam ]
    foreach no $fesl {
       if { $no==[ lindex $fesl [expr [llength $fesl]-1] ] } {
          readFes "$dir/$nam" festmp
       } else {
          readFes "$dir/$nam.$no" festmp
       }   
       lappend fest $n
       lappend fesval [calcFesDiff festmp]
       incr n
    }

    # And plot the data
    set convplot [multiplot -title "Free energy as a function of time" -xlabel "Frame" -nostats]
    $convplot add $fest $fesval -legend "Free energy" -lines -marker circle -radius 2 -fillcolor blue -color blue -nostats
    $convplot replot
}

proc fesTools::plotFesDiff {canv args} {
    variable fesdata
    variable fescol
    variable kT
      
    set dir [getargs $args "-directory" "Call to countFes without -directory flag"]
    set nam [ getFesFilename -directory $dir ]
    set no  [getargs $args "-number" "Call to plotFes without -number flag"]
    set datatag [getargs $args "-datatag" "Call to plotFes without -datatag flag"]
          
    # Get rid of any old data
    ::gtPlot::unplotData $canv -datatag $datatag

    set fesl [countFes -directory $dir -name $nam]
    if { $no==[ lindex $fesl [expr [llength $fesl]-1] ] } {
       tk_message -icon error -type ok -title Message -message "Plotting difference between current fes and final fes is pointless given current fes is final one"
    } else {
       set filename "$dir/$nam.$no"
    }

    # Read the fess
    readFes $filename fes1
    readFes "$dir/$nam" fes2
    # Compute the difference between free energies
    if { [fesDiff fes1 fes2 fesd]=="error" } { return }
    if { $fesdata([list fesd dim])==1 } {
        ::gtPlot::readData -datatag $datatag -npoints $fesdata([list fesd np]) -xdata $fesdata([list fesd x]) -ydata $fesdata([list fesd z]) -colors "black"
        ::gtPlot::plotData $canv -datatag $datatag -plotstyle "line" -linecolor "black"
    } elseif { $fesdata([list fesd dim])==2 } {
        ::gtPlot::readData -datatag $datatag -npoints $fesdata([list fesd np]) -xdata $fesdata([list fesd x]) -ydata $fesdata([list fesd y]) -zdata $fesdata([list fesd z])
        ::gtPlot::plotData $canv -datatag $datatag -plotstyle "surface" -surfacecolors $fescol ; #-ncolors $nfescolors
    } else {
        puts "Yes this is unlikely to ever be implemented by me - if you are a genius go for it!"
        return
    }

    if { ![info exists kT(kT)] } { kTenquiry }
    set fesdata(plotexists) 1
}

proc fesTools::unplotFes {canv args} {
    variable fesdata

    set datatag [getargs $args "-datatag" "Call to plotFes without -datatag flag"]
           
    # Get rid of the plotted data
    ::gtPlot::unplotData $canv -datatag $datatag
    # Set it so that free energy binds don't do anything anything
    if { [info exists fesdata(plotexists)] } { unset fesdata(plotexists) }
}

proc fesTools::plotFes {canv args} {
    variable fesdata
    variable fescol
    variable kT

    set dir [getargs $args "-directory" "Call to countFes without -directory flag"]
    set nam [ getFesFilename -directory $dir ]
    set no  [getargs $args "-number" "Call to plotFes without -number flag"]
    set datatag [getargs $args "-datatag" "Call to plotFes without -datatag flag"]

    # Get rid of any old data
    ::gtPlot::unplotData $canv -datatag $datatag

    # Get the filname to read in 
    set fesl [countFes -directory $dir -name $nam]    
    if { $no==[ lindex $fesl [expr [llength $fesl] - 1] ] } {
        set filename "$dir/$nam"
    } else { 
        set filename "$dir/$nam.$no"
    }

    # Read the file
    readFes $filename fes
    if { $fesdata([list fes dim])==1 } {
        ::gtPlot::readData -datatag $datatag -npoints $fesdata([list fes np]) -xdata $fesdata([list fes x]) -ydata $fesdata([list fes z]) -colors "black"
        ::gtPlot::plotData $canv -datatag $datatag -plotstyle "line" -linecolor "black"
    } elseif { $fesdata([list fes dim])==2 } {
        ::gtPlot::readData -datatag $datatag -npoints $fesdata([list fes np]) -xdata $fesdata([list fes x]) -ydata $fesdata([list fes y]) -zdata $fesdata([list fes z])
        ::gtPlot::plotData $canv -datatag $datatag -plotstyle "surface" -surfacecolors $fescol ; # -ncolors $nfescolors 
    } else {
        puts "Yes this is unlikely to ever be implemented by me - if you are a genius go for it!"
        return
    }

    if { ![info exists kT(kT)] } { kTenquiry }
    set fesdata(plotexists) 1
}

proc fesTools::clearSelection {canv args} {
    variable fesint
    variable fesdata    

    if { ![info exists fesdata(plotexists)] } { return }

    if { [info exists fesint(num)] } {
       unset fesint(num) 
       $canv delete fesstuff 
    }
}

proc fesTools::integrate {canv args} {
   variable fesint
   variable fesdata
   variable kT

   if { ![info exists fesdata(plotexists)] } { return }
   if { ![info exists fesint(num)] } { set fesint(num) 0 }

   set datatag [getargs $args "-datatag" "Call to integrate without -datatag flag"]
   # Get the z data for this fes
   set zdata [::gtPlot::getLasso -datatag $datatag -data "zdata"]

   # Now do the integration itself
   set total 0
   foreach z $zdata { set total [ expr $total + exp(-$z/$kT(kT)) ] }
   set fes [ expr -$kT(kT)*log($total) ]

   # We must write the free energy on the canvas here 
   if { $fesint(num)==0 } {
      set style "solid"
      set fesint($fesint(num)) $fes
      set fesint([list $fesint(num) box]) [::gtPlot::getLasso -datatag $datatag -data "box" ]
      $canv create text 10 [ expr [ winfo height $canv ] - 20 ] -anchor sw -text "G(solid box) = [format %.4g $fesint(0)] \n\n" -tag "fesstuff festext"
   } elseif { $fesint(num)==1} {
      set style "dashed"
      set fesint($fesint(num)) $fes
      set fesint([list $fesint(num) box]) [::gtPlot::getLasso -datatag $datatag -data "box"]
      $canv delete festext
      $canv create text 10 [ expr [ winfo height $canv ] - 20 ] -anchor sw -text " G(solid box) = [format %.4g $fesint(0)] \n G(dashed box) = [format %.4g $fesint(1)] \n difference = [format %.4g [expr $fesint(0)-$fesint(1)]]" -tag "fesstuff" -justify left
   } else {
      clearSelection $canv
      set fesint(num) 0
      set style "solid"
      set fesint($fesint(num)) $fes
      set fesint([list $fesint(num) box]) [::gtPlot::getLasso -datatag $datatag -data "box"]
      $canv create text 10 [ expr [ winfo height $canv ] - 20 ] -anchor sw -text "G(solid box) = [format %.4g $fesint(0)] \n\n" -tag "fesstuff festext"
   }
   ::gtPlot::drawShape $canv -shape "rectangle" -corners $fesint([list $fesint(num) box]) -linestyle $style -shapetag fesstuff  
   incr fesint(num)    ; # This keeps track of what we are doing next time this routine is called
}

proc fesTools::countFes {args} {
    set dir [getargs $args "-directory" "Call to countFes without -directory flag"]
    set nam [ getFesFilename -directory $dir ]

    set nfes [expr [llength [glob -dir $dir $nam*] ] - 1 ]
    for { set i 1 } { $i<=[expr $nfes + 1] } { incr i } { lappend fesframes $i }  
    return $fesframes
}

proc fesTools::appearance {canv wind args} {
    frame $wind -padx 1m -pady 1m

    labelframe $wind.afram -relief ridge -bd 2 -text "FES appearance" -padx 2m -pady 2m
    label $wind.afram.clab -text "colorscale" -anchor e
    ::gtPlot::colorMenu $wind.afram.c -textvariable fesTools::fescol

    grid $wind.afram.clab -row 1 -column 1 -sticky e
    grid $wind.afram.c -row 1 -column 2 -sticky e
    pack $wind.afram -in $wind -side top -fill both

    pack [ ::gtPlot::zscale_window $canv $wind.zscale ] -in $wind -side top -fill both
    return $wind
}

proc fesTools::calcFesDiff {tag} {
    variable fesdata
    variable fesint
    variable kT
  
    if { ![info exists fesdata(plotexists)] } { return }
    if { ![info exists fesint(num)] } {
       tk_messageBox -icon error -type ok -title Message -message "You must select a free energy difference to compute to calculate convergence"
       return 
    }

    set np $fesdata([list $tag np])
    set dim $fesdata([list $tag dim]) 

    if { $fesint(num)==1 } {
       set xmin [lindex $fesint([list 0 box]) 0]
       set xmax [lindex $fesint([list 0 box]) 2]
       if { $dim==2 } {
         set ymin [lindex $fesint([list 0 box]) 1]
         set ymax [lindex $fesint([list 0 box]) 3]
       }
    } elseif { $fesint(num)==2 } {
       set xmin [lindex $fesint([list 0 box]) 0]
       set xmax [lindex $fesint([list 0 box]) 2]
       set xmin2 [lindex $fesint([list 1 box]) 0]
       set xmax2 [lindex $fesint([list 1 box]) 2]
       if { $dim==2 } {
         set ymin [lindex $fesint([list 0 box]) 1]
         set ymax [lindex $fesint([list 0 box]) 3]
         set ymin2 [lindex $fesint([list 1 box]) 1]
         set ymax2 [lindex $fesint([list 1 box]) 3]
       }     
    } else {
       tk_messageBox -icon error -type ok -title Message -message "If this error comes up then I don't understand my own programming"
       return
    }   

    set total1 0 ; set total2 0
    for {set i 0} { $i<$np } { incr i } {
        if { $dim==1 } {
           set x [lindex $fesdata([list $tag x]) $i]
           if { $x>$xmin && $x<$xmax } { set total1 [ expr $total1 + exp(-[lindex $fesdata([list $tag z]) $i]/$kT(kT)) ] } 
           if { $fesint(num)!=1 } {
              if { $x>$xmin2 && $x<$xmax2 } { set total2 [ expr $total2 + exp(-[lindex $fesdata([list $tag z]) $i]/$kT(kT)) ] }
           }
        } elseif { $dim==2 } {
           set x [lindex $fesdata([list $tag x]) $i]
           set y [lindex $fesdata([list $tag y]) $i]
           if { $x>$xmin && $x<$xmax && $y>$ymin && $y<$ymax } { set total1 [ expr $total1 + exp(-[lindex $fesdata([list $tag z]) $i]/$kT(kT)) ]  }
           if { $fesint(num)!=1 } {
              if { $x>$xmin2 && $x<$xmax2 && $y>$ymin2 && $y<$ymax2 } { set total2 [ expr $total2 + exp(-[lindex $fesdata([list $tag z]) $i]/$kT(kT)) ]  }
           }
        } 
    }
    if { $fesint(num)==1 } { 
       return [ expr -$kT(kT)*log($total1) ]
    } else {
       return [ expr -$kT(kT)*log($total1) + $kT(kT)*log($total2) ]
    }
}

proc fesTools::fesDiff {tag1 tag2 ntag} {
    variable fesdata

    # Check for stupid errors
    if { $fesdata([list $tag1 dim])!=$fesdata([list $tag2 dim]) } {
       tk_messageBox -icon error -type ok -title Message -message "Dimension mismatch for fes difference"
       return "error"
    } elseif { $fesdata([list $tag1 np])!=$fesdata([list $tag2 np]) } {
       tk_messageBox -icon error -type ok -title Message -message "Size mistmatch for fes difference"
       return "error"
    }

    set dim $fesdata([list $tag1 dim])
    set np $fesdata([list $tag1 np]) 
    set fesdata([list $ntag x]) {}   ; set fesdata([list $ntag y]) {} ; set fesdata([list $ntag z]) {}
    for {set i 0} { $i<$np } { incr i } {
        # Copy x data
        set x1 [lindex $fesdata([list $tag1 x]) $i]
        set x2 [lindex $fesdata([list $tag2 x]) $i]
        if { $x1!=$x1 } {
           tk_messageBox -icon error -type ok -title Message -message "Grid mismatch for fes difference"
           return "error"
        } 
        lappend fesdata([list $ntag x]) $x1

        # Get z data
        set z1 [lindex $fesdata([list $tag1 z]) $i]
        set z2 [lindex $fesdata([list $tag2 z]) $i]
        lappend fesdata([list $ntag z]) [expr $z1 - $z2]

        # Copy y data
        if { $dim==2 } {
          set y1 [lindex $fesdata([list $tag1 y]) $i]
          set y2 [lindex $fesdata([list $tag2 y]) $i]
          if { $x1!=$x1 } { 
             tk_messageBox -icon error -type ok -title Message -message "Grid mismatch for fes difference"
             return "error"
          } 
          lappend fesdata([list $ntag y]) $y1
        }  
    }
    set fesdata([list $ntag dim]) $dim
    set fesdata([list $ntag np]) $np 
}

proc fesTools::readFes {filename tag} {
    variable fesdata

    set fd [open $filename r]
    set np 0
    # Get first line of file (is this a one or 2d fes)
    gets $fd line ; set sline [regexp -inline -all -- {\S+} $line]
    
    # Establish whether this is a one or 2d fes
    if { [llength $sline]==3 } {
       set fesdim 2 ; incr np
       lappend x [lindex $sline 0] ; lappend y [lindex $sline 1] ; lappend z [lindex $sline 2]
    } elseif { [llength $sline]==2 } {
       set fesdim 1 ; incr np
       # Note we set y here to an empty list just so that later thing will work
       lappend x [lindex $sline 0] ; lappend z [lindex $sline 1] ; set y {}
    } elseif { [llength $sline]!=0 } {
        tk_messageBox -icon error -type ok -title Message -message "Found a weird line in fes"
        return
    }

    while {[gets $fd line] != -1} {
       set sline [regexp -inline -all -- {\S+} $line]
       if { [llength $sline]==3 && $fesdim==2 } {
          lappend x [lindex $sline 0]
          lappend y [lindex $sline 1]
          lappend z [lindex $sline 2]
          incr np
       } elseif { [llength $sline]==2 && $fesdim==1 } {
          lappend x [lindex $sline 0] ; lappend z [lindex $sline 1]
          incr np
       } elseif { [llength $sline]!=0 } {
           tk_messageBox -icon error -type ok -title Message -message "Found a weird line in fes"
           return
       }
    }
    close $fd

    # Save everthing - don't read these in so that we can overwrite old data
    set fesdata([list $tag np ]) $np  ; set fesdata([list $tag dim]) $fesdim
    set fesdata([list $tag x]) $x
    set fesdata([list $tag y]) $y
    set fesdata([list $tag z]) $z
}

proc fesTools::deleteTopMolData {args} {
  variable feslist

  if { ![info exists feslist] } { return }
 
  # Hopefully this should delete everything that has a 
  # first element that is equal to molid 
  set molid [findMolecule]
  foreach {l1 val} [array get feslist] {
    if { [lindex $l1 0]==$molid } { unset feslist($l1) }
  }
}

proc fesTools::updateMenu {mnu args} {
   variable feslist
   variable xcoord
   variable ycoord

   set xcoord [getargs $args "-xcoord" "call to updateMenu without specification of xcoord"]
   set ycoord [getargs $args "-ycoord" "call to updateMenu without specification of ycoord"]
   set var [getargs $args "-variable" "call to updateMenu without specification of variable"]

   # Update the menu which tells us what free energy data we have for this system
   $mnu.menu delete 0 end
   $mnu configure -state disabled
   # Create none swiches
   $mnu configure -state normal
   $mnu.menu add radiobutton -label "none" -value "none" -variable $var
   set molid [findMolecule]
   if { [info exists feslist([list $molid $xcoord $ycoord name])] } {
       set n 0
       foreach fesname $feslist([list $molid $xcoord $ycoord name]) {
         set fesdir [lindex $feslist([list $molid $xcoord $ycoord dname]) $n]
         $mnu configure -state normal
         $mnu.menu add radiobutton -value $fesdir -label "$fesname" -variable $var
         incr n
       }
   }
}

proc fesTools::kTenquiry {args} {
   variable kT
   # Set defaults for kT (300K and kJ mol-1) 
   set kT(T)                300     
   set kT(units)            0.00831447

   # Get rid of old kt enquiry window
   if [winfo exists .kt] { destroy .kt }

   set tw [toplevel ".kt"]
   wm title $tw "Set value of kT"
   wm resizable $tw 0 0

   frame $tw.top -padx 2m -pady 2m
   pack [ kTMenu $tw.top.stuff ] -side top -fill both
   pack [ button $tw.ok -text "OK" -relief raised -command { set fesTools::kT(kT) [expr $fesTools::kT(T)*$fesTools::kT(units)] ; destroy .kt } ] -in $tw.top -side bottom
   pack $tw.top -fill both
}

proc fesTools::kTMenu {w args} {
   variable kT
   
   frame $w
   label $w.tlab -text "Temperature" -anchor e
   entry $w.t -width 5 -textvariable fesTools::kT(T)
   label $w.ulab -text "Units" -anchor e   
   menubutton $w.u -relief raised -bd 2 -direction flush -width 5 \
            -textvariable fesTools::kT(units) -menu $w.u.menu
   menu $w.u.menu
   $w.u.menu delete 0 end 
   $w.u configure -state disabled

   $w.u configure -state normal
   $w.u.menu add radiobutton -label "kJ/mol" -value 0.00831447 -variable fesTools::kT(units)
   $w.u configure -state normal
   $w.u.menu add radiobutton -label "kcal/mol" -value 0.001987191 -variable fesTools::kT(units)
   $w.u configure -state normal
   $w.u.menu add radiobutton -label "eV" -value 0.00008617343 -variable fesTools::kT(units)
   $w.u configure -state normal
   $w.u.menu add radiobutton -label "dlpoly" -value 0.831451115 -variable fesTools::kT(units)
   $w.u configure -state normal
   $w.u.menu add radiobutton -label "Rydbergs" -value 0.0000063363125 -variable fesTools::kT(units)
   $w.u configure -state normal
   $w.u.menu add radiobutton -label "Natural" -value 1.0 -variable fesTools::kT(units)

   grid $w.tlab -row 1 -column 1 -sticky e
   grid $w.t -row 1 -column 2 -sticky e
   grid $w.ulab -row 2 -column 1 -sticky w
   grid $w.u -row 2 -column 2 -sticky e
   return $w
}

# This gets arguments or returns a warning   
proc fesTools::getargs {arglist tag {warning 0} } {
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

proc fesTools::findMolecule {args} {
   # Get the new molid
   foreach m [molinfo list] {
      if { [ molinfo $m get top ] } { set topmol $m }
   }
   return $topmol                
}  
 



