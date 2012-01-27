#
# Driver Interface v 1.0
#
# Gareth Tribello
#
# To do 
# (1) Needs better getargs in here

package provide driverInterface 1.0

namespace eval ::driverInterface {
# Procedures to export
  namespace export openLenUnitDialogue
  namespace export openTorsionsWindow
  namespace export openDRMSDWindow
  namespace export openRDFWindow
  namespace export openADFWindow
  namespace export readPlumed
# Variables to export
  namespace export status
# Other variables
  variable len            ;# Handle to the length unit window
  variable tor            ;# Handle to the torsion angle window
  variable drmsd          ;# Handle to the drmsd window
  variable rdf            ;# Handle to the rdf window
  variable adf            ;# Handle to the adf window
  variable tmpd           ;# The name of the tempory directory
  variable filename       ;# The naeme of the plumed input file
  variable tfilename      ;# The name of the topology file
  variable lenunit        ;# Units of length in plumed input
  variable backbone       ;# Do backbone torsions
  variable sidechain      ;# Do sidechain torsions 
  variable aselect1       ;# Storage space for an atom selection
  variable aselect2       ;# Storage space for an atom selection       
  variable aselect3       ;# Storage space for an atom selection
  variable lbound         ;# Lower bound for rdf calculation
  variable ubound         ;# Upper bound for rdf calculation
  variable window         ;# Size of windows in rdf/adf calculation
  variable gwidth         ;# Thicknesses of gaussians for rdf/adf calculation 
  variable switchf        ;# Parameters for ADF switching function
  variable status         ;# Allows us to report when stuff is read in to other functions
  variable broadcastvar   ;# Tells us who we are broadcasting to
}

proc driverInterface::openLenUnitDialogue {args} {
   variable len
   variable tmpd
   variable filename
   variable tfilename
   variable lenunit
   variable broadcastvar

   set tmpd [getargs $args "-tmpdir" {}]
   set filename [getargs $args "-filename" {}]
   set broadcastvar [getargs $args "-broadcast" {}]

   set tfilename "none"

   if [winfo exists .len] {
     wm deiconify $len
     return 
   }

   # Open len dialogue
   set len [toplevel ".len"]
   wm title $len "length units"
   wm resizable $len 0 0

   frame $len.top -padx 1m -pady 1m
   pack [label $len.lab -text "Please specify the length units"] -in $len.top -side top -fill both
   frame $len.up -padx 1m -pady 1m
   label $len.up.lab -text "units"
   menubutton $len.up.m -relief raised -bd 2 -direction flush -width 5 \
              -textvariable driverInterface::lenunit -menu $len.up.m.menu
   menu $len.up.m.menu
   grid $len.up.lab -row 1 -column 1 -sticky e
   grid $len.up.m   -row 1 -column 2 -sticky e
   pack $len.up -in $len.top -side top -fill both   

   pack [label $len.flab -text "topology file containing masses and charges"] -in $len.top -side top -fill both
   frame $len.file -padx 1m -pady 1m
   entry $len.file.fent -width 5 -textvariable driverInterface::tfilename
   button $len.file.browse -text "Browse" -relief raised -command {set driverInterface::tfilename [tk_getOpenFile -initialdir [pwd] -title "Select the topology file"] }
   grid $len.file.fent -row 1 -column 1 -sticky e
   grid $len.file.browse -row 1 -column 2 -sticky e
   pack $len.file -in $len.top -side top -fill both

   # Setup the menu
   $len.up.m.menu delete 0 end
   $len.up.m configure -state disabled
   $len.up.m configure -state normal
   $len.up.m.menu add radiobutton -label "Angstroms" -value 1.0 -variable driverInterface::lenunit
   $len.up.m configure -state normal
   $len.up.m.menu add radiobutton -label "nm" -value 10 -variable driverInterface::lenunit
   $len.up.m configure -state normal
   $len.up.m.menu add radiobutton -label "au" -value 0.5292 -variable driverInterface::lenunit

   # Setup the buttons
   frame $len.but -padx 3m -pady 3m
   pack [button $len.but.cancel -text "Cancel" -relief raised -command {destroy .len ; return} ] -in $len.but -side right
   pack [button $len.but.ok -text "OK" -relief raised -command {[namespace code driverInterface::readPlumed] ; destroy .len } ] -in $len.but -side right
   pack $len.but -in $len.top -side top -fill both
   pack $len.top -fill both

   # Assume that lenghts are in Angstroms
   set lenunit 1.0
}

proc driverInterface::openTorsionsWindow {args} {
  variable tor
  variable aselect1
  variable backbone
  variable sidechain
  variable tmpd 
  variable broadcastvar 

  set tmpd [getargs $args "-tmpdir" {}]
  set broadcastvar [getargs $args "-broadcast" {}]
  # Set defaults
  set aselect1 "protein"  ; set backbone 1 ; set sidechain 0

  if [winfo exists .tor] {
    wm deiconify $tor
    return 
  }

  # We now create a window asking the user what he would like
  set tor [toplevel ".tor"]
  wm title $tor "Select torsions"
  wm resizable $tor 0 0

  frame $tor.top -padx 1m -pady 1m
  pack [label $tor.lab -text "Please select which torsions you wish to use" -anchor s] -in $tor.top -side top
  frame $tor.sel -padx 1m -pady 1m
  pack [label $tor.sel.lab -text "In selection" -anchor s] -in $tor.sel -side left
  pack [entry $tor.sel.ent -width 20 -textvariable driverInterface::aselect1 ] -in $tor.sel -side right
  pack $tor.sel -in $tor.top -side top -fill both

  #torsion Options
  frame $tor.ang -padx 1m -pady 1m
  label $tor.ang.lback -text "backbone" -anchor s
  checkbutton $tor.ang.back -variable driverInterface::backbone   
  label $tor.ang.lside -text "sidechain" -anchor s  
  checkbutton $tor.ang.side -variable driverInterface::sidechain
  grid $tor.ang.lback -row 1 -column 1 
  grid $tor.ang.back  -row 1 -column 2
  grid $tor.ang.lside -row 1 -column 3 
  grid $tor.ang.side  -row 1 -column 4
  pack $tor.ang -in $tor.top -side top -fill both

  # The buttons
  frame $tor.but -padx 1m -pady 1m
  pack [button $tor.but.canc -text "Cancel" -relief raised -command {destroy .tor } ] -in $tor.but -side left
  pack [button $tor.but.ok -text "OK" -relief raised -command {destroy .tor ; [namespace code driverInterface::calcTorsions] } ] -in $tor.but -side right
  pack $tor.but -in $tor.top -side top -fill both
  pack $tor.top -fill both
}

proc driverInterface::calcTorsions {args} {
  variable backbone
  variable sidechain
  variable aselect1
  variable tmpd

  # Check user isn't a retard
  if { $backbone==0 && $sidechain==0 } {
     tk_messageBox -icon error -type ok -title Message -message "No torsional angles selected" 
     return
  }

  # Open our plumed output file
  set od [open "$tmpd/plumed.dat" w]
  puts "The tempory directory is $tmpd"
  puts $od "PRINT W_STRIDE 1"

  # Get the list of residues that we are interested in
  set sel [atomselect top "$aselect1 and alpha"]
  set reslist [ $sel get resid ]

  if { $backbone==1 } {
     for { set i 0 } { $i<[llength $reslist] } { incr i } {
         # This is the first torsion
         set sel [atomselect top "resid [expr [lindex $reslist $i] - 1] and type C"]
         set at1 [ expr [$sel get index] + 1 ]
         set sel [atomselect top "resid [lindex $reslist $i] and type N"]
         set at2 [ expr [$sel get index] + 1 ]
         set sel [atomselect top "resid [lindex $reslist $i] and alpha"]
         set at3 [ expr [$sel get index] + 1 ]
         set sel [atomselect top "resid [lindex $reslist $i] and type C"]
         set at4 [ expr [$sel get index] + 1 ]
         puts $od "TORSION LIST $at1 $at2 $at3 $at4"      
         # This is the second torsion
         set sel [atomselect top "resid [lindex $reslist $i] and type N"]
         set at1 [ expr [$sel get index] + 1 ]
         set sel [atomselect top "resid [lindex $reslist $i] and alpha"]
         set at2 [ expr [$sel get index] + 1 ]
         set sel [atomselect top "resid [lindex $reslist $i] and type C"]
         set at3 [ expr [$sel get index] + 1 ]
         set sel [atomselect top "resid [expr [lindex $reslist $i] + 1] and type N"]
         set at4 [ expr [$sel get index] + 1 ]
         puts $od "TORSION LIST $at1 $at2 $at3 $at4"
     }
  } 
  if { $sidechain==1 } {
     puts "I will implement this if I ever need it"
#     for { set i 0 } { $i<[llength $reslist] } { incr i } {
#     }
  }
  close $od

  runDriver   ; # Actually run driver
  return
}

proc driverInterface::openDRMSDWindow {args} {
  variable drmsd
  variable aselect1 
  variable aselect2
  variable tmpd
  variable broadcastvar
  variable switchf

  set tmpd [getargs $args "-tmpdir" {}]
  set broadcastvar [getargs $args "-broadcast" {}]

  # Set default values
  set aselect1 "protein" ; set aselect2 "type CA" ; 
  set switchf(d0) 0;   set switchf(r0) 0
  set switchf(nn) 6;   set switchf(mm) 12

  if [winfo exists .drmsd] {
    wm deiconify $drmsd
    return 
  }

  set drmsd [toplevel ".drmsd"]
  wm title $drmsd "Contact map controls"
  wm resizable $drmsd 0 0
   
  frame $drmsd.top -padx 1m -pady 1m
  pack [label $drmsd.lab -text "Select the atoms to involve the in the contact map"] -in $drmsd.top -side top -fill both  

  frame $drmsd.at -padx 1m -pady 1m
  label $drmsd.at.lab1 -text "Take atoms from"
  entry $drmsd.at.ent1 -width 20 -textvariable driverInterface::aselect1
  label $drmsd.at.lab2 -text "Atom types" 
  entry $drmsd.at.ent2 -width 20 -textvariable driverInterface::aselect2
  label $drmsd.at.sflab -text "Switching function parameters:"
  label $drmsd.at.r0l -text "r0"
  entry $drmsd.at.r0 -width 5 -textvariable driverInterface::switchf(r0)
  label $drmsd.at.d0l -text "d0"
  entry $drmsd.at.d0 -width 5 -textvariable driverInterface::switchf(d0)
  label $drmsd.at.nnl -text "nn"
  entry $drmsd.at.nn -width 5 -textvariable driverInterface::switchf(nn)
  label $drmsd.at.mml -text "mm"
  entry $drmsd.at.mm -width 5 -textvariable driverInterface::switchf(mm) 
  grid $drmsd.at.lab1 -row 1 -column 1 -columnspan 2
  grid $drmsd.at.ent1 -row 1 -column 3 -columnspan 6
  grid $drmsd.at.lab2 -row 2 -column 1 -columnspan 2
  grid $drmsd.at.ent2 -row 2 -column 3 -columnspan 6
  grid $drmsd.at.sflab -row 3 -column 1 -columnspan 8
  grid $drmsd.at.r0l -row 4 -column 1
  grid $drmsd.at.r0  -row 4 -column 2
  grid $drmsd.at.d0l -row 4 -column 3
  grid $drmsd.at.d0  -row 4 -column 4
  grid $drmsd.at.nnl -row 4 -column 5
  grid $drmsd.at.nn  -row 4 -column 6
  grid $drmsd.at.mml -row 4 -column 7
  grid $drmsd.at.mm  -row 4 -column 8
  pack $drmsd.at -in $drmsd.top -side top -fill both
   
  # And the buttons
  frame $drmsd.buttons -padx 1m -pady 1m
  pack [button $drmsd.buttons.cancel -text "cancel" -relief raised -command { destroy .drmsd ; return } ] -in $drmsd.buttons -side left
  pack [button $drmsd.buttons.ok -text "ok" -relief raised -command { destroy .drmsd; [namespace code driverInterface::doDRMSD] } ] -in $drmsd.buttons -side right
  pack $drmsd.buttons -in $drmsd.top -side top -fill both
  pack $drmsd.top -fill both
}

proc driverInterface::doDRMSD {args} {
  variable aselect1 
  variable aselect2
  variable tmpd
  variable switchf

  puts "The tempory directory is $tmpd"

  # Get the list of atoms
  set sel [atomselect top "$aselect1 and $aselect2"]
  set calist [$sel get index]
  #puts "CA list $calist length [llength $calist]" 

  set od [open "$tmpd/plumed.dat" w]
  puts $od "PRINT W_STRIDE 1"
  for {set i 1} { $i<[llength $calist] } { incr i } {
     for {set j 0} { $j<$i } { incr j } { 
         if { $switchf(r0)==0 } {
              puts $od "DISTANCE LIST [expr [lindex $calist $i] + 1] [expr [lindex $calist $j] + 1]"
         } else {
              puts $od "COORD LIST [expr [lindex $calist $i] + 1] [expr [lindex $calist $j] + 1] R_0 $switchf(r0) D_0 $switchf(d0) NN $switchf(nn) MM $switchf(mm)" 
         }
     }
  }
  close $od
  runDriver
  return
}

proc driverInterface::openRDFWindow {args} {
  variable rdf
  variable aselect1
  variable aselect2
  variable lbound
  variable ubound
  variable window
  variable gwidth
  variable tmpd
  variable broadcastvar

  set tmpd [getargs $args "-tmpdir" {}]
  set broadcastvar [getargs $args "-broadcast" {}]

  # Set default values for everything
  set aselect1 "all" ; set aselect2 "all"
  set lbound 0 ; set ubound 5 ; set window 0.25 ; set gwidth 0.125

  if [winfo exists .rdf] {
    wm deiconify $rdf
    return 
  }

  # Open RDF dialogue
  set rdf [toplevel ".rdf"]
  wm title $rdf "RDF CV controls"
  wm resizable $rdf 0 0

  frame $rdf.top -padx 1m -pady 1m
  pack [label $rdf.lab -text "Please select how you would like to calculate the RDF"] -in $rdf.top -side top -fill both
  frame $rdf.atoms -padx 1m -pady 1m
  label $rdf.atoms.lcent -text "Central atom selection"  
  entry $rdf.atoms.ecent -width 10 -textvariable driverInterface::aselect1
  label $rdf.atoms.loth -text "Other atom selection"
  entry $rdf.atoms.eoth  -width 10 -textvariable driverInterface::aselect2
  grid $rdf.atoms.lcent -row 1 -column 1 -sticky e  
  grid $rdf.atoms.ecent -row 1 -column 2 -sticky e
  grid $rdf.atoms.loth  -row 1 -column 3 -sticky e
  grid $rdf.atoms.eoth  -row 1 -column 4 -sticky e
  pack $rdf.atoms -in $rdf.top -side top -fill both

  frame $rdf.param -padx 1m -pady 1m
  label $rdf.param.lrange -text "RDF range"
  entry $rdf.param.lbound -width 5 -textvariable driverInterface::lbound
  entry $rdf.param.ubound -width 5 -textvariable driverInterface::ubound
  label $rdf.param.lwind -text "Window size"
  entry $rdf.param.window -width 5 -textvariable driverInterface::window
  label $rdf.param.lgauss -text "Gaussian width"
  entry $rdf.param.gauss -width 5 -textvariable driverInterface::gwidth
  grid $rdf.param.lrange -row 1 -column 1 -sticky e 
  grid $rdf.param.lbound -row 1 -column 2 -sticky e 
  grid $rdf.param.ubound -row 1 -column 3 -sticky e
  grid $rdf.param.lwind  -row 1 -column 4 -sticky e
  grid $rdf.param.window -row 1 -column 5 -sticky e
  grid $rdf.param.lgauss -row 1 -column 6 -sticky e 
  grid $rdf.param.gauss  -row 1 -column 7 -sticky e
  pack $rdf.param -in $rdf.top -side top -fill both

  # And the buttons
  frame $rdf.buttons -padx 1m -pady 1m
  pack [button $rdf.buttons.cancel -text "cancel" -relief raised -command { destroy .rdf ; return } ] -in $rdf.buttons -side left
  pack [button $rdf.buttons.ok -text "ok" -relief raised -command { destroy .rdf; [namespace code driverInterface::doRDF] } ] -in $rdf.buttons -side right  
  pack $rdf.buttons -in $rdf.top -side top -fill both
  pack $rdf.top -fill both
}

proc driverInterface::doRDF {args} {
  variable aselect1
  variable aselect2
  variable lbound
  variable ubound
  variable window
  variable gwidth
  variable tmpd

  puts "The tempory directory is $tmpd" 
  set od [open "$tmpd/plumed.dat" w]
  puts $od "PRINT W_STRIDE 1"

  # Get the number of bins in the distribution
  set nbins [ expr int( ($ubound - $lbound)/$window ) ]
  if { [ expr $ubound - $lbound - $nbins*$window ] > 0 } { incr nbins }

  # Sort out the atom selection
  if { $aselect1==$aselect2 } {
     set listout "<list1> <list1>" 
     set sel [atomselect top $aselect1]
     set alist [$sel get index]
     puts $od "list1->"  
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     } 
     puts $od "list1<-"
  } else {
     set listout "<list1> <list2>"
     set sel [atomselect top $aselect1]
     set alist [$sel get index]
     puts $od "list1->"  
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }
     puts $od "list1<-"
     set sel [atomselect top $aselect2]
     set alist [$sel get index]
     puts $od "list2->"
     puts $alist  
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }
     puts $od "list2<-" 
  }

  for { set i 0 } { $i < $nbins } { incr i } {
    puts $od "RDF RDF_LABEL 1 LIST $listout RANGE [ expr $lbound + $i*$window ] [ expr $lbound + ($i + 1)*$window ] WIDTH $gwidth"
  }

  close $od
  runDriver   ;  # Now actually run driver 
  return
}

proc driverInterface::openADFWindow {args} {
  variable adf
  variable aselect1
  variable aselect2
  variable aselect3
  variable window
  variable gwidth
  variable switchf
  variable filename
  variable tmpd
  variable broadcastvar  

  set broadcastvar [getargs $args "-broadcast" {}]
  set tmpd [getargs $args "-tmpdir" {}] 

  # Set default values for everything
  set aselect1 "all" ; set aselect2 "all"  ; set aselect3 "all"
  set window 0.175 ; set gwidth 0.09
  set switchf(d0) 0;   set switchf(r0) 3.0
  set switchf(nn) 6;   set switchf(mm) 12

  if [winfo exists .adf] {
     wm deiconify $adf
     return 
  }

  # Open ADF dialogue
  set adf [toplevel ".adf"]
  wm title $adf "ADF CV controls"
  wm resizable $adf 0 0

  frame $adf.top -padx 1m -pady 1m
  pack [label $adf.lab -text "Please select how you would like to calculate the ADF"] -in $adf.top -side top -fill both
  frame $adf.atoms -padx 1m -pady 1m
  label $adf.atoms.lcent -text "Central atom selection"
  entry $adf.atoms.ecent -width 10 -textvariable driverInterface::aselect1
  label $adf.atoms.loth1 -text "atom 1 selection"
  entry $adf.atoms.eoth1  -width 10 -textvariable driverInterface::aselect2
  label $adf.atoms.loth2 -text "atom 2 selection"
  entry $adf.atoms.eoth2  -width 10 -textvariable driverInterface::aselect3
  grid $adf.atoms.lcent -row 1 -column 1 -sticky e
  grid $adf.atoms.ecent -row 1 -column 2 -sticky e
  grid $adf.atoms.loth1  -row 2 -column 1 -sticky e
  grid $adf.atoms.eoth1  -row 2 -column 2 -sticky e
  grid $adf.atoms.loth2  -row 3 -column 1 -sticky e
  grid $adf.atoms.eoth2  -row 3 -column 2 -sticky e
  pack $adf.atoms -in $adf.top -side top -fill both

  labelframe $adf.func -relief ridge -bd 2 -text "switching function" -padx 2m -pady 2m
  pack [ label $adf.func.ld0 -text "d_0" ] -in $adf.func -side left 
  pack [ entry $adf.func.d0 -width 5 -textvariable driverInterface::switchf(d0) ] -in $adf.func -side left
  pack [ label $adf.func.lr0 -text "r_0" ] -in $adf.func -side left
  pack [ entry $adf.func.r0 -width 5 -textvariable driverInterface::switchf(r0) ] -in $adf.func -side left
  pack [ label $adf.func.lnn -text "nn" ] -in $adf.func -side left
  pack [ entry $adf.func.nn -width 5 -textvariable driverInterface::switchf(nn) ] -in $adf.func -side left
  pack [ label $adf.func.lmm -text "mm" ] -in $adf.func -side left
  pack [ entry $adf.func.mm -width 5 -textvariable driverInterface::switchf(mm) ] -in $adf.func -side left
  pack $adf.func -in $adf.top -side top -fill both

  # Controls on the windows
  frame $adf.param -padx 1m -pady 1m
  label $adf.param.lwind -text "Window size"
  entry $adf.param.window -width 5 -textvariable driverInterface::window
  label $adf.param.lgauss -text "Gaussian width"
  entry $adf.param.gauss -width 5 -textvariable driverInterface::gwidth
  grid $adf.param.lwind  -row 1 -column 1 -sticky e
  grid $adf.param.window -row 1 -column 2 -sticky e
  grid $adf.param.lgauss -row 1 -column 3 -sticky e
  grid $adf.param.gauss  -row 1 -column 4 -sticky e
  pack $adf.param -in $adf.top -side top -fill both 

  # Add the buttons
  frame $adf.buttons -padx 1m -pady 1m
  pack [button $adf.buttons.cancel -text "cancel" -relief raised -command { destroy .adf ; return } ] -in $adf.buttons -side left
  pack [button $adf.buttons.ok -text "ok" -relief raised -command { destroy .adf; [namespace code driverInterface::doADF] } ] -in $adf.buttons -side right
  pack $adf.buttons -in $adf.top -side top -fill both
  pack $adf.top -fill both
}

proc driverInterface::doADF {args} {
  variable aselect1
  variable aselect2
  variable aselect3
  variable switchf
  variable window
  variable gwidth
  variable tmpd

  puts "The tempory directory is $tmpd"
  set od [open "$tmpd/plumed.dat" w]
  puts $od "PRINT W_STRIDE 1"

  # Get the number of bins in the distribution
  set nbins [ expr int( 3.142 / $window ) ]
  if { [ expr 3.142 - $nbins*$window ] > 0 } { incr nbins }

  if { $aselect1==$aselect2 && $aselect1==$aselect3 } {
     set listout "<list1> <list1> <list1>"
     set sel [atomselect top $aselect1]
     set alist [$sel get index]
     puts $od "list1->" 
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }
     puts $od "list1<-"
  } elseif { $aselect1==$aselect2 } {
     set listout "<list1> <list1> <list2>"
     set sel [atomselect top $aselect1]
     set alist [$sel get index]
     puts $od "list1->" 
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }
     puts $od "list1<-"
     set sel [atomselect top $aselect3]
     set alist [$sel get index]
     puts $od "list2->"
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }
     puts $od "list2<-"
  } elseif { $aselect1==$aselect3 } {
     set listout "<list1> <list2> <list1>"
     set sel [atomselect top $aselect1]
     set alist [$sel get index]
     puts $od "list1->"
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }
     puts $od "list1<-"
     set sel [atomselect top $aselect2]
     set alist [$sel get index]
     puts $od "list2->"
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }   
     puts $od "list2<-"
  } elseif { $aselect2==$aselect3 } {
     set listout "<list1> <list2> <list2>"
     set sel [atomselect top $aselect1]
     set alist [$sel get index]
     puts $od "list1->"
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }
     puts $od "list1<-"
     set sel [atomselect top $aselect2]
     set alist [$sel get index]
     puts $od "list2->"
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }   
     puts $od "list2<-"
  } else {
     set listout "<list1> <list2> <list3>"
     set sel [atomselect top $aselect1]
     set alist [$sel get index]
     puts $od "list1->"
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }
     puts $od "list1<-"
     set sel [atomselect top $aselect2]
     set alist [$sel get index]
     puts $od "list2->"
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }
     puts $od "list2<-"
          set sel [atomselect top $aselect3]
     set alist [$sel get index]
     puts $od "list3->"
     for { set i 0 } { $i<[llength $alist] } { incr i } {
         puts $od "[expr [lindex $alist $i] + 1]"
     }
     puts $od "list3<-"
  }

  for { set i 0 } { $i < $nbins } { incr i } {
    puts $od "ADF RDF_LABEL 1 LIST $listout RANGE [ expr $i*$window ] [ expr ($i + 1)*$window ] WIDTH $gwidth R_0 $switchf(r0) D_0 $switchf(d0) NN $switchf(nn) MM $switchf(mm)"
  }
  close $od
  runDriver   ;  # Now actually run driver 
  return
}

proc driverInterface::getargs {arglist tag defl {n 1}} {
   set pos [lsearch $arglist $tag]
   if {$pos<0}  { return $defl}
   return [join [ lrange $arglist [expr $pos+$n ] [expr $pos+$n ] ] ]
}

proc driverInterface::readPlumed {args} {
   variable filename
   variable tfilename
   variable tmpd
   variable lenunit 

   # open a file and setup the tmpdir to which we are copying the file
   set fd [open $filename r]
   set od [open "$tmpd/plumed.dat" w]

   # We must print colvars at every trajectory step
   puts $od "PRINT W_STRIDE 1"

   while { [gets $fd line] != -1 } {
      set sline [regexp -inline -all -- {\S+} $line]

      # Deal with CV units and extra files
      if { !([lsearch $sline "COORD"] < 0) || !([lsearch $sline "HBONDS"] < 0) || \
           !([lsearch $sline "WATERBRIDGE"] < 0) || !([lsearch $sline "ELSTPOT"] < 0) } {
         set flag 0
         set newcomm {}
         foreach elem $sline {
            if { $flag==1 } {
               lappend newcomm [expr $elem * $lenunit]
               set flag 0
            } else {
               lappend newcomm $elem
               if { $elem=="R_0" || $elem=="D_0" } { set flag 1 }
            }
         }
         puts $od [join $newcomm]
      } elseif { ![lsearch $sline "ENERGY"] < 0 } {
         tk_messageBox -icon error -type ok -title Message -message "It is not possible to calculate the energy using driver"
         return
      } elseif { !([lsearch $sline "ALPHARMSD"] < 0) || !([lsearch $sline "ANTIBETARMSD"] < 0) || !([lsearch $sline "PARABETARMSD"] < 0) } {
          set flag 0
          set newcomm {}
          foreach elem $sline {
            if { $flag==2 } {
               lappend newcomm [expr $elem * $lenunit]
               set flag 0
            } elseif { $flag==1 } {
               set flag 0
            } else {
               if { $elem=="ANGSTROM SCALE" } {
                  set flag 1
               } elseif { $elem=="R_0" } {
                  set flag 2
                  lappend newcomm $ele
               } else {
                  lappend newcomm $ele
               }
            }
          }
          puts $od [join $newcomm]
      } elseif { !([lsearch $sline "S_PATH"] < 0) || !([lsearch $sline "Z_PATH"] < 0) || !([lsearch $sline "TARGETED"] < 0) } {
          if { !([lsearch -all $sline "CMAP"] < 0) } {
             if { $lenunit==1.0 } {
                set pos [lsearch $sline "INDEX"]
                set fname [join [ lindex $sline [expr $pos+1 ] ] ]
                foreach filen [glob $fname*] { file copy -force $filen $tmpd/$filen }
                set pos [lsearch $sline "MAP"]
                set fname [join [ lindex $sline [expr $pos+1 ] ] ]
                foreach filen [glob $fname*] { file copy -force $filen $tmpd/$filen } 
                set pos [lsearch $sline "GROUP"]
                if { !($pos < 0) } {
                   set fname [join [ lindex $sline [expr $pos+1 ] ] ]
                   foreach filen [glob $fname*] { file copy -force $filen $tmpd/$filen }
                }
             } else {
                tk_messageBox -icon error -type ok -title Message -message "Cannot interpret PATH CV with CMAP if plumed.dat file is not in Angstroms"
                return
             }
             puts $od $line
          } else {
             set pos [lsearch $sline "FRAMESET"]
             set fname [join [ lindex $sline [expr $pos+1 ] ] ]
             foreach filen [glob $fname*] { file copy -force $filen $tmpd/$filen }
             set flag 0
             set newcomm {}
             foreach elem $sline {
               if { $flag==1 } { 
                  lappend newcomm [expr $elem / ( $lenunit * $lenunit ) ]
                  set flag 0
               } else {
                  lappend newcomm $elem
                  if { $elem=="LAMBDA" } { set flag 1 }
               }
             }
             puts $od [join $newcomm]
          }
      } elseif { !([lsearch $sline "BESPOKE"] < 0) } {
         # GAT note the problems that might occur with distance cvs and this if
         #     we are not careful and ensure to take this into account when we
         #     develop this further
         set pos [ lsearch $sline "FILENAME" ]
         if { $pos>0 } {
            set fname [join [ lindex $sline [expr $pos+1 ] ] ]
            file copy -force $fname $tmpd
         } else {
            tk_messageBox -icon error -type ok -title Message -message "Syntax for bespoke cv is wrong - no filename"
            return
         }
         puts $od $line
      } elseif { !([lsearch $sline "PCA"] < 0) } {
         set pos [lsearch  $sline "FRAME" ]
         if { $pos>0 } {
            set fname [join [ lindex $sline [expr $pos+1 ] ] ]
            file copy -force $fname $tmpd
         } else {
            tk_messageBox -icon error -type ok -title Message -message "Syntax for pca cv is wrong - no frame"
            return
         }
         set pos [lsearch  $sline "EIGENVEC" ]
         if { $pos>0 } {
            set fname [join [ lindex $sline [expr $pos+1 ] ] ]
            file copy -force $fname $tmpd
         } else {
            tk_messageBox -icon error -type ok -title Message -message "Syntax for pca cv is wrong - no eigenvec"
            return
         }
      } elseif { !([lsearch $sline "RDF"] < 0) || !([lsearch $sline "ZDIST"] <0) } {
         set flag 0
         set newcomm {}
         set bounds {}
         foreach elem $sline {
            if { $flag!=0 } {
               lappend newcomm [expr $elem * $lenunit]
               if { $flag==2 } {
                  lappend bounds [expr $elem * $lenunit]
               } elseif { $flag==1 && [ llength $bounds ]==1 } {
                  lappend bounds [expr $elem * $lenunit]
               }
               set flag [expr $flag-1]
            } else {
               lappend newcomm $elem
               if { $elem=="RANGE" } {
                  set flag 2
               } elseif { $elem=="WIDTH" } {
                  set flag 1
               }
            }
         }
         puts $od [join $newcomm]
      } elseif { [lsearch $sline "PRINT"]<0 && [lsearch $sline "RECONNAISSANCE"]<0 && \
           [lsearch $sline "ONIONS"]<0 && [lsearch $sline "BASINS"]<0 && \
           [lsearch $sline "CLUSTER"]<0 && [lsearch $sline "HILLS"]<0 && \
           [lsearch $sline "WELLTEMPERED"]<0 } {
         # CVS dealt with here  : RGYR    DISTANCE  TORSION    POSITION    MINDIST    ANGLE
         #                        DIPOLE  DIHCOR    ALPHABETA  RMSDTOR     PUCKERING  HELIX
         puts $od $line
         # GAT note that syntax for walls etc is not OK  -- this should be dealt with perhaps?
      }
      # GAT missing CVs   ::  CMAP
   }
   close $fd   ; close $od    ; # Close the files 
  
   # Copy the topology file if it is setup
   if { $tfilename!="none" && $tfilename!="" } { file copy -force $tfilename $tmpd/topol.pdb }

   runDriver   ; # Now run driver
   return
}

proc driverInterface::runDriver {args} {
  variable tmpd
  variable status
  variable broadcastvar

  global env

  set owd [pwd] ; cd $tmpd

  # Check that driver is where it is supposed to be 
  if { ![file exists "$env(plumedir)/utilities/driver/driver"] } {
     tk_messageBox -icon error -type ok -title Message -message "Could not find driver should be in $env(plumedir)/utilities/driver/driver"
     cd $owd
     return
  }

  # We must write out a pdb file for analysis
  if { ![file exists topol.pdb] } {
    tk_messageBox -icon error -type ok -title Message -message "Warning - you did not specify a pdb file containing masses and charges \n please check masses and charges in $tmpd/topol.pdb"
    [atomselect top all] set occupancy [ [atomselect top all] get mass ]
    [atomselect top all] set beta [ [atomselect top all] get charge ]
    [atomselect top all] writepdb topol.pdb
  }

  # Write out a dcd file for the whole trajectory
  animate write dcd temp.dcd waitfor all top

  # Get the periodic boundary conditions from vmd or by asking user
  set pbc [PeriodicBoundaryCondition]
  if { [lsearch $pbc -cell]<0 && [lsearch $pbc -nopbc]<0 } {
      tk_messageBox -icon warning -type ok -title Message -message "Sorry I couldn't find pbcs inside your input trajectory \n I have output all the file to run driver in $tmpd \n
 Please run driver manually and load resulting COLVAR file"
      cd $owd
      return
  }
  # Now run driver
  puts "Executing: $env(plumedir)/utilities/driver/driver -dcd temp.dcd -pdb topol.pdb -plumed plumed.dat $pbc"
  catch {eval exec $env(plumedir)/utilities/driver/driver -dcd temp.dcd -pdb topol.pdb -plumed plumed.dat $pbc} driver_stdout

  # Check it worked
  if { ![file exists COLVAR] } {
     puts $driver_stdout
     puts "-------------"
     tk_messageBox -icon error -type ok -title Message -message "Something went wrong when running driver.  Please see messages in log and files in tempory directory [pwd]"
     cd $owd
     return
  }
  cd $owd
  set status($broadcastvar) 1    ; # This broadcasts to whoever is interested that we have new cv data
  return
}

proc driverInterface::PeriodicBoundaryCondition {args} {
  if { [molinfo top get alpha]!=90 || [molinfo top get beta]!=90 || [molinfo top get gamma]!=90 } {
     tk_messageBox -icon error -type ok -title Message -message "Can't run driver with non-orthorhombic unit cells"       
     return -1
  }   
  if { [molinfo top get a]==0 } {
     set ans [tk_messageBox -icon question -type yesno -title Message -message "Are there really no pbcs?"]
     switch -- $ans {
         yes { return -nopbc }
         no  { return -1 }
     } 
  } else {
     return "-cell [molinfo top get a] [molinfo top get b] [molinfo top get c]"
  }
}
