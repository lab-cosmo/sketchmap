#
# Graphical user interface for PLUMED
#
# $Id: plumed_gui.tcl,v 1.0 2010/11/15
# 

# To do list
# (1) Need to do a lot of testing of using gui in tandem with multiple molecules in vmd (this needs to be more robust methinks)
# (2) It would be nice if this could recognise when new molecules are loaded and if it could add an element to cvdata([list $molid havecvs])
# (3) Separate reconnaissance tools to separate namespace?
# (4) We need to be ablet to recognise which molecule we have plotted data for
# (5) Draw all colvars doesn't behave well with axis adjustment tools
# (6) We want all other windows to be destroyed if we shut down the gui
# (9) Control over when we rescale axis is very unsatisfactory - is there a neater way to do this?
# (10) Make read fes work with files that are not called fes.dat - i.e. we want to rename all the files we copy fes.dat
# (11) Note to self molecule changes are still a bit of a mess

source $env(VMDDIR)/scripts/init.d/driver_interface.tcl
source $env(VMDDIR)/scripts/init.d/colvarTools.tcl
source $env(VMDDIR)/scripts/init.d/cvlist.tcl

package provide plumed_gui 2.0
package require cvlist 1.0
package require driverInterface 1.0
package require cv_traj 1.0
package require gtPlot 2.0
package require fesTools 1.0

vmd_install_extension plumed_gui plumedVis_tk "Analysis/GISMO"

proc plumedVis_tk {} {
    ::plumedVis::create_gui
    return $plumedVis::w
}

namespace eval ::plumedVis {
  variable w                    ;# handle to window
  variable menumol              ;# the molecule that is currently being used for the menus
  variable cvwindow             ;# handle to cvlook window
  variable feswindow            ;# handle to feslook window
  variable axiswindow           ;# handle to axis properties window
  variable fps                  ;# handle to window that controls all landmark point extraction stuff
  variable dimred               ;# handle to dimensionality reduction window
  variable osample              ;# handle to out of sample embedding window
  variable pfilew               ;# handle to projection file creating window
  variable appearance           ;# Various controls over the way things look
  variable filename             ;# filename to read in
  variable fcoord               ;# What free energy surface are we showing
  variable xcoord               ;# The coord on the x axis
  variable ycoord               ;# The coord on the y axis
  variable scoord               ;# The size of the points
  variable ccoord               ;# The color of the points
  variable showfesdiff          ;# Are we showing the difference from the final fes?
  variable hideColvars          ;# Are we showing the basins
  variable allcolvars           ;# A binary flag that tells us if we are plotting all colvars or just one
  variable fesdisplay           ;# The number of the currently displayed fes
}

proc plumedVis::setDefaults {args} {
   variable appearance
   variable fesdisplay
   variable allcolvars
   set appearance(canvasSize)     400           ; #The size of the canvas
   set appearance(canvasPadding)  60            ; #The ammount of white space around the canvas
   set appearance(font)          "Helvetica 8"  ; #Font for axis labels and so on 
   set showfesdiff         0                    ; #We are showing the fes not the difference
   set hideColvars         0                    ; #We are showing the colvar positions
   set fesdisplay           1                   ; #Initialize all stuff for free energy surfaces 
   set allcolvars          0                    ; #We are not plotting all the colvars by default               
   ::cv_traj::setDefaults -datatag colvar       ; #Set defaults for stuff in cv_traj
}

proc plumedVis::create_gui {args} {
  # Just create the window and initialize data structures
  # No molecule has been selected yet
  # Also set up traces on VMD variables

  variable w
  variable appearance
  variable xcoord
  variable ycoord
  variable fesdisplay

  global vmd_frame
  global env

  # Check plumedir is defined (if this is not the case we should crash
  if { $env(plumedir) == "" } { tk_messageBox -icon error -type ok -title Message -message "I have not found plumedir in env, update .bashrc" }

  # If already initialized, destroy it
  if [winfo exists .bplot] {
    wm deiconify $w
    return
  }

  set w [toplevel ".bplot"]
  wm title $w "G.I.S.M.O - Graphical Interface for Sketch-Map"
  wm resizable $w 1 1
  bind $w <Destroy> plumedVis::destroy

  # Set all the defaults
  setDefaults

  # Create a tempory directory in which to do calculations 
  set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]"
  file mkdir $tmpd 

  frame $w.top -padx 1m -pady 1m

  pack [ createMenubar $w.top.menubar ] -padx 1 -fill x -side top

  # This is our graph object
  frame $w.fr
  gtPlot::gtcreate $w.fr.canv -size $appearance(canvasSize) -font $appearance(font) -padding $appearance(canvasPadding)
  grid $w.fr.canv -sticky news
  grid rowconfigure $w.fr 0 -weight 1
  grid columnconfigure $w.fr 0 -weight 1
  pack $w.fr -in $w.top -fill both -expand true -side top

  # Add axis controls to appearance menu (do we really want this here GAT)  
  $w.top.menubar.look.menu add checkbutton -label "Fix z scale" -variable "::gtPlot::bounds(fixedz)"

  pack [ createControls $w.top.pmenus ] -in $w.top -padx 1 -fill x -side top
  pack $w.top -fill both -expand true                      ; # Complete the packing
  updateMenus                                              ; # Sets up all the menus and sets what's plotted to none
  updateFesMenu                                            ; # Set up the fes menu

  # This should bind every colvar point to go to the frame
  $w.fr.canv bind colvar: <Button-1> [namespace code { cv_traj::Goto %W -datatag colvar } ]
  $w.fr.canv bind fps: <Button-1> [namespace code { cv_traj::Goto %W -datatag fps } ] 
 
  # This stuff allows us to select points 
  trace add variable ::gtPlot::doneLasso write [namespace code { cv_traj::addSelection $w.fr.canv -datatag colvar} ] 
  # This stuff allows us to integrate free energies
  bind $w.fr.canv <Button-1> { ::cv_traj::clearSelection %W -datatag colvar ; ::fesTools::clearSelection %W }
  trace add variable ::gtPlot::doneLasso write [namespace code { fesTools::integrate $w.fr.canv -datatag free_energy} ]

  # This controls the plotter
  trace add variable plumedVis::xcoord write [namespace code plumedVis::updateColvars ] 
  trace add variable plumedVis::ycoord write [namespace code plumedVis::updateColvars ]
  trace add variable plumedVis::scoord write [namespace code plumedVis::updateColvars ]
  trace add variable plumedVis::ccoord write [namespace code plumedVis::updateColvars ]
  trace add variable plumedVis::fcoord write [namespace code plumedVis::updateFes]
  trace add variable plumedVis::showfesdiff write [namespace code plumedVis::drawGraph]
  trace add variable plumedVis::fesdisplay write [namespace code plumedVis::drawGraph]
  trace add variable plumedVis::hideColvars write [namespace code plumedVis::updateColvars ]
  trace add variable vmd_frame write [namespace code updateHighlight]
  # This ensures colvars are read when they are calculated
  trace add variable ::driverInterface::status(plumedVis) write [namespace code plumedVis::readColvar]
  # This ensures fes is plotted after sum hills
  trace add execution ::fesTools::sumHills leave [namespace code plumedVis::updateFesMenu] 
  # This kills the dimensionality reduction window when it needs to be killed
  trace add variable ::cvlist::dimredFinished write {destroy .dimred}
  # This kills the osample window when it needs to be killed
  trace add variable ::cvlist::osampleFinished write {destroy .osample}
  # This kills the projection file window when it needs to be killed
  trace add variable ::cvlist::projectionFileFinished write {destroy .pfilew}
  # This makes sure menus are updated when we have new dimensionality reduction data
  trace add execution ::cvlist::addDataCol leave [namespace code plumedVis::updateMenus]
}

proc plumedVis::destroy {args} {
   variable allcolvars
   variable xcoord
   variable ycoord

   global vmd_frame

   # Unset the plots (hopefully)
   set allcolvars 0
   set xcoord "none"

   # Remove the traces   
   trace remove variable plumedVis::xcoord write [namespace code plumedVis::updateColvars ] 
   trace remove variable plumedVis::ycoord write [namespace code plumedVis::updateColvars ] 
   trace remove variable plumedVis::scoord write [namespace code plumedVis::updateColvars ]
   trace remove variable plumedVis::ccoord write [namespace code plumedVis::updateColvars ]
   trace remove variable plumedVis::hideColvars write [namespace code plumedVis::updateColvars ]
   trace remove variable plumedVis::fcoord write [namespace code plumedVis::updateFes]
   trace remove variable plumedVis::showfesdiff write [namespace code plumedVis::drawGraph] 
   trace remove variable plumedVis::fesdisplay write [namespace code plumedVis::drawGraph]
   trace remove variable vmd_frame write [namespace code updateHighlight]
   trace remove variable ::driverInterface::status(plumedVis) write [namespace code plumedVis::readColvar]
   trace remove variable ::gtPlot::doneLasso write [namespace code { cv_traj::addSelection $w.fr.canv -datatag colvar} ]
   trace remove variable ::gtPlot::doneLasso write [namespace code { fesTools::integrate $w.fr.canv -datatag free_energy} ]
   trace remove execution ::fesTools::sumHills leave [namespace code plumedVis::updateFesMenu] 
   trace remove variable ::cvlist::dimredFinished write {destroy .dimred}
   trace remove variable ::cvlist::osampleFinished write {destroy .osample}
   trace remove variable ::cvlist::projectionFileFinished write {destroy .pfilew}
   trace remove execution ::cvlist::addDataCol leave [namespace code plumedVis::updateMenus]
}

proc plumedVis::initosample {args} {
  variable osample

  if { ![ ::cvlist::exists -datatag pgui ] } { tk_messageBox -icon error -type ok -title Message -message "Can't do dimensionality reduction out of sample embedding without cv data" ; return }

  set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]/osample"
  file mkdir $tmpd

  if [winfo exists .osample] {
     wm deiconify $osample
     return
  }

  set osample [toplevel ".osample"]
  wm title $osample "Out of sample embedding controls"
  wm resizable $osample 0 0
  bind $osample <Destroy> cvlist::destroyOutOfSampleWindow

  frame $osample.top -padx 1m -pady 1m
  pack [ ::cvlist::outOfSampleWindow $osample.top.fr -datatag pgui -tmpd $tmpd ] -in $osample.top -side top -fill both
  # And the buttons to make it run
  frame $osample.top.buttons -padx 1m -pady 1m
  pack [ button $osample.top.buttons.cancel -text "cancel" -relief raised -command {destroy .osample} ] -in $osample.top.buttons -side left
  pack [ button $osample.top.buttons.ok -text "run" -relief raised -command [namespace code { cvlist::outOfSampleEmbedding -datatag pgui } ] ] -in $osample.top.buttons -side right
  pack $osample.top.buttons -in $osample.top -side top -fill both -expand 1
  pack $osample.top -fill both
}

proc plumedVis::initdimred {args} {
  variable dimred

  if { ![ ::cvlist::exists -datatag pgui ] } { tk_messageBox -icon error -type ok -title Message -message "Can't do dimensionality reduction without cv data" ; return }

  set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]/dimred"
  file mkdir $tmpd

  if [winfo exists .dimred] {
     wm deiconify $dimred
     return
  }

  set dimred [toplevel ".dimred"]
  wm title $dimred "Dimensionality reduction controls"
  wm resizable $dimred 0 0

  frame $dimred.top -padx 1m -pady 1m
  pack [ ::cvlist::dimensionalityReductionWindow $dimred.top.fr -datatag pgui -tmpd $tmpd ] -in $dimred.top -side top -fill both
  pack $dimred.top -fill both
}

proc plumedVis::initfps {args} {
   variable fps
   if { ![ ::cvlist::exists -datatag pgui ] } { tk_messageBox -icon error -type ok -title Message -message "Can't extract landmarks without cv data" ; return }

   set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]/FPS"     
   file mkdir $tmpd

   if [winfo exists .fps] {
     wm deiconify $fps
     return
   }

   set fps [toplevel ".fps"]
   wm title $fps "Get landmarks controls"
   wm resizable $fps 0 0
 
   frame $fps.top -padx 1m -pady 1m
   pack [ ::cvlist::landmarkWindow $fps.top.fr -datatag pgui -tmpd $tmpd ] -in $fps.top -side top -fill both
   frame $fps.d -padx 3m -pady 3m
   pack [ button $fps.d.ok -text "Show landmarks" -relief raised -command [namespace code plumedVis::plotLandmarks] ] -in $fps.d -side left
   pack [ button $fps.d.get -text "Store landmarks" -relief raised -command [namespace code plumedVis::storeLandmarks] ] -in $fps.d -side left
   pack [ button $fps.d.dismiss -text "Dismiss" -relief raised -command { destroy .fps} ] -in $fps.d -side right
   pack $fps.d -in $fps.top -side top
   pack $fps.top -fill both
}

proc plumedVis::plotLandmarks {args} {
  variable w
  variable xcoord
  variable ycoord
  variable hideColvars

  set landmarks [ ::cvlist::createLandmarks -datatag pgui ]
  if { [llength $landmarks]==0 } { return }
  if { $xcoord!="none" && $ycoord!="none" } { 
     ::cv_traj::clearSelection $w.fr.canv -datatag colvar
     ::cv_traj::createSelection $w.fr.canv -datatag colvar -selection $landmarks -color "blue"
  } else {
     tk_messageBox -icon error -type ok -title Message -message "You cannot plot landmarks unless data is plotted.  Set x and y coordinate menus and try again"
     return
  }
}

proc plumedVis::storeLandmarks {args} {
  set landmarks [ ::cvlist::createLandmarks -datatag pgui ]
  if { [llength $landmarks ]==0 } { return }
  set newmol [ storeSelection -selection $landmarks ]
  ::cvlist::addDataCol -datatag pgui -molecule $newmol -name "weights" -data [ ::cvlist::getLandmarkWeights -datatag pgui ]
}

proc plumedVis::openProjectionFileWindow {args} {
  variable pfilew
  if { ![ ::cvlist::exists -datatag pgui ] } { 
     tk_messageBox -icon error -type ok -title Message -message "Can't create a projection withou cv data" 
     return 
  }

  if [winfo exists .pfilew] {
     wm deiconify $pfielw
     return
  }

  set pfilew [toplevel ".pfilew"]
  wm title $pfilew "Save projection file window"
  wm resizable $pfilew 0 0
  bind $pfilew <Destroy> cvlist::destroyProjectionFileWindow

  frame $pfilew.top -padx 1m -pady 1m
  pack [ ::cvlist::projectionFileWindow $pfilew.top.fr -datatag pgui ] -in $pfilew.top -side top -fill both
  
  # These are the buttons
  frame $pfilew.d -padx 3m -pady 3m
  pack [ button $pfilew.d.cancel -text "Cancel" -relief raised -command { destroy .pfilew } ] -in $pfilew.d -side left
  pack [ button $pfilew.d.ok -text "OK" -relief raised -command [namespace code {::cvlist::printProjectionFile -datatag pgui}] ] -in $pfilew.d -side right
  pack $pfilew.d -in $pfilew.top -side top
  pack $pfilew.top -fill both
}

proc plumedVis::updateHighlight {args} {
  variable w
  variable menumol
  variable xcoord
  variable ycoord
  variable hideColvars
  variable allcolvars

  # Create new menus if molecule has changed (also deletes the plot) 
  set molid [findMolecule]
  if { $molid!=$menumol } { updateMenus ; return }

  # Again we need something to tell us what molecule is plotted -- GAT

  # If we are in allcv plotting mode then we just redraw the whole bar chart.
  if { $allcolvars==1 } { plotAllCVs ; return } 
  
  # This does tne update for colvar highlights
  if { $hideColvars!=1 && $xcoord!="none" && $ycoord!="none" } { ::cv_traj::updateHighlight $w.fr.canv -datatag colvar }     
}

# This plots a bar chart of the cvs in the current frame
proc plumedVis::plotAllCVs {args} {
   variable w
   variable allcolvars
   variable xcoord
   variable ycoord
   variable ccoord
   variable scoord

   global vmd_frame

   # Get rid of any plots of cvs in 2d
   if { $xcoord!="none" } { set xcoord "none" }
   if { $ycoord!="none" } { set ycoord "none" }
   if { $ccoord!="none" } { set ccoord "none" }
   if { $scoord!="none" } { set scoord "none" }
   # Get rid of any free energy surface from the plot
   ::gtPlot::unplotData $w.fr.canv -datatag free_energy

   if { ![ ::cvlist::exists -datatag pgui] } { return }

   # Extract the colvar data for this particular frame
   set molid [findMolecule]
   set ydata [::cvlist::getDataRow -datatag pgui -row $vmd_frame($molid) ]
   set n 0 
   foreach item $ydata { lappend xdata $n ; incr n }

   # Remove the old allcolvars data
   ::gtPlot::unplotData $w.fr.canv -datatag allcolvars
   # Pass data to gtplot
   ::gtPlot::readData -datatag allcolvars -npoints [llength $ydata] -xdata $xdata -ydata $ydata -colors "black" 
   # Draw the graph
   ::gtPlot::plotData $w.fr.canv -datatag allcolvars -plotstyle "bar" -stride 1 -linecolor black 

   set allcolvars 1   ;  # This tells all other routines that all the colvars are plotted
}

proc plumedVis::drawGraph {args} {
   variable w
   variable xcoord
   variable ycoord
   variable fcoord
   variable fesdisplay
   variable hideColvars
   variable showfesdiff

   # This should also work with allcolvar mode as in it currently doesn't
   if { $fcoord!="none" } {
      if { $showfesdiff==1 } {
        ::fesTools::plotFesDiff $w.fr.canv -directory $fcoord -datatag free_energy -number $fesdisplay  
      } else {
        ::fesTools::plotFes $w.fr.canv -directory $fcoord  -datatag free_energy -number $fesdisplay   
      }  
   }
   # Use the traj plotter to plot the fps data
   if { $xcoord!="none" && $ycoord!="none" } {
     # Use the traj plotter to plot the data 
     if { $hideColvars!=1 } { ::cv_traj::plotData $w.fr.canv -datatag colvar }
   }
}

proc plumedVis::updateFes {args} {
  variable fesdisplay
  variable xcoord
  variable ycoord
  variable fcoord
  variable allcolvars
  variable w  

  # Delete the bar charts of all the colvars 
  ::gtPlot::unplotData $w.fr.canv -datatag allcolvars ; set allcolvars 0

  # Add something to delete the fes here GAT
  if { $fcoord=="none" } { 
    ::fesTools::unplotFes $w.fr.canv -datatag free_energy
    return 
  }
  ::fesTools::configureSpinbox $w.top.pmenus.fram.spinb -directory $fcoord -textvariable plumedVis::fesdisplay
}

proc plumedVis::updateColvars {args} {
   variable xcoord
   variable ycoord
   variable scoord
   variable ccoord
   variable fcoord
   variable allcolvars
   variable w

   # Update the fes menu
   ::fesTools::updateMenu $w.top.pmenus.faxis.m -xcoord $xcoord -ycoord $ycoord -variable plumedVis::fcoord

   # Delete the bar charts of all the colvars 
   ::gtPlot::unplotData $w.fr.canv -datatag allcolvars ; set allcolvars 0

   # Delete the old datapoints
   ::cv_traj::delete $w.fr.canv -datatag colvar
   ::cv_traj::delete $w.fr.canv -datatag fps

   # Delete any old free energy surfaces
   set fcoord "none"

   # Return if we are not plotting anything
   if { $xcoord!="none" && $ycoord!="none" } { 
      # Check that there is data for the molecule on top
      if { ![ ::cvlist::exists -datatag pgui ] } {
         updateMenus
         tk_messageBox -icon error -type ok -title Message -message "No cv information available for this molecule"
         return
      }
      set xdata [ ::cvlist::getDataCol -datatag pgui -name $xcoord ]
      set ydata [ ::cvlist::getDataCol -datatag pgui -name $ycoord ]
      if { $xdata=="error" || $ydata=="error" } {
         updateMenus; tk_messageBox -icon error -type ok -title Message -message "No cv information available or number of colvar values does not match number of frames for the molecule"; return;
      }

      # Get color data
      if { $ccoord!="none" } {
         set cdata [ ::cvlist::getDataCol -datatag pgui -name $ccoord ]   
      } else {
         set cdata "white"
      }
      
      # Get point size data
      if { $scoord!="none" } {
         set weights [ ::cvlist::getDataCol -datatag pgui -name $scoord ] 
      } else {
         for { set i 0 } { $i< [molinfo top get numframes] } { incr i } { lappend weights 1 }
      }
      # Set up the frames
      for { set i 0 } { $i< [molinfo top get numframes] } { incr i } { lappend frames $i }

      # Pass the data to the traj plotter
      ::cv_traj::readData -datatag colvar -npoints [molinfo top get numframes] -xdata $xdata -ydata $ydata -cdata $cdata -weights $weights -frames $frames
   }

   # Redraw the graph
   drawGraph
}

proc plumedVis::updateFesMenu {args} {
   variable xcoord
   variable ycoord
   variable fcoord
   variable w
   ::fesTools::updateMenu $w.top.pmenus.faxis.m -xcoord $xcoord -ycoord $ycoord -variable plumedVis::fcoord
}

proc plumedVis::storeSelection {args} {
   variable w
   variable xcoord
   variable ycoord
   variable recondata 

   global vmd_frame

   set frames [getargs $args "-selection" 0]
   if { [llength $frames]==1 && $frames==0 } {
      if { $xcoord=="none" || $ycoord=="none" } { return }
      set frames [ ::cv_traj::getSelection -datatag colvar ] 
      ::cv_traj::clearSelection $w.fr.canv -datatag colvar
   } 

   if { ![ ::cvlist::exists -datatag pgui ] } { return }
 
   # Find the tmp directory in which to do i/o
   set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]"

   # Copy the cvdata (and coordinates) of the old molecule to the new molecule 
   set newmol [ ::cvlist::clone -tmpdir $tmpd -gtag pgui -selection $frames ]
   # Make sure the menus are up to date
   updateMenus

   return $newmol
}

# These routines calculate various kinds of colvars
proc plumedVis::getColvars {type} {
   set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]"

   if { $type=="torsions" } {
       ::driverInterface::openTorsionsWindow -tmpdir $tmpd -broadcast "plumedVis"
   } elseif { $type=="drmsd" } {
       ::driverInterface::openDRMSDWindow -tmpdir $tmpd -broadcast "plumedVis"
   } elseif { $type=="rdf" } {
       ::driverInterface::openRDFWindow -tmpdir $tmpd -broadcast "plumedVis"
   } elseif { $type=="adf" } {
       ::driverInterface::openADFWindow -tmpdir $tmpd -broadcast "plumedVis"
   } else {
       tk_messageBox -icon error -type ok -title Message -message "There is no colvar shortcut of type $type"
   }
}

# These routines read various types of file
proc plumedVis::browse {todo title} {
   variable w
   variable filename
   variable fcoord
   variable xcoord
   variable ycoord

   set tmp [ tk_getOpenFile -initialdir [pwd] -title $title ]

   if { $todo=="readcolvars" && $tmp != "" } {
      set filename $tmp 
      readColvar 
   } elseif { $todo=="calccolvars" && $tmp != "" } {
      set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]"
      set filename "$tmpd/COLVAR"
      ::driverInterface::openLenUnitDialogue -filename $tmp -tmpdir $tmpd -broadcast "plumedVis"
   } elseif { $todo=="readfes" && $tmp != "" } {
      set filename $tmp
      set fesstat [::fesTools::fesCreator -style "readfes" -tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]" -filen $filename]
      if { $fesstat=="error" } { return }
      set fcoord $fesstat
      ::fesTools::updateMenu $w.top.pmenus.faxis.m -xcoord $xcoord -ycoord $ycoord -variable plumedVis::fcoord
   } elseif { $todo=="findbasins" && $tmp != "" } {
      set filename $tmp
      openBasinsDialogue
   } elseif { $todo=="sumhills" && $tmp != "" } {
      set filename $tmp
      ::fesTools::fesCreator -style "sumhills" -tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]" -filen $filename 
   } elseif { $todo=="readonions" && $tmp != "" } {
      set filename $tmp
      readOnions
   } else {
      puts "No mode specified for after in browse"
      return
   }
}

proc plumedVis::drawConvergencePlot {args} {
   variable fcoord

   if { $fcoord=="none" } {
     tk_messageBox -icon error -type ok -title Message -message "There is no free energy surface plotted"
     return
   }
   ::fesTools::drawConvergencePlot -directory $fcoord                   
}

proc plumedVis::readColvar {args} {
   variable xcoord
   variable ycoord
   variable filename
   ::cvlist::readColvar -filename $filename -datatag pgui
}

proc plumedVis::reset {args} {
# This will delete all the data in the cvlist, 
# all free energy surfaces and all plotted data.
# It only does this for the current top molecule though
  variable xcoord
  variable ycoord
  variable ccoord
  variable sccord

  # First delete the plot
  set xcoord "none" ; set ycoord "none" 
  set ccoord "none" ; set scoord "none"

  # Now delete all free energy surface data for this molecule
  ::fesTools::deleteTopMolData 

  # Now delete everything from the cvlists
  ::cvlist::deleteTopMolData -datatag pgui  

  # And update the menus
  updateMenus
}

proc plumedVis::updateMenus {args} {
   variable menumol
   variable xcoord
   variable ycoord
   variable ccoord
   variable scoord
   variable w

   set menumol [findMolecule]
   ::cvlist::updateMenu $w.top.pmenus.xaxis.m -datatag pgui -variable plumedVis::xcoord
   ::cvlist::updateMenu $w.top.pmenus.yaxis.m -datatag pgui -variable plumedVis::ycoord 
   ::cvlist::updateMenu $w.top.pmenus.saxis.m -datatag pgui -variable plumedVis::scoord 
   ::cvlist::updateMenu $w.top.pmenus.caxis.m -datatag pgui -variable plumedVis::ccoord
}

proc plumedVis::printCanvas {args} {
# This retarded routine is necessary as I can't pass the window through the button
  variable w
  ::gtPlot::printCanvas $w.fr.canv
}

# These are the tools for drawing the parts of the window
proc plumedVis::createMenubar {bar} {
  frame $bar -relief raised -bd 2

  # File menu
  menubutton $bar.file -text "File   " -underline 0 -menu $bar.file.menu
  $bar.file config -width 5
  pack $bar.file -side left
  menu $bar.file.menu -tearoff no
  $bar.file.menu add command -label "Load colvars" -command [namespace code {plumedVis::browse "readcolvars" "Please select the COLVAR file to read in"} ]
  $bar.file.menu add command -label "Load fes" -command [namespace code {plumedVis::browse "readfes" "Please select the 2d FES to read in"} ]
  $bar.file.menu add command -label "Load onions" -command [namespace code {plumedVis::browse "readonions" "Select the ONIONS file from a recon metad simuiation"} ]
  $bar.file.menu add command -label "Print graph" -command [namespace code plumedVis::printCanvas]    
  $bar.file.menu add command -label "Save projection file" -command [namespace code plumedVis::openProjectionFileWindow] 
  $bar.file.menu add command -label "Reset molecule" -command [namespace code plumedVis::reset]   

  # Calculate menu
  menubutton $bar.calc -text "Calculate   " -underline 0 -menu $bar.calc.menu
  $bar.calc config -width 10
  pack $bar.calc -side left
  menu $bar.calc.menu -tearoff no
  $bar.calc.menu add command -label "Calculate colvars" -command [namespace code {plumedVis::browse "calccolvars" "Please select a plumed input file"} ]
  $bar.calc.menu add command -label "Calculate torsional angles" -command [namespace code {plumedVis::getColvars "torsions"} ] 
  $bar.calc.menu add command -label "Calculate contact map" -command [namespace code {plumedVis::getColvars "drmsd"} ]
  $bar.calc.menu add command -label "Calculate discretized rdf" -command [namespace code {plumedVis::getColvars "rdf"} ]
  $bar.calc.menu add command -label "Calculate discretized adf" -command [namespace code {plumedVis::getColvars "adf"} ]
  $bar.calc.menu add command -label "Save selection" -command [namespace code plumedVis::storeSelection]
  $bar.calc.menu add command -label "Extract landmark points" -command [namespace code plumedVis::initfps]
  $bar.calc.menu add command -label "Dimensionality reduction" -command [namespace code plumedVis::initdimred]
  $bar.calc.menu add command -label "Out of sample embedding" -command [namespace code plumedVis::initosample]
  $bar.calc.menu add command -label "Sum hills" -command [namespace code {plumedVis::browse "sumhills" "Please select a HILLS file"} ]
  $bar.calc.menu add command -label "Convergence of selected free energy" -command [namespace code plumedVis::drawConvergencePlot]
  $bar.calc.menu add command -label "Find reconnaissance basins" -command [namespace code {plumedVis::browse "findbasins" "Please select a BASINS file"} ]

  # Appearance menu
  menubutton $bar.look -text "Appearance   " -underline 0 -menu $bar.look.menu
  $bar.look config -width 15
  pack $bar.look -side left
  menu $bar.look.menu -tearoff no
  $bar.look.menu add command -label "Axis appearance" -command [namespace code plumedVis::openAxisLook]
  $bar.look.menu add command -label "Colvar appearance" -command [namespace code plumedVis::openCVlook]
  $bar.look.menu add command -label "Fes appearance" -command [namespace code plumedVis::openfeslook]
  $bar.look.menu add command -label "Automatic axis" -command [namespace code plumedVis::updateColvars]

  # Help menu
  menubutton $bar.help -text "Help   " -menu $bar.help.menu
  $bar.help config -width 5
  pack $bar.help -side right
  menu $bar.help.menu -tearoff no
  $bar.help.menu add command -label "G.I.S.M.O. Help..." -command "vmd_open_url http://sketchmap.berlios.de"

  return $bar
}

proc plumedVis::createControls {w} {
  variable fesdisplay
  frame $w -relief raised -bd 2
  axisMenu $w.xaxis "x-axis" plumedVis::xcoord 5
  axisMenu $w.yaxis "y-axis" plumedVis::ycoord 5
  axisMenu $w.saxis "size" plumedVis::scoord 5 
  axisMenu $w.caxis "color" plumedVis::ccoord 5
  button $w.allcvs -text "all cvs" -relief raised -command plumedVis::plotAllCVs
  axisMenu $w.faxis "fes data" plumedVis::fcoord 18 

  # This creates show basins
  frame $w.bas
  label $w.bas.lab -text "hide \n colvars" -anchor e
  checkbutton $w.bas.box -variable plumedVis::hideColvars
  grid $w.bas.lab -row 1 -column 1
  grid $w.bas.box -row 1 -column 2

  # This bit creates the fes frame spinbox
  frame $w.fram
  label $w.fram.slab -text "fes frame" -anchor e
  ::fesTools::createSpinbox $w.fram.spinb -textvariable plumedVis::fesdisplay
  label $w.fram.ldiff -text "show \n difference" -anchor e
  checkbutton $w.fram.diff -variable plumedVis::showfesdiff
  grid $w.fram.slab -row 1 -column 1
  grid $w.fram.spinb -row 1 -column 2
  grid $w.fram.ldiff -row 1 -column 3
  grid $w.fram.diff -row 1 -column 4

  # Put all the controls into the frame now 
  grid $w.xaxis  -row 1 -column 1
  grid $w.yaxis  -row 1 -column 2
  grid $w.saxis  -row 2 -column 1
  grid $w.caxis  -row 2 -column 2
  grid $w.allcvs -row 1 -column 3
  grid $w.bas -row 2 -column 3
  grid $w.faxis  -row 1 -column 4
  grid $w.fram   -row 2 -column 4

  return $w
}

proc plumedVis::openAxisLook {args} {
  variable w
  variable axiswindow
  variable xcoord
  variable ycoord
  variable fcoord

  if { $xcoord=="none" || $ycoord=="none" } {
     if { $fcoord=="none" } { return }
  }

  if [winfo exists .pvis_axis] {
     wm deiconify $axiswindow
     return      
  }

  set axiswindow [toplevel ".pvis_axis"]
  wm title $axiswindow "Axis Appearance"
  wm resizable $axiswindow 0 0 

  frame $axiswindow.top -padx 1m -pady 1m
  pack [ ::gtPlot::axis_window $w.fr.canv $axiswindow.top.fr ] -in $axiswindow.top -side top -fill both
  frame $axiswindow.d -padx 3m -pady 3m
  pack [ button $axiswindow.d.ok -text "Apply" -relief raised -command [namespace code plumedVis::drawGraph ] ] -in $axiswindow.d -side left
  pack [ button $axiswindow.d.dismiss -text "Dismiss" -relief raised -command { destroy .pvis_axis } ] -in $axiswindow.d -side right
  pack $axiswindow.d -in $axiswindow.top -side top
  pack $axiswindow.top -fill both
}

proc plumedVis::openfeslook {args} {
  variable w
  variable feswindow

  if [winfo exists .pvis_feslook] {
    wm deiconify $feswindow 
    return
  } 
  set feswindow [toplevel ".pvis_feslook"]
  wm title $feswindow "Free energy appearance"
  wm resizable $feswindow 0 0

  frame $feswindow.top -padx 1m -pady 1m
  pack [ ::fesTools::appearance $w.fr.canv $feswindow.top.fr ] -in $feswindow.top -side top -fill both
  frame $feswindow.d -padx 3m -pady 3m
  pack [ button $feswindow.d.ok -text "Apply" -relief raised -command [namespace code plumedVis::drawGraph] ] -in $feswindow.d -side left
  pack [ button $feswindow.d.dismiss -text "Dismiss" -relief raised -command { destroy .pvis_feslook} ] -in $feswindow.d -side right
  pack $feswindow.d -in $feswindow.top -side top
  pack $feswindow.top -fill both
}

proc plumedVis::openCVlook {args} {
  variable w
  variable cvwindow
  variable allcolvars
  variable xcoord
  variable ycoord

  if { $allcolvars==1 } { return }
  if { $xcoord=="none" || $ycoord=="none" } { return }

  if [winfo exists .pvis_cvlook] {
     wm deiconify $cvwindow
     return
  }
  set cvwindow [toplevel ".pvis_cvlook"]
  wm title $cvwindow "CV Plotter Appearance"
  wm resizable $cvwindow 0 0

  frame $cvwindow.top -padx 1m -pady 1m
  pack [ ::cv_traj::cvlook_window $cvwindow.top.fr ] -in $cvwindow.top -side top -fill both
  frame $cvwindow.d -padx 3m -pady 3m
  pack [ button $cvwindow.d.ok -text "Apply" -relief raised -command [namespace code plumedVis::updateColvars] ] -in $cvwindow.d -side left 
  pack [ button $cvwindow.d.dismiss -text "Dismiss" -relief raised -command { destroy .pvis_cvlook} ] -in $cvwindow.d -side right
  pack $cvwindow.d -in $cvwindow.top -side top
  pack $cvwindow.top -fill both
}

proc plumedVis::axisMenu {w textlabel var width} {
   frame $w
   label $w.lab -text $textlabel -anchor w
   menubutton $w.m -relief raised -bd 2 -direction flush -width $width \
         -textvariable $var -menu $w.m.menu
   menu $w.m.menu  

   grid $w.lab -row 1 -column 1 -sticky we
   grid $w.m -row 1 -column 2 -sticky we

   return $w
}

# This finds the current molecule that is on top
proc plumedVis::findMolecule {args} {
   # Get the new molid
   foreach m [molinfo list] {
      if { [molinfo $m get top ] } { set topmol $m ; molinfo $m set active 1 }
      # There should only ever be one molecule drawn at a time (otherwise we can't make sense of colvars)
      if { [molinfo $m get drawn] && ![molinfo $m get top] } { molinfo $m set active 0 }
   }      
   return $topmol 
}   

# This returns a tempory directory
proc plumedVis::tmpdir {args} {
   global tcl_platform
   switch $tcl_platform(platform) {
       unix {
           set tmpdir /tmp   ;# or even $::env(TMPDIR), at times.
       } macintosh {
           set tmpdir $::env(TRASH_FOLDER)  ;# a better place?
       } default {
           set tmpdir [pwd]
           catch {set tmpdir $::env(TMP)}
           catch {set tmpdir $::env(TEMP)}
       }
   }
}

proc plumedVis::getargs {arglist tag {warning 0} } {
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


