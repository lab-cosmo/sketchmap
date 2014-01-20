#
# Utilities for plotting CVs over a trajectory and linking them to 
# the currently shown frame
# 
# $Id: pcolvarTools.tcl, v 1.0 2011/04/24
# 
# (1) Add functionality to control the size and colors on input points - so we can do color according to and size according to weight
# (4) Weights must be a number between zero and one (actually I need to think about this some more)

source $env(VMDDIR)/scripts/init.d/gtplot.tcl

package provide cv_traj 1.0
package require gtPlot 2.0

namespace eval ::cv_traj {
# Routines we export
  namespace export setDefaults
  namespace export readData
  namespace export plotData
  namespace export delete
  namespace export Goto
  namespace export addLandmarkSelection 
  namespace export addSelection
  namespace export clearSelection
  namespace export getSelection
  namespace export cvlook_window
# Variables we export
  namespace export highlight_nbefore
  namespace export highlight_nafter
  namespace export cvstride
  namespace export trajcol
  namespace export boxSize
  namespace export cvcolors
  namespace export pscale
# Declaration of variables 
  variable cvstride          1  ;# stride between displayed cv values
  variable highlight_nbefore 0  ;# Number of frames to highlight before the current point
  variable highlight_nafter  0  ;# Number of frames to highlight after the current point
  variable highlight_lines   1  ;# Draw lines between data points in highlight
  variable trajcol      "grey"  ;# This is the colorscale used for drawing the trajectory (grey/jet)
  variable cvcolors     "hot"   ;# This is the colorscale used for drawing colvars 
  variable boxSize          3   ;# Default size of boxes
  variable pscale           3   ;# Default ratio between smallest and largest points
  variable framemap             ;# The map that tells us how to go from a trajectory frame to a point
  variable selection            ;# The list of selected points
}

proc cv_traj::setDefaults {args} {
  variable cvstride
  variable highlight_nbefore 
  variable highlight_nafter  
  variable highlight_lines
  variable trajcol      
  variable boxSize 
  variable selection 

  set cvstride 1
  set highlight_nbefore 0
  set highlight_nafter  0
  set highlight_lines   1
  set trajcol      "grey"
  set boxSize           3
  set selection        {}
} 

proc cv_traj::readData {args} {
  variable framemap
  variable cvcolors

  set datatag [getargs $args "-datatag" "In cv_traj readData no frames data"]
  set framemap($datatag) [getargs $args "-frames" "In cv_traj no frames data"]
  set ndata [getargs $args "-npoints" "In cv_traj no npoints data"]
  set xdata [getargs $args "-xdata" "In cv_traj no xdata"]
  set ydata [getargs $args "-ydata" "In cv_traj no ydata"]
  set cdata [getargs $args "-cdata" "In cv_traj no cdata"]
  set weights [getargs $args "-weights" "In cv_traj no weights"]

  # This does stuff for the weights
  #if { [llength $framemap($datatag)]==$ndata } {
  #    set weights [getargs $args "-weights" "In cv_traj no weights"]
  #} else {
  #    # Calculate the weight of every point from the framemap
  #    for { set i 0 } { $i<$ndata } { incr i } { set w($i) 0 }
  #    foreach index $framemap($datatag) { incr w($index) }
  #    for { set i 0 } { $i<$ndata } { incr i } { lappend weights $w($i) }
  #} 

  # Pass the data to the graph plotter
  ::gtPlot::readData -datatag $datatag: -npoints $ndata -xdata $xdata -ydata $ydata -weights $weights -colors $cdata -cscale $cvcolors
}

proc cv_traj::delete {canv args} {
  set datatag [getargs $args "-datatag" "In cv_traj delete did not specify a datatag for data"] 
  ::gtPlot::unplotData $canv -datatag $datatag:
}

proc cv_traj::plotData {canv args} {
  variable cvstride
  variable boxSize
  variable trajcol 
  variable pscale

  global vmd_frame

  # Get the datatag for the data 
  set datatag [getargs $args "-datatag" "In cv_traj plotData did not specify a datatag for data"]

  # Plot the data 
  ::gtPlot::plotData $canv -datatag $datatag: -plotstyle "scatter" -stride $cvstride -pointsize $boxSize -pointscaling $pscale
  updateHighlight $canv -datatag $datatag   ; # And sort out the highlighting  
}

proc cv_traj::Goto {w args} {
  variable framemap
  set datatag [getargs $args "-datatag" "In cv_traj Goto did not specify a datatag for data"]

  set id [ $w find withtag current]
  set taglist [$w gettags $id]
  set listindex [lsearch -glob $taglist $datatag:*]
  if { $listindex < 0 } { return }
  set frametag [lindex $taglist $listindex]
  lassign [split $frametag :] foo point

  # Now use the framemap to get the frame
  set frame 0
  foreach fp $framemap($datatag) {
    if { $fp==$point } { break }
    incr frame 
  } 
  animate goto $frame
}

proc cv_traj::updateHighlight {canv args} {
  variable highlight_nbefore
  variable highlight_nafter
  variable highlight_lines 
  variable trajcol
  variable framemap

  global vmd_frame

  set datatag [getargs $args "-datatag" "In cv_traj updateHighlight not specify a datatag for data"]

  # Check for non insane values of highlight_nbefore and highlight_nafter
  if { ![string is integer -strict $highlight_nbefore] || $highlight_nbefore<0 } { return }
  if { ![string is integer -strict $highlight_nafter] || $highlight_nafter<0 } { return }

  # Find the top molecule
  foreach m [molinfo list] {
     if { [molinfo $m get top] } { set topmol $m }
  }

  # We now find the frame we are currently on and map it to the points
  set center [ lindex $framemap($datatag) $vmd_frame($topmol)]
  set beforelist {}       ;    set afterlist {}
  # Get the frames after
  for { set i 1 } { $i<=$highlight_nafter } { incr i } { 
      set ind [ expr $vmd_frame($topmol) + $i ]
      if { $ind<[molinfo top get numframes] } { 
         set tmp [lindex $framemap($datatag) $ind]  
         if { $tmp!=$center } {
            set flag 0
            for { set j 0 } { $j<[llength $afterlist] } { incr j } {
                if { [lindex $afterlist $j]==$tmp } { set flag 1 } 
            }
            if { $flag==0 } { lappend afterlist $tmp }
         } 
      }
  }
  
  # Get the frames before
  for { set i 1 } { $i<=$highlight_nbefore } { incr i } {
      set ind [expr $vmd_frame($topmol) - $i]
      if { $ind>=0 } {
         set tmp [lindex $framemap($datatag) $ind]
         if { $tmp!=$center } {
            set flag 0
            for { set j 0 } { $j<[llength $afterlist] } { incr j } {
                if { [lindex $afterlist $j]==$tmp } { set flag 1 }
            }
            for { set j 0 } { $j<[llength $beforelist] } { incr j } {
                if { [lindex $beforelist $j]==$tmp } { set flag 1 }
            } 
            if { $flag==0 } { lappend beforelist $tmp }
         }  
      }
  }  

  # This gets rid of any selections we might have made
  clearSelection $canv -datatag $datatag
#  ::gtPlot::highlight $canv -type "onecolor" -datatag $datatag: -listp {} -color "red"
  # This actually draws the highlight
  ::gtPlot::highlight $canv -type "scaled" -datatag $datatag: -center $center -lines $highlight_lines -before $beforelist -after $afterlist -colors $trajcol
}

proc cv_traj::addSelection {canv args} {
   variable selection
   set datatag [getargs $args "-datatag" "In cv_traj addSelection did not specify a datatag for data"]
   # Find the contents of the lasso
   set tmp [ ::gtPlot::getLasso -datatag $datatag: -data indices ]
   # Add the contents to the selection 
   for { set i 0 } { $i<[llength $tmp] } { incr i } { lappend selection [lindex $tmp $i] }
   # And draw the selection on the graph
   createSelection $canv -datatag $datatag -color "red" -selection $selection
}

proc cv_traj::createSelection {canv args} {
  variable selection
  set datatag [getargs $args "-datatag" "In cv_traj createSelection did not specify a datatag for data"]
  set color [getargs $args "-color" "In cv_traj createSelection did not specify a color for selection"]
  set selection [getargs $args "-selection" "In cv_traj createSelection did not specify a selection"]
  # Now highlight the selection
   ::gtPlot::highlight $canv -type "onecolor" -datatag $datatag: -listp $selection -color $color  
}

proc cv_traj::clearSelection {canv args} {
  variable selection
  set datatag [getargs $args "-datatag" "In cv_traj clearSelection did not specify a datatag for data"] 
  set selection {}
  ::gtPlot::highlight $canv -type "onecolor" -datatag $datatag: -listp $selection -color "red"
}

proc cv_traj::getSelection {args} {
  variable selection
  variable framemap

  if { [llength $selection]==0 } { 
     tk_messageBox -icon error -type ok -title Message -message "No selection has been made"
     return 
  }
   
  set datatag [getargs $args "-datatag" "In cv_traj getSelection did not specify a datatag for data"]
  set framelist {}  ;  set frame 0
  foreach p $selection {
     set frame 0
     foreach fp $framemap($datatag) {
        if { $fp==$p } {
           lappend framelist $frame
           break
        }
        incr frame
     }
  }
  return $framelist
}

proc cv_traj::cvlook_window {wind args} {

  frame $wind -padx 1m -pady 1m

  labelframe $wind.afram -relief ridge -bd 2 -text "plot appearance" -padx 2m -pady 2m
  label $wind.afram.lab -text "Plot stride" -anchor e
  entry $wind.afram.ent -width 5 -textvariable cv_traj::cvstride
  label $wind.afram.slab -text "Point size" -anchor e
  eval spinbox $wind.afram.sspin -textvariable cv_traj::boxSize -from 1 -to 10 -increment 1 -width 5
  label $wind.afram.clab -text "Point scaling" -anchor e
  eval spinbox $wind.afram.cspin -textvariable cv_traj::pscale -from 1 -to 10 -increment 1 -width 5
  label $wind.afram.col -text "CV color scale" -anchor e
  ::gtPlot::colorMenu $wind.afram.cscal -textvariable cv_traj::cvcolors

  grid $wind.afram.lab -row 1 -column 1 -sticky e
  grid $wind.afram.ent -row 1 -column 2 -sticky e
  grid $wind.afram.slab -row 2 -column 1 -sticky e
  grid $wind.afram.sspin -row 2 -column 2 -sticky e
  grid $wind.afram.clab -row 3 -column 1 -sticky e
  grid $wind.afram.cspin -row 3 -column 2 -sticky e
  grid $wind.afram.col -row 4 -column 1 -sticky e
  grid $wind.afram.cscal -row 4 -column 2 -sticky e
  pack $wind.afram -in $wind -side top -fill both

  labelframe $wind.hfram -relief ridge -bd 2 -text "highlight appearance" -padx 2m -pady 2m
  label $wind.hfram.lnb4 -text "Highlight frames before current"
  entry $wind.hfram.enb4 -width 5 -textvariable cv_traj::highlight_nbefore
  label $wind.hfram.lna4 -text "Highlight frames after current"
  entry $wind.hfram.ena4 -width 5 -textvariable cv_traj::highlight_nafter
  label $wind.hfram.clab -text "Colorscale for highlight"
  ::gtPlot::colorMenu $wind.hfram.c -textvariable cv_traj::trajcol 
  label $wind.hfram.llla -text "Draw line through highlighted points"
  checkbutton $wind.hfram.lbox -variable cv_traj::highlight_lines

  grid $wind.hfram.lnb4 -row 1 -column 1 -sticky e
  grid $wind.hfram.enb4 -row 1 -column 2 -sticky w
  grid $wind.hfram.lna4 -row 2 -column 1 -sticky e
  grid $wind.hfram.ena4 -row 2 -column 2 -sticky w
  grid $wind.hfram.clab -row 3 -column 1 -sticky e
  grid $wind.hfram.c -row 3 -column 2 -sticky w
  grid $wind.hfram.llla -row 4 -column 1 -sticky e
  grid $wind.hfram.lbox -row 4 -column 2 -sticky w
  pack $wind.hfram -in $wind -side top -fill both

  return $wind
}

# This gets arguments or returns a warning   
proc cv_traj::getargs {arglist tag {warning 0} } {
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

