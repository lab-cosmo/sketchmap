# GT plot v 2.0
# 
# Gareth Tribello
#
# To do list
#
# (1) If you are feeling keen - if you have two windows: if you plot data in window 1 and then attempt to zoom on window 2 this won't work (all the data will disapear)

package provide gtPlot 2.0

namespace eval ::gtPlot {
# Procedures to export
    namespace export gtcreate                ; # Creates a gtPlot object 
    namespace export readData                ; # Reads in a data set
    namespace export plotData                ; # Plot the data on a graph
    namespace export unplotData              ; # Removes plotted data
    namespace export highlight               ; # Highlights a set of data points
    namespace export createScalez            ; # Draw a colorscale for 2d plots
    namespace export axis_window             ; # Opens the axis properties dialogue window 
    namespace export zscale_window           ; # Controls the appearance of the zscale
    namespace export printCanvas             ; # Print out the current canvas
    namespace export colorMenu               ; # Creates a menu that allows one to control plots 
    namespace export getLasso                ; # Passes the list of points in a selection to an external routine
    namespace export drawShape               ; # Draw a shape on the canvas
# Variables to export
    namespace export bounds
    namespace export doneLasso
    variable ncolors       100  ; # Number of colors to use to plot free energy surfaces and colored colvars  
    variable maxgrid       100  ; # Maximum number of points in a grid that you allow to plot (depends on how long you are willing to wait)
    variable ntics_d         3  ; # Default number of major tics
    variable paddding           ; # Controls the width of the padding around the graph
    variable font               ; # What font do we use for writing numbers
    variable bounds             ; # The automatic boundaries for the graph
    variable zoomID             ; # The device used by zoom
    variable lassoID            ; # Stuff used by lasso
    variable doneLasso          ; # Pragmatic solution to not working trace execute LassoEnd
}

proc gtPlot::setDefaultAxes {canv args} {
   variable bounds
   set bounds([list $canv xfixtics])       0
   set bounds([list $canv yfixtics])       0
   set bounds([list $canv xmajticl])       4.0
   set bounds([list $canv ymajticl])       4.0
   set bounds([list $canv xminticl])       2.0
   set bounds([list $canv yminticl])       2.0
   set bounds([list $canv xnmintic])       1
   set bounds([list $canv ynmintic])       1
   set bounds([list $canv txmin])          0
   # Stuff for z
   set bounds([list $canv zticspacing])    0
   set bounds([list $canv zmajticl])       4.0
   set bounds([list $canv zminticl])       2.0
   set bounds([list $canv znmintic])       1
   set bounds([list $canv zmin])           0
   set bounds([list $canv zmax])           0
}

# This creates the graph
proc gtPlot::gtcreate {canv args} {
    variable font
    variable padding

    # Interpret the arguments
    set size [getargs $args "-size" "-size flag missing in call to gtcreate"]
    set font [getargs $args "-font" "-font flag missing in call to gtcreate"]
    set padding [getargs $args "-padding" "-padding flag missing in call to gtcreate"]

    setDefaultAxes $canv    ; # Set everything for axes to its default setting

    # Creates a canvas with a border
    eval { canvas $canv -highlightthickness 0 -borderwidth 0 -background white } \
      -width [ expr $size + 2 * $padding ] -height [ expr $size + 2 * $padding ]

    # This will set up our zoom tool
    bind $canv <Shift-Button-1> [namespace code { gtPlot::ZoomStart %W %x %y} ]           ; # Set the initial point
    bind $canv <Shift-B1-Motion> [namespace code { gtPlot::ZoomMove %W %x %y} ]           ; # Adjust the box size 
    bind $canv <Shift-B1-ButtonRelease> [namespace code { gtPlot::ZoomEnd %W %x %y} ]     ; # Delete the zoom box and zoom
    # This is for the lasso points tool
    bind $canv <Control-Button-1> [namespace code { gtPlot::LassoStart %W %x %y} ]        ; # Set the initial point
    bind $canv <Control-B1-Motion> [namespace code { gtPlot::LassoMove %W %x %y} ]        ; # Adjust the box size
    bind $canv <Control-B1-ButtonRelease> [namespace code { gtPlot::LassoEnd %W %x %y} ]  ; # Delete the box and grab points inside
    # Other tools
    bind $canv <Motion> [namespace code { gtPlot::writeCoord %W %x %y} ]                  ; # This is for the little thing that keeps track of our position
    bind $canv <Configure> [namespace code { gtPlot::drawGraph %W } ]                     ; # This redraws the graph when the window is resized  

    return
}

proc gtPlot::readData {args} {
    variable data
    variable bounds
    variable ncolors

    # Get the number of data points
    set tag [getargs $args "-datatag" "-datatag flag not specified in call to gtPlot::readData"]
    set data([list $tag np]) [getargs $args "-npoints" "-npoints flag not specified in call to gtPlot::readData"]

    # Store the tag names 
    set taglist {}
    foreach keyitem [array names data] {
       # The != ensures we overwrite any data with the same tag
       if { [ lsearch $keyitem tag ] != -1 && $data($keyitem) != $tag } { lappend taglist $data($keyitem) }
    }
    set nsets [llength $taglist]
    for { set i 0 } { $i<$nsets } { incr i } { set data([list tag $i]) [lindex $taglist $i] }
    set data([list tag $nsets]) $tag

    # Get x data
    set data([list $tag x]) [getargs $args "-xdata" "-xdata flag not specified in call to gtPlot::readData"]
    getRange $tag x
     
    # Get y data
    set data([list $tag y]) [getargs $args "-ydata" "-ydata flag not specified in call to gtPlot::readData"]
    getRange $tag y

    # Get weights data
    set warray [getargs $args "-weights" 0]
    if { $warray==0 } { 
       set data([list $tag w]) {}
       for { set i 0 } { $i<$data([list $tag np]) } { incr i } { lappend data([list $tag w]) 1.0 }
    } else {
       # Find range in weights
       set wmax [lindex $warray 0]
       set wmin [lindex $warray 0]
       for {set i 1} { $i<$data([list $tag np]) } { incr i } {
           set wtmp [lindex $warray $i]
           if { $wtmp>$wmax } { set wmax $wtmp }
           if { $wtmp<$wmin } { set wmin $wtmp }
       }
       set wrange [expr 1.0*($wmax - $wmin)] 
       # Now create the weights
       if { $wrange==0 } {
          set data([list $tag w]) {}
          for {set i 0} { $i<$data([list $tag np]) } { incr i } { lappend data([list $tag w]) 0 }
       } else {
          set data([list $tag w]) {}
          for {set i 0} { $i<$data([list $tag np]) } { incr i } {
              lappend data([list $tag w]) [ expr 1.0*( [lindex $warray $i] - $wmin ) / double($wrange) ]
          }
       } 
    } 

    # Get z data (if its there)   -- Must put this back in once we have 3d data
    set data([list $tag z]) [getargs $args "-zdata" 0] 
    if { $data([list $tag z])==0 } {
        set data([list $tag is3d]) 0
    } else { 
        set data([list $tag is3d]) 1
        # Get zrange if we are not using one from previous calculation
        getRange $tag z    
        # set bounds([list $canv zmin]) 0  ; set bounds([list $canv zmax]) 0   # This should be done but it is boring and hard GAT  
    } 

    # Get color data
    set incdata [getargs $args "-colors" 0]
    if { [llength $incdata]==1 } {
       set data([list $tag color]) {}
       for { set i 0 } { $i<$data([list $tag np]) } { incr i } { lappend data([list $tag color]) [lindex $incdata 0] }
    } elseif { [llength $incdata]==$data([list $tag np]) } {
       set data([list $tag color]) {}
       set cmin [lindex $incdata 0]
       set cmax [lindex $incdata 0]
       for {set i 1} { $i<$data([list $tag np]) } { incr i } {
           set ctmp [lindex $incdata $i]
           if { $ctmp>$cmax } { set cmax $ctmp }
           if { $ctmp<$cmin } { set cmin $ctmp }
       }
       # Create a colorscale
       set cmapt [getargs $args "-cscale" "reading data with color scale but scale was not specified use -scale"]
       set color_incr [ expr ( $cmax - $cmin )  / double($ncolors) ]
       set colormap [ createColormap $cmapt $ncolors] 
       for { set i 0 } { $i<$data([list $tag np]) } { incr i } { 
           set kk [ expr int( ( [lindex $incdata $i] - $cmin ) / $color_incr ) ]
           lappend data([list $tag color]) [ lindex $colormap $kk ]
       }
    } elseif { $data([list $tag is3d])!=1 } {
       tk_messageBox -icon error -type ok -title Message -message "Missing -colors flag in call to read data"
    }

    # Currently there are no highlights
    set data([list $tag arehighlights]) 0
    # We want to see every point
    set data([ list $tag stride ]) 1

    if { [llength $data([list $tag color])]!=$data([list $tag np]) || [llength $data([list $tag w])]!=$data([list $tag np]) } {
        tk_messageBox -icon error -type ok -title Message -message "Either length or color array is wrong size [llength $data([list $tag color])] [llength $data([list $tag w])]"
    }

}

proc gtPlot::plotData {canv args} {
    variable data
    variable bounds

    # Get compulsory arguments
    set tag [getargs $args "-datatag" "-datatag flag not specified in call to gtPlot::plotData"]

    # GAT this can probably go in the future
    # set pos [lsearch $args "-color"]
    # if { $pos<0 } { 
    #    set color "white" 
    # } else {
    #    set color [ lindex $args [expr $pos+1 ] ]
    # }

    # Suss out what we need to do based on the type
    set data([list $tag drawpoints]) 0 ; set data([list $tag drawline]) 0 ; set data([list $tag drawbars]) 0
    set type [getargs $args "-plotstyle" 0]
    if { $type=="scatter" } {
        set data([list $tag drawpoints]) 1
        set data([ list $tag psize]) [getargs $args "-pointsize" "-pointscale is not specified in call to gtPlot::plotData"]
        set data([ list $tag pscale]) [getargs $args "-pointscaling" "-pointscaling is not specified in call to gtPlot::plotData"]
        set data([list $tag stride]) [getargs $args "-stride" 0]
        if { $data([list $tag stride])==0 } { set data([list $tag stride]) 1 }
    } elseif { $type=="surface" } {
        if { $data([ list $tag is3d ])!=1 } {
           tk_messageBox -icon error -type ok -title Message -message "Cannot plot a surface with 2d data"
           return 
        } 
        set data([ list $tag fescolors]) [getargs $args "-surfacecolors" "-surfacecolors is not specified in call to gtPlot::plotData"]
        #set data([ list $tag ncolors]) [getargs $args "-ncolors" "-ncolors is not specified in call to gtPlot::plotData"]
        setupGrid $tag   ; # Put the data onto a 2d grid
    } elseif { $type=="bar" } {
        set data([list $tag drawbars]) 1
        set data([list $tag color]) [getargs $args "-linecolor" "-linecolor flag is not specified in call to gtPlot::plotData"]
        # For bar charts we have to extend the x axis from the range in the data
        set data([ list $tag x max]) [expr $data([ list $tag x max ]) + [lindex $data([list $tag x]) 1] - [lindex $data([list $tag x]) 0] ]
    } elseif { $type=="line" } {
        set data([list $tag drawline]) 1
        set data([list $tag color]) [getargs $args "-linecolor" "-linecolor flag is not specified in call to gtPlot::plotData"]
    } else {
        tk_messageBox -icon error -type ok -title Message -message "Unsuported graph type in call to gtPlot::plotData"
    }    

    # Delete any old plot with this tag
    $canv delete $tag

    # This keeps track of which data is plotted in this window
    set data([ list $tag $canv displayed ]) 1

    # If user has not set axis using dialogue then get them automatically
    if { $bounds([list $canv txmin])==0 } {
       set bounds([list $canv txmin]) [getAutoBound $canv "min" "x"]  ;   set bounds([list $canv txmax]) [getAutoBound $canv "max" "x"]
       set bounds([list $canv tymin]) [getAutoBound $canv "min" "y"]  ;   set bounds([list $canv tymax]) [getAutoBound $canv "max" "y"] 
    }

    # Set the axes 
    setAxis $canv $bounds([list $canv txmin]) $bounds([list $canv txmax]) $bounds([list $canv tymin]) $bounds([list $canv tymax]) 

    # Draw the graph
    drawGraph $canv
}

proc gtPlot::drawGraph {canv} {
    variable data

    set flag 0 
    foreach keyitem [array names data] {
       # Look for all the tag names
       if { [ lsearch $keyitem tag ] != -1 } {
          set tag $data($keyitem)
          # Look to see if this tag is displayed in this canvas
          if { [info exists data([ list $tag $canv displayed ])] } { set flag 1 }  
       }
    }
    if { $flag==0 } { return }

    # Delete old axes and old data points
    $canv delete all

    # Draw the new axis here
    drawAxis $canv "x"
    drawAxis $canv "y"

    set flag 0
    foreach keyitem [array names data] {
       if { [lsearch $keyitem tag ] != -1 } {
          set tag $data($keyitem)
          if { [info exists data([ list $tag $canv displayed ])] } {
               if { $data([ list $tag drawbars ])==1 } {
                   barplotter $canv $tag   ; # This plots a bar chart
               } elseif { $data([ list $tag is3d ]) == 0 } {
                   2dplotter $canv $tag    ; # This plots data points
               } elseif { $data([ list $tag is3d ]) == 1 && $flag == 0 } {
                   3dplotter $canv $tag    ; # This plots free energy surfaces  
                   set flag 1
               } else {
                   tk_messageBox -icon error -type ok -title Message -message "Can't plot more than one fes at once"
                   return
               }
          }
       }
    }
    # Put the fes at the bottom of the stack - this may yet go back GAT
    if { $flag>0 } { $canv lower surf displayed ; $canv lower surf axis }
}

# This removes data from the plot
proc gtPlot::unplotData {canv args} {
    variable data
    variable bounds

    set tag [getargs $args "-datatag" "-datatag missing in call to gtplot::unplotData"]     

    # Check there is data plotted with this tag
    if { [llength [$canv find withtag $tag] ]!=0 } {
       # Delete the data with this tag
       $canv delete $tag
       $canv delete basinsize   ; # GAT check this

       # Keep track of what data is displayed
       unset data([ list $tag $canv displayed ])

       # If this is a fes delete the zticspacing
       if { $data([list $tag is3d])==1 } { set bounds([list $canv zticspacing]) 0 } 
    }
    # Reset the axis
    set bounds([list $canv txmin]) 0
    # set bounds([list $canv zmin]) 0 
    # set bounds([list $canv zmax]) 0

    # Check if any data is left (if there is none delete the axis)
    set flag 0
    foreach keyitem [array names data] {
       # Look for all the tag names
       if { [ lsearch $keyitem tag ] != -1 } {
            set tag $data($keyitem)
            # Check if data is displayed
            if { [info exists data([ list $tag $canv displayed ])] } { set flag 1 }
       }
    }
    if { $flag == 0 } {
       $canv delete axis
       $canv delete whereIam
    }
}

proc gtPlot::drawShape {canv args} {
   set type [getargs $args "-shape" "No shape given in call to drawShape"]
   if { $type=="rectangle"} {
      set corners [getargs $args "-corners" "No corners given in call to drawShape/rectangle"]
      set lstyle [getargs $args "-linestyle" "No linestyle given in call to drawShape/rectangle"]
      set tag [getargs $args "-shapetag" "No tag given for shape"]

      drawRectangle $canv $corners $lstyle $tag
   } else {
      puts "No procedure to draw shape $type"
   } 
}

proc gtPlot::drawRectangle {canv corners lstyle tag} {
   variable bounds

   set x1 [ point2canvas $canv "x" [lindex $corners 0] ]
   set y1 [ point2canvas $canv "y" [lindex $corners 1] ] 
   set x2 [ point2canvas $canv "x" [lindex $corners 2] ]
   set y2 [ point2canvas $canv "y" [lindex $corners 3] ]

   if { $lstyle=="solid" } {
      $canv create rectangle $x1 $y1 $x2 $y2 -tag $tag 
   } elseif { $lstyle=="dashed" } {
      $canv create rectangle $x1 $y1 $x2 $y2 -tag $tag -dash { 4 4 }
   } else {
      puts "Unsuported linestyle in call to drawRectangle"
   }
}

proc gtPlot::axis_window {canv window} {
    variable bounds

    set bounds([list $canv txmin]) $bounds([list $canv xmin])   ; set bounds([list $canv txmax]) $bounds([list $canv xmax])
    set bounds([list $canv tymin]) $bounds([list $canv ymin])   ; set bounds([list $canv tymax]) $bounds([list $canv ymax])

    frame $window -padx 1m -pady 1m

    labelframe $window.x -relief ridge -bd 2 -text "x axis" -padx 2m -pady 2m
    label $window.x.minl -text "Min" -anchor e
    entry $window.x.min -width 5 -textvariable gtPlot::bounds([list $canv txmin])
    label $window.x.maxl -text "Max" -anchor e
    entry $window.x.max -width 5 -textvariable gtPlot::bounds([list $canv txmax])
    grid $window.x.minl -row 1 -column 1 -sticky e
    grid $window.x.min -row 1 -column 2 -sticky e
    grid $window.x.maxl -row 1 -column 3 -sticky e
    grid $window.x.max -row 1 -column 4 -sticky e

    labelframe $window.x.majfram -relief ridge -bd 1 -text "Major Tics" -padx 1m -pady 1m
    label $window.x.majfram.spal -text "spacing" -anchor e
    entry $window.x.majfram.spa -width 5 -textvariable gtPlot::bounds([list $canv xticspacing])
    label $window.x.majfram.lenl -text "tic length" -anchor e
    spinbox $window.x.majfram.len -textvariable gtPlot::bounds([list $canv xmajticl]) -from 1.0 -to 8.0 -increment 0.25 -width 5
    label $window.x.majfram.fixtl -text "Fix tic spacing" -anchor e
    checkbutton $window.x.majfram.fixt -variable gtPlot::bounds([list $canv xfixtics])
    grid $window.x.majfram.spal -row 1 -column 1 -sticky e
    grid $window.x.majfram.spa  -row 1 -column 2 -sticky e
    grid $window.x.majfram.lenl -row 1 -column 3 -sticky e
    grid $window.x.majfram.len  -row 1 -column 4 -sticky e
    grid $window.x.majfram.fixtl -row 2 -column 2 -columnspan 2 -sticky e
    grid $window.x.majfram.fixt -row 2 -column 4 -sticky e

    labelframe $window.x.minfram -relief ridge -bd 1 -text "Minor Tics" -padx 1m -pady 1m
    label $window.x.minfram.nl -text "n. tics" -anchor e
    spinbox $window.x.minfram.n -textvariable gtPlot::bounds([list $canv xnmintic]) -from 0 -to 5 -increment 1 -width 5
    label $window.x.minfram.lenl -text "tic length" -anchor e
    spinbox $window.x.minfram.len -textvariable gtPlot::bounds([list $canv xminticl]) -from 1.0 -to 4.0 -increment 0.125 -width 5
    grid $window.x.minfram.nl -row 1 -column 1 -sticky e
    grid $window.x.minfram.n -row 1 -column 2 -sticky e
    grid $window.x.minfram.lenl -row 1 -column 3 -sticky e
    grid $window.x.minfram.len -row 1 -column 4 -sticky e

    grid $window.x.majfram -row 2 -column 1 -columnspan 4 -sticky e
    grid $window.x.minfram -row 3 -column 1 -columnspan 4 -sticky e
    pack $window.x -in $window -side left -fill both

    # Add y axis controls
    labelframe $window.y -relief ridge -bd 2 -text "y axis" -padx 2m -pady 2m
    label $window.y.minl -text "Min" -anchor e
    entry $window.y.min -width 5 -textvariable gtPlot::bounds([list $canv tymin])
    label $window.y.maxl -text "Max" -anchor e
    entry $window.y.max -width 5 -textvariable gtPlot::bounds([list $canv tymax])
    grid $window.y.minl -row 1 -column 1 -sticky e
    grid $window.y.min -row 1 -column 2 -sticky e
    grid $window.y.maxl -row 1 -column 3 -sticky e
    grid $window.y.max -row 1 -column 4 -sticky e

    labelframe $window.y.majfram -relief ridge -bd 1 -text "Major Tics" -padx 1m -pady 1m
    label $window.y.majfram.spal -text "spacing" -anchor e
    entry $window.y.majfram.spa -width 5 -textvariable gtPlot::bounds([list $canv yticspacing])
    label $window.y.majfram.lenl -text "tic length" -anchor e
    spinbox $window.y.majfram.len -textvariable gtPlot::bounds([list $canv ymajticl]) -from 1.0 -to 8.0 -increment 0.25 -width 5
    label $window.y.majfram.fixtl -text "Fix tic spacing" -anchor e
    checkbutton $window.y.majfram.fixt -variable gtPlot::bounds([list $canv yfixtics])
    grid $window.y.majfram.spal -row 1 -column 1 -sticky e
    grid $window.y.majfram.spa  -row 1 -column 2 -sticky e
    grid $window.y.majfram.lenl -row 1 -column 3 -sticky e
    grid $window.y.majfram.len  -row 1 -column 4 -sticky e
    grid $window.y.majfram.fixtl -row 2 -column 2 -columnspan 2 -sticky e
    grid $window.y.majfram.fixt -row 2 -column 4 -sticky e

    labelframe $window.y.minfram -relief ridge -bd 1 -text "Minor Tics" -padx 1m -pady 1m
    label $window.y.minfram.nl -text "n. tics" -anchor e
    spinbox $window.y.minfram.n -textvariable gtPlot::bounds([list $canv ynmintic]) -from 0 -to 5 -increment 1 -width 5
    label $window.y.minfram.lenl -text "tic length" -anchor e
    spinbox $window.y.minfram.len -textvariable gtPlot::bounds([list $canv yminticl]) -from 1.0 -to 4.0 -increment 0.125 -width 5
    grid $window.y.minfram.nl -row 1 -column 1 -sticky e
    grid $window.y.minfram.n -row 1 -column 2 -sticky e
    grid $window.y.minfram.lenl -row 1 -column 3 -sticky e
    grid $window.y.minfram.len -row 1 -column 4 -sticky e

    grid $window.y.majfram -row 2 -column 1 -columnspan 4 -sticky e
    grid $window.y.minfram -row 3 -column 1 -columnspan 4 -sticky e
    pack $window.y -in $window -side right -fill both

    return $window
}

proc gtPlot::drawAxis {canv axis} {
    variable padding
    variable bounds

    if { $axis == "x" } {
        set csize [expr [ winfo width $canv ] - 2.0*$padding ]
        $canv create line $padding $bounds([list $canv yorigin_canvas]) [ expr $csize + $padding ] $bounds([list $canv yorigin_canvas]) -tag "axis"
    } else {
        set csize [expr [ winfo height $canv ] - 2.0*$padding ]
        $canv create line $bounds([list $canv xorigin_canvas]) $padding $bounds([list $canv xorigin_canvas]) [ expr $csize + $padding ] -tag "axis"
    }
    drawTics $canv $axis
}

proc gtPlot::drawTics {canv axis} {
    variable data
    variable bounds
    variable padding
    variable font

    # Check that data is plotted
    set flag 0
    foreach keyitem [array names data] {
       # Look for all the tag names
       if { [ lsearch $keyitem tag ] != -1 } {
            set tag $data($keyitem)
            # Check if data is displayed
            if { [info exists data([ list $tag $canv displayed ])]  } { set flag 1 }
       }
    }
    if { $flag == 0 } { return }

    $canv delete tics$axis

    if { $axis == "x" } {
       set csize [expr [ winfo width $canv ] - 2.0*$padding ]
       set min $bounds([list $canv xmin])
       set max $bounds([list $canv xmax])
       set ticdir 1
       set zero $bounds([list $canv xorigin_points])
       set oaxis $bounds([list $canv yorigin_canvas])
       set majticl $bounds([list $canv xmajticl])
       set minticl $bounds([list $canv xminticl])
       set nmtics $bounds([list $canv xnmintic])
       set ticspacing $bounds([list $canv xticspacing])
    } else {
       set csize [expr [ winfo height $canv ] - 2.0*$padding ]
       set min $bounds([list $canv ymin])
       set max $bounds([list $canv ymax])
       set ticdir -1
       set zero $bounds([list $canv yorigin_points])
       set oaxis $bounds([list $canv xorigin_canvas])
       set majticl $bounds([list $canv ymajticl])
       set minticl $bounds([list $canv yminticl])
       set nmtics $bounds([list $canv ynmintic])
       set ticspacing $bounds([list $canv yticspacing])
    }

    # Check that major ticspacing is a number
    if { ![string is double -strict $ticspacing] } { return }

    # Draw minor tics between origin and first major tic
    for { set j 1 } { $j <=$nmtics } { incr j } {
       set mlab [ expr $zero + ( $j * $ticspacing/($nmtics+1) ) ]
       set mp [ point2canvas $canv $axis $mlab ]
       if { $mlab<$max && $axis=="x" } {
          $canv create line $mp $oaxis $mp [ expr $oaxis + $minticl ] -tag "axis"
       } elseif { $mlab<$max && $axis=="y" } {
          $canv create line $oaxis $mp [ expr $oaxis - $minticl ] $mp -tag "axis"
       }
    }

    # Draw tics on +ve axis
    set ntics [ expr int( ( $max - $zero ) / $ticspacing) ]
    for { set i 1 } { $i <= $ntics } { incr i } {
       set lab [ format %.4g [ expr $zero + ( $i * $ticspacing ) ] ]
       set p [ point2canvas $canv $axis $lab ]
       if { $lab<$max && $axis=="x" } {
          $canv create line $p $oaxis $p [ expr $oaxis + $majticl ] -tag "axis"
          $canv create text $p [ expr $oaxis + 2 * $majticl ] -font $font -anchor n -text $lab -tag "axis"
       } elseif { $lab<$max && $axis=="y" } {
          $canv create line $oaxis $p [ expr $oaxis - $majticl ] $p -tag "axis"
          $canv create text [ expr $oaxis - 2 * $majticl ] $p -font $font -anchor e -text $lab -tag "axis"
       }
       for { set j 1 } { $j <=$nmtics } { incr j } {
          set mlab [ expr $lab + ( $j * $ticspacing/($nmtics+1) ) ]
          set mp [ point2canvas $canv $axis $mlab ]
          if { $mlab<$max && $axis=="x" } {
             $canv create line $mp $oaxis $mp [ expr $oaxis + $minticl ] -tag "axis"
          } elseif { $mlab<$max && $axis=="y" } {
             $canv create line $oaxis $mp [ expr $oaxis - $minticl ] $mp -tag "axis"
          }
       }
    }

    # Draw minor tics between origin and first major tic
    for { set j 1 } { $j <=$nmtics } { incr j } {
       set mlab [ expr $zero - ( $j * $ticspacing/($nmtics+1) ) ]
       set mp [ point2canvas $canv $axis $mlab ]
       if { $mlab<$max && $axis=="x" } {
          $canv create line $mp $oaxis $mp [ expr $oaxis + $minticl ] -tag "axis"
       } elseif { $mlab<$max && $axis=="y" } {
          $canv create line $oaxis $mp [ expr $oaxis - $minticl ] $mp -tag "axis"
       } 
    }

    # Draw tics on -ve axis
    set ntics [expr int( ( $zero - $min ) / $ticspacing) ]
    for { set i 1 } { $i <= $ntics } { incr i } {
       set lab [ format %.4g [ expr $zero - ( $i * $ticspacing ) ] ]
       set p [ point2canvas $canv $axis $lab ]
       if { $lab>$min && $axis=="x" } {
          $canv create line $p $oaxis $p [ expr $oaxis + $majticl ] -tag "axis"
          $canv create text $p [ expr $oaxis + 2 * $majticl ] -font $font -anchor n -text $lab -tag "axis"
       } elseif { $lab>$min && $axis=="y" } {
          $canv create line $oaxis $p [ expr $oaxis - $majticl ] $p -tag "axis"
          $canv create text [ expr $oaxis - 2 * $majticl ] $p -font $font -anchor e -text $lab -tag "axis"
       }
       for { set j 1 } { $j <=$nmtics } { incr j } {
          set mlab [ expr $lab - ( $j * $ticspacing/($nmtics+1) ) ]
          set mp [ point2canvas $canv $axis $mlab ]
          if { $lab>$min && $axis=="x" } {
             $canv create line $mp $oaxis $mp [ expr $oaxis + $minticl ] -tag "axis"
          } elseif { $lab>$min && $axis=="y" } {
             $canv create line $oaxis $mp [ expr $oaxis - $minticl ] $mp -tag "axis"
          }
       }
    }
}

proc gtPlot::setAxis {canv xmin xmax ymin ymax} {
    variable bounds
    variable ntics_d

    # Check the axis are sane
    if { $xmax<$xmin || $ymax<$ymin } { 
       set bounds([list $canv txmin]) $bounds([list $canv xmin])   ;   set bounds([list $canv txmax]) $bounds([list $canv xmax])
       set bounds([list $canv tymin]) $bounds([list $canv ymin])   ;   set bounds([list $canv tymax]) $bounds([list $canv ymax])
       return 
    }

    # Now set the axis for real
    set bounds([list $canv xmin]) $xmin ; set bounds([list $canv xmax]) $xmax
    set bounds([list $canv ymin]) $ymin ; set bounds([list $canv ymax]) $ymax

    # Setup other stuff for drawing points from range
    set bounds([list $canv xcenter]) [expr ( $xmin + $xmax ) / 2.0 ]
    set bounds([list $canv ycenter]) [expr ( $ymin + $ymax ) / 2.0 ]
    set bounds([list $canv xrange]) [expr $xmax - $xmin ]
    set bounds([list $canv yrange]) [expr $ymax - $ymin ]

    # Find the origin (on the canvas)
    findOrigin $canv "x"    
    findOrigin $canv "y"

    if { $bounds([list $canv xfixtics])==0 } {
       if { [expr $xmax - $bounds([list $canv xorigin_points]) ] > [ expr $bounds([list $canv xorigin_points]) - $xmin ] } {
          set bounds([list $canv xticspacing]) [ format %1.0e [ expr ( $xmax - $bounds([list $canv xorigin_points]) ) / $ntics_d ] ]
       } else {
          set bounds([list $canv xticspacing]) [ format %1.0e [ expr ( $bounds([list $canv xorigin_points]) - $xmin ) / $ntics_d ] ]
       }
    }
    if { $bounds([list $canv yfixtics])==0 } {
       if { [expr $ymax - $bounds([list $canv yorigin_points]) ] > [ expr $bounds([list $canv yorigin_points]) - $ymin ] } {
          set bounds([list $canv yticspacing]) [ format %1.0e [ expr ( $ymax - $bounds([list $canv yorigin_points]) ) / $ntics_d ] ]
       } else {
          set bounds([list $canv yticspacing]) [ format %1.0e [ expr ( $bounds([list $canv yorigin_points]) - $ymin ) / $ntics_d ] ]
       }
    }
}

proc gtPlot::findOrigin {canv axis} {
    variable padding
    variable bounds

    if { $axis == "x" } {
        set csize [expr [ winfo width $canv ] - 2.0*$padding ]
        set min $bounds([list $canv xmin])
        set max $bounds([list $canv xmax])
    } else {
        set csize [expr [ winfo height $canv ] - 2.0*$padding ]
        set min $bounds([list $canv ymin])
        set max $bounds([list $canv ymax])
    }
   
    if { $min<0 && $max>0 && $axis=="x" } {
        set bounds([list $canv xorigin_points]) 0
        set bounds([list $canv xorigin_canvas]) [ point2canvas $canv $axis 0.0 ]
    } elseif { $min<0 && $max>0 && $axis=="y" } {
        set bounds([list $canv yorigin_points]) 0
        set bounds([list $canv yorigin_canvas]) [ point2canvas $canv $axis 0.0 ]
    } elseif { $max<=0 && $axis=="x" } {
        set bounds([list $canv xorigin_points]) $max
        set bounds([list $canv xorigin_canvas]) [ expr $padding + $csize ]
    } elseif { $min>=0 && $axis=="x" } {
        set bounds([list $canv xorigin_points]) $min
        set bounds([list $canv xorigin_canvas]) $padding
    } elseif { $max<=0 && $axis=="y" } {
        set bounds([list $canv yorigin_points]) $max
        set bounds([list $canv yorigin_canvas]) $padding
    } elseif { $min>=0 && $axis=="y" } {
        set bounds([list $canv yorigin_points]) $min
        set bounds([list $canv yorigin_canvas]) [ expr $padding + $csize ]
    } else {
        tk_messageBox -icon error -type ok -title Message -message "Couldn't work out find origin"
    }
}

# These are all used to set up axis
proc gtPlot::setupGrid {tag} {
    variable data

    # Find the number of points in the grid in the x direction
    set kk 0
    while { [lindex $data([list $tag x]) $kk] == [lindex $data([list $tag x]) [expr $kk+1] ] } { incr kk }
    incr kk

    # Setup the grid
    set np $data([list $tag np])
    set i 0
    for { set nn 0 } { $nn<$np } { set nn [expr $nn + $kk] } {
        for { set j 0 } { $j<$kk } { incr j } {
            set data([list $tag fes $i $j]) \
                [list [lindex $data([list $tag x]) [expr $nn+$j] ] \
                      [lindex $data([list $tag y]) [expr $nn+$j] ] \
                      [lindex $data([list $tag z]) [expr $nn+$j] ] ]
        }
        incr i
    }

#    # We have now transfered the data onto a grid so we 
#    # can get rid of the data we read in
#    unset data([list $tag x])
#    unset data([list $tag y])
#    unset data([list $tag z])

    # Setup x and y grid size
    set data([list $tag y ngrid]) $kk
    set data([list $tag x ngrid]) $i
    # Setup x and y pixel widths
    set data([list $tag x pixwidth]) [expr [lindex $data([list $tag fes 1 0]) 0] - [lindex $data([list $tag fes 0 0]) 0] ]
    set data([list $tag y pixwidth]) [expr [lindex $data([list $tag fes 0 1]) 1] - [lindex $data([list $tag fes 0 0]) 1] ]

    # Some debugging - lets test how we do with the grid
    #for { set nn 0 } { $nn<$i } { incr nn } {
    #    for { set j 0 } { $j<$kk } { incr j } {
    #         puts "Grid point at gx, gy $nn $j x [lindex $data([list $tag fes $nn $j]) 0] y [lindex $data([list $tag fes $nn $j]) 1] z [lindex $data([list $tag fes $nn $j]) 2]"
    #    }
    #}
}

proc gtPlot::getRange {tag axis} {
    variable data
    set np $data([list $tag np])   ; # Retrieve number of data points

    # We have to check the data has the correct number of elements
    if { [llength $data([ list $tag $axis])]!=$np } {
       tk_messageBox -icon error -type ok -title Message -message "Not enough $axis data"
       return
    }

    for { set i 0} { $i<$np } { incr i } {
        set p [lindex $data([list $tag $axis]) $i]
        # Find maximum and minimum in the data
        if { $i == 0 } {
            set min $p;  set max $p;
        }
        if { $p > $max } { set max $p }
        if { $p < $min } { set min $p }
    }
   
    if { $min == $max } {
        # If the data is all on one point reset the x and y boundaries artificially
        set data([ list $tag $axis min ]) [ expr $min - 1 ]; set data([ list $tag $axis max ]) [ expr $max + 1 ];
    } else {
        # Otherwise change the boundaries by 5% of the range
        set data([ list $tag $axis min ]) [ expr $min - ( 0.05 *( $max - $min ) ) ]  ;
        set data([ list $tag $axis max ]) [ expr $max + ( 0.05 *( $max - $min ) ) ]  ;
    }
}

proc gtPlot::getAutoBound {canv searchtype axis} {
    variable data

    set flag 0
    foreach keyitem [array names data] {
       # Look for all the tag names
       if { [ lsearch $keyitem tag ] != -1 } {
          set tag $data($keyitem)
          # Check the data is displayed 
          if { [info exists data([ list $tag $canv displayed ])] } {
              if { $flag == 0 } {
                   set result $data([ list $tag $axis $searchtype ]) ; incr flag
              }
              if { $searchtype == "min" && $data([ list $tag $axis $searchtype ])<$result } {
                   set result $data([ list $tag $axis $searchtype ])
              } elseif { $searchtype == "max" && $data([ list $tag $axis $searchtype ])>$result } {
                   set result $data([ list $tag $axis $searchtype ])
              }
          }
       }
    }
    return $result
}

proc gtPlot::highlight {canv args} {
   variable data

   set tag [getargs $args "-datatag" "No datatag in call to gtPlot highlight"]
   set type [getargs $args "-type" "No type datatag in call to gtPlot highlight"]
   if { $type=="scaled" } {
     set data([list $tag arehighlists]) 1
     set data([list $tag highlight center]) [getargs $args "-center" "No center specified in call to gtPlot highlight"]
     set data([list $tag highlight isline]) [getargs $args "-lines" "Needs -lines flag in call to gtPlot highlight"]
     set data([list $tag highlight before]) [getargs $args "-before" "Need -before flag in call to gtPlot highlight"]
     set data([list $tag highlight afterc]) [getargs $args "-after" "Need -after flag in call to gtPlot highlight"]
     set data([list $tag trajcolors])       [getargs $args "-colors" "-colors flag missing from call to gtPlot highlight"]
     drawHighlight $canv $tag
  } elseif { $type=="onecolor" } {
     set listp [getargs $args "-listp" "-listp flag missing from call to gtPlot highlight"]
     set color [getargs $args "-color" "-color flag missing from call to gtPlot highlight"]
     drawMonoHighlight $canv $tag $listp $color
  } else {
     tk_messageBox -icon error -type ok -title Message -message "Invalid type spcified in call to gtPlot highlight"
     return
  }
}

proc gtPlot::drawMonoHighlight {canv tag listp color} {
   variable data
   variable bounds

   # This removes the old highlight
   foreach id [$canv find withtag mhighlight] {
      set taglist [$canv gettags $id]
      set listindex [lsearch -glob $taglist $tag*]
      set frametag [lindex $taglist $listindex]
      lassign [split $frametag :] foo num

      # Check if point is plotted (if so reconfigure and dtag
      if { [expr $num%$data([list $tag stride])] == 0 && [lsearch $taglist displayed]>=0 } {
         $canv itemconfigure $id -outline black -fill [lindex $data([list $tag color]) $num]
         $canv dtag $id mhighlight
      }
   }
   $canv delete mhighlight 

   # Now draw the highlight
   foreach val $listp {
     if { [expr $val%$data([list $tag stride])] == 0 } {
        set id [$canv find withtag $tag$val ]
        $canv itemconfigure $id -outline black -fill $color
        $canv raise $id displayed
        $canv addtag mhighlight withtag $tag$val
     } else {
        set x [lindex $data([list $tag x]) $val]
        set y [lindex $data([list $tag y]) $val]
        set w [ expr $data([ list $tag psize]) * ( 1  + [lindex $data([list $tag w]) $val] * $data([ list $tag pscale ]) ) ]
        if { $x >= $bounds([list $canv xmin]) && $x <= $bounds([list $canv xmax]) && $y >= $bounds([list $canv ymin]) && $y <= $bounds([list $canv ymax]) } {
            drawPoint $canv $x $y $w $color "black" "$tag$val $tag mhighlight"
        } else {
            drawPoint $canv $x $y $w $color "black" "$tag$val $tag mhighlight"
        }
     }
   }
}

proc gtPlot::drawHighlight {canv tag} {
   variable data
   variable bounds

   $canv delete beforeline
   $canv delete afterline

   # Get rid of highlighting of drawn points  
   # for super clever version this would only get rid of highlights on selected set
   foreach id [$canv find withtag highlighted] {
      # Recover the frame number this corresponds to
      set taglist [$canv gettags $id]
      set listindex [lsearch -glob $taglist $tag*]
      set frametag [lindex $taglist $listindex]
      lassign [split $frametag :] foo num

      # Check if point is plotted (if so reconfigure and dtag
      if { [expr $num%$data([list $tag stride])] == 0 && [lsearch $taglist displayed]>=0 } {
         $canv itemconfigure $id -outline black -fill [ lindex $data([list $tag color]) $num ]
         $canv dtag $id highlighted
      } 
   }

   # Get rid of old highlight stuff
   $canv delete highlighted

   # Setup a colorscale
   set nbefore [llength $data([list $tag highlight before])]
   set nafter [llength $data([list $tag highlight afterc])]
   if { $nbefore==0 && $nafter==0 } {
      set colormap [createColormap $data([list $tag trajcolors]) 2]
   } elseif { $nbefore>$nafter } {
      set colormap [createColormap $data([list $tag trajcolors]) [expr $nbefore + 1] ]
   } else {
      set colormap [createColormap $data([list $tag trajcolors]) [expr $nafter + 1] ]
   }

   # Add highlight on central point
   set val $data([list $tag highlight center])
   if { [expr $val%$data([list $tag stride])] == 0 } {
      set id [$canv find withtag $tag$val ]
      $canv itemconfigure $id -outline black -fill [lindex $colormap 0]
      $canv raise $id displayed
      $canv addtag highlighted withtag $tag$val
   } else {
      set x [lindex $data([list $tag x]) $val]
      set y [lindex $data([list $tag y]) $val]
      set w [ expr $data([ list $tag psize]) * ( 1  + [lindex $data([list $tag w]) $val] * $data([list $tag pscale]) ) ]
      if { $x >= $bounds([list $canv xmin]) && $x <= $bounds([list $canv xmax]) && $y >= $bounds([list $canv ymin]) && $y <= $bounds([list $canv ymax]) } {
          drawPoint $canv $x $y $w [lindex $colormap 0] "black" "$tag$val $tag highlighted" 
      } else {
      }
   }

   for { set i 0 } { $i<$nbefore } { incr i } {
      set val [ lindex $data([list $tag highlight before]) $i ]
      if { [expr $val%$data([list $tag stride])] == 0 } {
         set id [$canv find withtag $tag$val ]
         $canv itemconfigure $id -outline black -fill [lindex $colormap [expr $i+1] ]
         $canv raise $id displayed
         $canv addtag highlighted withtag $tag$val
      } else {
         set x [lindex $data([list $tag x]) $val]
         set y [lindex $data([list $tag y]) $val]
         set w [ expr $data([ list $tag psize]) * ( 1  + [lindex $data([list $tag w]) $val] * $data([list $tag pscale]) ) ]
         if { $x >= $bounds([list $canv xmin]) && $x <= $bounds([list $canv xmax]) && $y >= $bounds([list $canv ymin]) && $y <= $bounds([list $canv ymax]) } {
             drawPoint $canv $x $y $w [ lindex $colormap [expr $i+1] ] "black" "$tag$val $tag highlighted" 
         } else {
             drawPoint $canv $x $y $w [ lindex $colormap [expr $i+1] ] "black" "$tag$val $tag highlighted" 
         }
      }
      if { $data([list $tag highlight isline]) == 1 } {
         if { $i==0 } {
            set gal $data([list $tag highlight center])
         } else {
            set gal [ lindex $data([list $tag highlight before]) [expr $i-1] ]
         }
         set x1 [lindex $data([list $tag x]) $val]
         set y1 [lindex $data([list $tag y]) $val]
         set x2 [lindex $data([list $tag x]) $gal]
         set y2 [lindex $data([list $tag y]) $gal]
         drawLine $canv $x1 $y1 $x2 $y2 [lindex $colormap $i] 0.1 "$tag highlighted beforeline" 
         # set lid [$canv find withtag bl:$i]
         # $canv lower $lid displayed
      }
   }
   $canv lower beforeline displayed

   for { set i 0 } { $i<$nafter } { incr i } {
      set val [ lindex $data([list $tag highlight afterc]) $i ]
      if { [expr $val%$data([list $tag stride])] == 0 } {
         set id [$canv find withtag $tag$val ]
         $canv itemconfigure $id -outline black -fill [lindex $colormap [expr $i+1] ]
         $canv raise $id displayed
         $canv addtag highlighted withtag $tag$val
      } else {
         set x [lindex $data([list $tag x]) $val]
         set y [lindex $data([list $tag y]) $val]
         set w [ expr $data([ list $tag psize]) * ( 1  + [lindex $data([list $tag w]) $val] * $data([list $tag pscale]) ) ]
         if { $x >= $bounds([list $canv xmin]) && $x <= $bounds([list $canv xmax]) && $y >= $bounds([list $canv ymin]) && $y <= $bounds([list $canv ymax]) } {
             drawPoint $canv $x $y $w [lindex $colormap [expr $i+1] ] "black" "$tag$val $tag highlighted" 
         } else {
             drawPoint $canv $x $y $w [lindex $colormap [expr $i+1] ] "black" "$tag$val $tag highlighted" 
         }
      }
      if { $data([list $tag highlight isline]) == 1 } {
         if { $i==0 } {
            set gal $data([list $tag highlight center])
         } else {
            set gal [ lindex $data([list $tag highlight afterc]) [expr $i-1] ]
         }
         set x1 [lindex $data([list $tag x]) $val]
         set y1 [lindex $data([list $tag y]) $val]
         set x2 [lindex $data([list $tag x]) $gal]
         set y2 [lindex $data([list $tag y]) $gal]
         drawLine $canv $x1 $y1 $x2 $y2 [lindex $colormap $i] 0.1 "$tag highlighted afterline" 
         #set lid [$canv find withtag al:$i]
         #$canv lower $lid displayed
      }
   }
   $canv lower afterline displayed
}

##################################################### GRAPH PLOTTING TOOLS #############################################################

proc gtPlot::2dplotter {canv tag} {
    variable bounds
    variable data

    # Get the number of points
    set np $data([ list $tag np ])

    # Draw data points
    for { set i 0 } { $i<$np } { incr i } {
        if { [expr $i%$data([list $tag stride])] == 0 } {
           set x [lindex $data([list $tag x]) $i]
           set y [lindex $data([list $tag y]) $i]
           set c [lindex $data([list $tag color]) $i]
           set w [ expr $data([ list $tag psize]) * ( 1  + [lindex $data([list $tag w]) $i] * $data([list $tag pscale]) ) ]
           if { $data([list $tag drawpoints ]) == 1 } {
              if { $x >= $bounds([list $canv xmin]) && $x <= $bounds([list $canv xmax]) && $y >= $bounds([list $canv ymin]) && $y <= $bounds([list $canv ymax]) } {
                  if { $c=="white" } {
                     drawPoint $canv $x $y $w $c "black" "$tag$i $tag displayed" 
                  } else {
                     drawPoint $canv $x $y $w $c $c "$tag$i $tag displayed"
                  }
              } else { 
                  drawPoint $canv $x $y $w "white" "white" "$tag$i $tag"
              }
           }
           if { $data([list $tag drawline]) == 1 && $i>0 } {
                set x0 [lindex $data([list $tag x]) [expr $i - 1] ]
                set y0 [lindex $data([list $tag y]) [expr $i - 1] ]
                drawLine $canv $x0 $y0 $x $y $data([list $tag color]) 0.1 "$tag displayed dataline" 
           }
        }
    }
    set val 0

    ; # Extra work to deal with highlights 
    if { $data([list $tag arehighlights]) == 1 } { drawHighlights $canv $tag }   
}

proc gtPlot::barplotter {canv tag} {
    variable bounds
    variable data

    # Get the number of points
    set np $data([ list $tag np ])

    # The base of the bars is the x-axis
    set base $bounds([list $canv yorigin_points])

    # Now draw the bar chart
    for { set i 0 } { $i<$np } { incr i } {
        set x1 [lindex $data([list $tag x]) $i]
        if { [expr $i + 1]==$np } {
           set x2 [ expr 2*[lindex $data([list $tag x]) $i] - [lindex $data([list $tag x]) [expr $i - 1] ] ]
        } else {
           set x2 [lindex $data([list $tag x]) [expr $i + 1] ]
        }
        set y [lindex $data([list $tag y]) $i]
        drawLine $canv $x1 $base $x1 $y    $data([list $tag color]) 0.1 "$tag displayed bar" 
        drawLine $canv $x1 $y    $x2 $y    $data([list $tag color]) 0.1 "$tag displayed bar" 
        drawLine $canv $x2 $y    $x2 $base $data([list $tag color]) 0.1 "$tag displayed bar" 
   }
}

proc gtPlot::3dplotter {canv tag} {
    variable bounds
    variable data
    variable ncolors  
    
    # Check number of colors is specified properly
    if { ![ string is integer -strict $ncolors ] } {
       tk_messageBox -icon error -type ok -title Message -message "Number of colors in scale not specified"
       return
    }

    # Find the portion of the grid is inside the x and y range
    set gridbounds [ range2grid $tag x $bounds([list $canv xmin]) $bounds([list $canv xmax]) ]
    set xgridstart [lindex $gridbounds 0]   ; set xgridend [lindex $gridbounds 1]
    set gridbounds [ range2grid $tag y $bounds([list $canv ymin]) $bounds([list $canv ymax]) ]
    set ygridstart [lindex $gridbounds 0]   ; set ygridend [lindex $gridbounds 1]

    if { $xgridstart=="error" || $ygridstart=="error" } { return } 

    # If the user has not set z axis size using dialogue or is a fool then we get it automatically from the data
    if { $bounds([list $canv zmin])>=$bounds([list $canv zmax]) } {
       set bounds([list $canv zmin]) $data([list $tag z min]) ; set bounds([list $canv zmax]) $data([list $tag z max])
    }
    set zmin $bounds([list $canv zmin])  ;  set zmax $bounds([list $canv zmax])

    # Setup colorscale
    set colormap [ createColormap $data([list $tag fescolors]) $ncolors]
    set color_incr [ expr ($zmax - $zmin) / $ncolors ]
    # Draw the color scale
    drawScale $canv $tag $data([list $tag fescolors]) $zmin $zmax

    # Now plot the free energy surface
    for { set i $xgridstart } { $i<$xgridend } { incr i } {
        for { set j $ygridstart } { $j<$ygridend } { incr j } {
            # Find the point on the grid
            set x1 [ lindex $data([list $tag fes $i $j]) 0 ]
            set y1 [ lindex $data([list $tag fes $i $j]) 1 ]
            set x2 [ lindex $data([list $tag fes [expr $i + 1] [expr $j + 1] ]) 0 ]
            set y2 [ lindex $data([list $tag fes [expr $i + 1] [expr $j + 1] ]) 1 ]

            # Average the value of the fes at the four corners of this pixel
            set dat1 [ lindex $data([list $tag fes $i $j]) 2 ]
            set dat2 [ lindex $data([list $tag fes $i [expr $j + 1] ]) 2 ]
            set dat3 [ lindex $data([list $tag fes [expr $i + 1] $j])  2 ]
            set dat4 [ lindex $data([list $tag fes [expr $i + 1] [expr $j + 1] ]) 2 ]
            set av [ expr ( $dat1 + $dat2 + $dat3 + $dat4 ) / 4. ]
            if { $av<$zmin } { 
               set av $zmin 
            } elseif { $av>$zmax } {
               set av $zmax
            }

            # Convert this average value into a point on our discretized color scale
            set kk [ expr int( ( $av - $zmin ) / $color_incr ) ]
            set color [ lindex $colormap $kk ]

            # Draw this point
            drawPixel $canv $x1 $y1 $x2 $y2 $color "pixel surf $tag"
        }
    }

    # Put the fes at the bottom of the stack - this may yet go back GAT
#    $canv lower surf axis ;  $canv lower surf displayed
}

##################################################END GRAPH PLOTTING TOOLS #############################################################

proc gtPlot::zscale_window {canv wind args} {
    variable bounds

    frame $wind 
    labelframe $wind.zfram -relief ridge -bd 2 -text "Z-scale length" -padx 2m -pady 2m
    label $wind.zfram.minl -text "min" -anchor e
    entry $wind.zfram.min -width 5 -textvariable gtPlot::bounds([list $canv zmin])
    label $wind.zfram.maxl -text "max" -anchor e
    entry $wind.zfram.max -width 5 -textvariable gtPlot::bounds([list $canv zmax])
    grid $wind.zfram.minl -row 1 -column 1 -sticky e
    grid $wind.zfram.min  -row 1 -column 2 -sticky e
    grid $wind.zfram.maxl -row 2 -column 1 -sticky e
    grid $wind.zfram.max  -row 2 -column 2 -sticky e
    pack $wind.zfram -in $wind -side top -fill both

    labelframe $wind.sfram -relief ridge -bd 2 -text "Major scale tics" -padx 2m -pady 2m
    label $wind.sfram.spal -text "spacing" -anchor e
    entry $wind.sfram.spa -width 5 -textvariable gtPlot::bounds([list $canv zticspacing])
    label $wind.sfram.len -text "tic length" -anchor e
    spinbox $wind.sfram.en -textvariable gtPlot::bounds([list $canv zmajticl]) -from 1.0 -to 8.0 -increment 0.25 -width 5
    grid $wind.sfram.spal -row 1 -column 1 -sticky e
    grid $wind.sfram.spa -row 1 -column 2 -sticky e
    grid $wind.sfram.len  -row 2 -column 1 -sticky e
    grid $wind.sfram.en   -row 2 -column 2 -sticky e
    pack $wind.sfram -in $wind -side top -fill both

    labelframe $wind.mfram -relief ridge -bd 1 -text "Minor scale Tics" -padx 1m -pady 1m
    label $wind.mfram.nl -text "n. tics" -anchor e
    spinbox $wind.mfram.n -textvariable gtPlot::bounds([list $canv znmintic]) -from 0 -to 5 -increment 1 -width 5
    label $wind.mfram.lenl -text "tic length" -anchor e
    spinbox $wind.mfram.len -textvariable gtPlot::bounds([list $canv zminticl]) -from 1.0 -to 4.0 -increment 0.125 -width 5
    grid $wind.mfram.nl   -row 1 -column 1 -sticky e
    grid $wind.mfram.n    -row 1 -column 2 -sticky e
    grid $wind.mfram.lenl -row 2 -column 1 -sticky e
    grid $wind.mfram.len  -row 2 -column 2 -sticky e
    pack $wind.mfram -in $wind -side top -fill both
    return $wind
}

proc gtPlot::drawScale {canv tag fescolors zmin zmax} {
    variable ntics_d
    variable ncolors
    variable padding
    variable bounds
    variable font

    if { $bounds([list $canv zticspacing])==0 } {
       set ticspacing [ format %1.0e [ expr ( $zmax - $zmin ) / $ntics_d ] ] 
       set bounds([list $canv zticspacing]) $ticspacing
    }

    set ticspacing $bounds([list $canv zticspacing])
    set nmtics $bounds([list $canv znmintic])
    set majticl $bounds([list $canv zmajticl])
    set minticl $bounds([list $canv zminticl])

    # Check that ticspacing is a double
    if {![string is double -strict $ticspacing] } { return }

    $canv delete colorscale

    set colormap [ createColormap $fescolors $ncolors]
    set iheight [expr 0.5*[winfo height $canv] / $ncolors ]
    set sstart [expr 0.75*[winfo height $canv] ]

    # Set the x position of the colorscale
    set xl [expr [winfo width $canv] - 0.8*$padding ]
    set xu [expr [winfo width $canv] - 0.6*$padding ]
 
    # Draw the colorscale
    for {set i 0 } { $i<$ncolors } { incr i } {
       set yl [expr $sstart - $i*$iheight] ; set yu [expr $yl - $iheight]
       set color [ lindex $colormap $i ]
       $canv create rectangle $xl $yl $xu $yu -tag "$tag colorscale notfes" -fill $color -outline $color
    }

    # Create a rectangle around the color scale
    set yu [expr $sstart - $ncolors*$iheight ] ;
    $canv create rectangle $xl $sstart $xu $yu -tag "$tag colorscale notfes" -outline black

    # Now create tics
    set oaxis [expr [winfo width $canv] - 0.6*$padding ]
    set ntics [ expr int( ( $zmax - $zmin ) / $ticspacing) ]
    set realspacing [expr ( 0.5*[winfo height $canv] / ( $zmax - $zmin ) ) * $ticspacing ]

   # Draw tics
   for { set i 0 } { $i<=$ntics } { incr i } {
      set lab [ format %.4g [ expr $zmin + ( $i * $ticspacing ) ] ]
      set p [ expr $sstart - $i*$realspacing  ]
      $canv create line $oaxis $p [ expr $oaxis + $majticl ] $p -tag "$tag colorscale notfes"
      $canv create text [ expr $oaxis + 2 * $majticl ] $p -font $font -anchor w -text $lab -tag "$tag colorscale"

      # Draw minor tics 
      if { $i != 0 } {
         for { set j 1 } { $j <=$nmtics } { incr j } {
             set mp [ expr $p + ($j * $realspacing/($nmtics+1) ) ]
             $canv create line $oaxis $mp [ expr $oaxis + $minticl ] $mp -tag "$tag colorscale notfes"
         }
      }
   }
}

proc gtPlot::range2grid {tag axis min max} {
    variable data
    variable maxgrid

    if { $axis == "x" } {
         set n [expr $data([list $tag x ngrid]) - 1 ]   ; # Last point is skiped because of averages
         if { $n > $maxgrid } {
           tk_messageBox -icon error -type ok -title Message -message "Grid is too dense to plot in a reasonable time"
           return [list error]
         }
         set gmax $n ; set minf 0
         for { set i 1 } { $i<$n } { incr i } {
            if { $minf==0 && [lindex $data([list $tag fes $i 0]) 0 ] >= $min } { set gmin [expr $i-1]; set minf 1 }
            if { $minf==1 && [lindex $data([list $tag fes $i 0]) 0] >= $max } { set gmax [expr $i]; break; }
         }
    } else {
        set n [expr $data([list $tag y ngrid]) - 1 ]  ; # Last point is skiped because of averages
        if { $n > $maxgrid } {
          tk_messageBox -icon error -type ok -title Message -message "Grid is too dense to plot in a reasonable time"
          return [list error]
        }
        set gmax $n ; set minf 0
        for { set i 1 } { $i<$n } { incr i } {
           if { $minf==0 && [lindex $data([list $tag fes 0 $i]) 1 ] >= $min } { set gmin [expr $i-1]; set minf 1 }
           if { $minf==1 && [lindex $data([list $tag fes 0 $i]) 1] >= $max } { set gmax [expr $i]; break; }
        }
    }

    return [list $gmin $gmax]
}

proc gtPlot::drawPixel {canv x1 y1 x2 y2 color tag} {
   # Tranform the input coordinates onto the canvas   
   set x1 [ point2canvas $canv "x" $x1 ]
   set y1 [ point2canvas $canv "y" $y1 ]
   set x2 [ point2canvas $canv "x" $x2 ]
   set y2 [ point2canvas $canv "y" $y2 ]

   $canv create rectangle $x1 $y1 $x2 $y2 -tag "pixel $tag" -outline $color -fill $color
}

proc gtPlot::drawPoint {canv x y s fcolor ocolor tag} {
   set xx [ point2canvas $canv "x" $x ]
   set yy [ point2canvas $canv "y" $y ]

   $canv create rectangle [expr $xx - $s] [expr $yy - $s] [expr $xx + $s] [expr $yy + $s] -tag "pixel $tag" -fill $fcolor -outline $ocolor
}

proc gtPlot::drawLine {canv x1 y1 x2 y2 color width tag} { 
   set x1 [ point2canvas $canv "x" $x1 ]
   set y1 [ point2canvas $canv "y" $y1 ]
   set x2 [ point2canvas $canv "x" $x2 ]
   set y2 [ point2canvas $canv "y" $y2 ]

   $canv create line $x1 $y1 $x2 $y2 -tag "line $tag" -fill $color -width $width
}

# This writes our little coordinate thing
proc gtPlot::writeCoord {canv x y} {
   variable data
   variable font

   # Delete the old label   
   $canv delete whereIam

   # A check that there is some data plotted
   set flag 0
   foreach keyitem [array names data] {
      # Look for all the tag names
      if { [ lsearch $keyitem tag ] != -1 } {
           set tag $data($keyitem)
           # Check if data is displayed
           if { [info exists data([ list $tag $canv displayed ])] } { set flag 1 }
      }
   }
   if { $flag == 0 } { return }

   # This gets our coordinates
   set px [format %.4g [canvas2point $canv x $x ] ]
   set py [format %.4g [canvas2point $canv y $y ] ]

   # Write the position on the canvas
   $canv create text [ expr [ winfo width $canv ] - 5] 15 -font $font -anchor se -text "($px, $py)" -tag "whereIam"
}

# Routines to convert between canvas coordinates and real coordinates
proc gtPlot::canvas2point {canv d x} {
   variable padding
   variable bounds
  
   if { $d == "x" } {
      set s [ expr [ winfo width $canv ] - 2.0*$padding ]
      return [ expr ( ( $x - $padding ) / $s - 0.5 ) * $bounds([list $canv xrange]) + $bounds([list $canv xcenter]) ]
   } else {
      set s [ expr [ winfo height $canv ] - 2.0*$padding ]
      return [ expr ( 1.0 - ( $x - $padding ) / $s - 0.5 ) * $bounds([list $canv yrange]) + $bounds([list $canv ycenter]) ]
   }
}

proc gtPlot::point2canvas {canv d x} {
   variable padding
   variable bounds

   if { $d == "x" } {
      set s [ expr [ winfo width $canv ] - 2.0*$padding ]
      return [ expr $s * ( 0.5 + ( $x - $bounds([list $canv xcenter]) ) / $bounds([list $canv xrange]) ) + $padding ]
   } else {
      set s [ expr [ winfo height $canv ] - 2.0*$padding ]
      return [ expr $s * ( 0.5 - ( $x - $bounds([list $canv ycenter]) ) / $bounds([list $canv yrange]) ) + $padding ]
   }
}

# Stuff for the zoom tool
proc gtPlot::ZoomStart {canv x y} {
   variable zoomID
   set zoomID(0) [$canv create rect $x $y $x $y -dash { 4 4 } -tag zoom ]
   set zoomID(1) [$canv create rect [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] -tag zoom ]
   set zoomID(2) [$canv create rect [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] -tag zoom ]
   set zoomID(xmin) $x     ;   set zoomID(xmax) $x
   set zoomID(ymin) $y     ;   set zoomID(ymax) $y                                              
}   

proc gtPlot::ZoomMove {canv x y} {
   variable zoomID
   $canv coords $zoomID(0) [lreplace [$canv coords $zoomID(0)] 2 3 $x $y]
   $canv coords $zoomID(2) [lreplace [$canv coords $zoomID(2)] 0 3 [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] ]
   set zoomID(xmax) $x    ;    set zoomID(ymin) $y
}

proc gtPlot::ZoomEnd {canv x y} {
   variable zoomID
   variable bounds

   $canv delete zoom

   # Get the coordinates of the zoom box
   set bounds([list $canv txmin]) [canvas2point $canv "x" $zoomID(xmin) ]
   set bounds([list $canv tymin]) [canvas2point $canv "y" $zoomID(ymin) ]
   set bounds([list $canv txmax]) [canvas2point $canv "x" $zoomID(xmax) ]
   set bounds([list $canv tymax]) [canvas2point $canv "y" $zoomID(ymax) ]
   set bounds([list $canv zmin]) 0  ; set bounds([list $canv zmax]) 0

   # Set the axis
   setAxis $canv $bounds([list $canv txmin]) $bounds([list $canv txmax]) $bounds([list $canv tymin]) $bounds([list $canv tymax])
   # Draw the graph
   drawGraph $canv
}

# Stuff for the lasso tool
proc gtPlot::LassoStart {canv x y} {
   variable lassoID
   set lassoID(0) [$canv create rect $x $y $x $y -dash { 4 4 } -tag lasso ]
   set lassoID(1) [$canv create rect [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] -tag lasso ]
   set lassoID(2) [$canv create rect [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] -tag lasso ]
   set lassoID(xmin) $x     ;   set lassoID(xmax) $x
   set lassoID(ymin) $y     ;   set lassoID(ymax) $y
}

proc gtPlot::LassoMove {canv x y} {
   variable lassoID
   $canv coords $lassoID(0) [lreplace [$canv coords $lassoID(0)] 2 3 $x $y]
   $canv coords $lassoID(2) [lreplace [$canv coords $lassoID(2)] 0 3 [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] ]
   set lassoID(xmax) $x    ;    set lassoID(ymin) $y
}

proc gtPlot::LassoEnd {canv x y} {
   variable lassoID
   variable data
   variable doneLasso

   $canv delete lasso       ; # Delete the lasso stuff

   # Get the proper coordinates of the lasso
   set lassoID(xmin) [canvas2point $canv "x" $lassoID(xmin)]
   set lassoID(ymin) [canvas2point $canv "y" $lassoID(ymin)] 
   set lassoID(xmax) [canvas2point $canv "x" $lassoID(xmax)]
   set lassoID(ymax) [canvas2point $canv "y" $lassoID(ymax)]

   # Now make a list of all the points in the lasso
   foreach keyitem [array names data] {
      if { [lsearch $keyitem tag ] != -1 } {
          set tag $data($keyitem)
          set lassoID([list $tag points]) {} 
          set np $data([ list $tag np ])
          for { set i 0 } { $i<$np } { incr i } {
              set x [lindex $data([list $tag x]) $i]
              set y [lindex $data([list $tag y]) $i]
              if { $x >= $lassoID(xmin) && $x <= $lassoID(xmax) && $y >= $lassoID(ymin) && $y <= $lassoID(ymax) } {
                 lappend lassoID([list $tag points]) $i
              }
          }
      }
   }

   set doneLasso 1
   return
}

proc gtPlot::getLasso {args} {
   variable data
   variable lassoID
   set tag [getargs $args "-datatag" "No datatag specified in call to gtPlot getLasso"] 
   set dtag [getargs $args "-data" "No -data flag specified in call to gtPlot getLasso"]

   # Check there is data plotted with this tag
   set flag 0
   foreach keyitem [array names data] {
      if { [lsearch $keyitem tag ] != -1 } {
         if { $data($keyitem)==$tag } { set flag 1 ; break }
      }
   }
   if { $flag==0 } { return }

   if { $dtag=="indices" } {
      return $lassoID([list $tag points])
   } elseif { $dtag=="zdata" } {
      if { $data([list $tag is3d])==0 } {
          set zlist {}
          foreach i $lassoID([list $tag points]) { lappend zlist [lindex $data([list $tag y]) $i] }
          return $zlist
      } else {
          set zlist {}
          foreach i $lassoID([list $tag points]) { lappend zlist [lindex $data([list $tag z]) $i] }
          return $zlist
      }
   } elseif { $dtag=="box" } {
      return [list $lassoID(xmin) $lassoID(ymin) $lassoID(xmax) $lassoID(ymax)]
   }
}

proc gtPlot::printCanvas {canv args} {
   set filen [tk_getSaveFile -initialfile "plumedplot.ps" -title "Print graph to file" -parent $canv \
    -filetypes [list {{Postscript files} {.ps}} {{All files} {*}}]]

   if { $filen!="" } { $canv postscript -file $filen }
   return
}

# This gets arguments or returns a warning
proc gtPlot::getargs {arglist tag {warning 0} } {
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

# This builds a menu to pick a colormap
proc gtPlot::colorMenu {window args} {
    set var [getargs $args "-textvariable" "-textvariable unset in call to gtPlot colorMenu"]

    menubutton $window -relief raised -bd 2 -direction flush -textvariable $var -menu $window.menu
    menu $window.menu
    $window configure -state disabled
    $window.menu delete 0 end
    $window configure -state normal
    $window.menu add radiobutton -label jet -value jet -variable $var
    $window configure -state normal
    $window.menu add radiobutton -label hsv -value hsv -variable $var 
    $window configure -state normal
    $window.menu add radiobutton -label hot -value hot -variable $var 
    $window configure -state normal
    $window.menu add radiobutton -label cool -value cool -variable $var 
    $window configure -state normal
    $window.menu add radiobutton -label grey -value grey -variable $var
    return $window
}

# This sets up our colormaps
proc gtPlot::createColormap {colorscale ncolors} {
  switch -- $colorscale {
    hsv {
        set hueStart     0.0
        set hueEnd     240.0
        set colorMap   {}

        for {set i 0} {$i <= $ncolors} {incr i} {
            set dh [expr {($hueStart - $hueEnd) / ($ncolors - 1)}]
            set hue  [expr {$hueStart - ($i * $dh)}]
            if {$hue < 0.0} {
                set hue  [expr {360.0 + $hue}]
            }
            set rgbList [Hsv2rgb $hue 1.0 1.0]
            set r    [expr {int([lindex $rgbList 0] * 65535)}]
            set g    [expr {int([lindex $rgbList 1] * 65535)}]
            set b    [expr {int([lindex $rgbList 2] * 65535)}]

            set color  [format "#%.4x%.4x%.4x" $r $g $b]
            lappend colorMap $color
        }
        return $colorMap
    }
    hot {
        set colorMap {}
        set nc1          [expr {int($ncolors * 0.33)}]
        set nc2          [expr {int($ncolors * 0.67)}]

        for {set i 0} {$i <= $ncolors} {incr i} {

            if {$i <= $nc1} {
                set fval  [expr { double($i) / (double($nc1)) } ]
                set r     [expr {int($fval * 65535)}]
                set g     0
                set b     0
            } else {
                if {$i <= $nc2} {
                    set fval  [expr { double($i-$nc1) / (double($nc2-$nc1)) } ]
                    set r     65535
                    set g     [expr {int($fval * 65535)}]
                    set b     0
                } else {
                    set fval  [expr { double($i-$nc2) / (double($ncolors-$nc2)) } ]
                    set r     65535
                    set g     65535
                    set b     [expr {int($fval * 65535)}]
                }
            }
            set color  [format "#%.4x%.4x%.4x" $r $g $b]
            lappend colorMap $color
        }
        return $colorMap
    }
    cool {
        set colorMap {}
        for {set i 0} {$i <= $ncolors} {incr i} {

            set fval  [expr { double($i) / (double($ncolors)-1) } ]
            set val   [expr { 1.0 - $fval }]

            set r    [expr {int($fval * 65535)}]
            set g    [expr {int($val * 65535)}]
            set b    65535

            set color  [format "#%.4x%.4x%.4x" $r $g $b]
            lappend colorMap $color
        }
        return $colorMap
    }
    grey {
        set colorMap {}
        for {set i 0} {$i <= $ncolors} {incr i} {

            set fval  [expr { double($i) / (double($ncolors)-1) } ]
            set val  [expr {0.4 + (0.5 * $fval) }]

            set r    [expr {int($val * 65535)}]
            set g    [expr {int($val * 65535)}]
            set b    [expr {int($val * 65535)}]

            set color  [format "#%.4x%.4x%.4x" $r $g $b]
            lappend colorMap $color
        }
        return $colorMap
    }
    jet {
         set hueStart   240.0
         set hueEnd       0.0
         set colorMap   {}

         for {set i 0} {$i <= $ncolors} {incr i} {

            set dh [expr {($hueStart - $hueEnd) / ($ncolors - 1)}]
            set hue  [expr {$hueStart - ($i * $dh)}]
            if {$hue < 0.0} {
                set hue  [expr {360.0 + $hue}]
            }
            set rgbList [Hsv2rgb $hue 1.0 1.0]
            set r    [expr {int([lindex $rgbList 0] * 65535)}]
            set g    [expr {int([lindex $rgbList 1] * 65535)}]
            set b    [expr {int([lindex $rgbList 2] * 65535)}]

            set color  [format "#%.4x%.4x%.4x" $r $g $b]
            lappend colorMap $color
         }
         return $colorMap
    }
  }
}

proc gtPlot::Hsv2rgb {h s v} {
    set v [expr {double($v)}]
    set r [set g [set b 0.0]]
    if {$h == 360} { set h 0 }
    # if you feed the output of rgb2hsv back into this
    # converter, h could have the value -1 for
    # grayscale colors.  Set it to any value in the
    # valid range.
    if {$h == -1} { set h 0 }
    set h [expr {$h/60}]
    set i [expr {int(floor($h))}]
    set f [expr {$h - $i}]
    set p1 [expr {$v*(1-$s)}]
    set p2 [expr {$v*(1-($s*$f))}]
    set p3 [expr {$v*(1-($s*(1-$f)))}]
    switch -- $i {
        0 { set r $v  ; set g $p3 ; set b $p1 }
        1 { set r $p2 ; set g $v  ; set b $p1 }
        2 { set r $p1 ; set g $v  ; set b $p3 }
        3 { set r $p1 ; set g $p2 ; set b $v  }
        4 { set r $p3 ; set g $p1 ; set b $v  }
        5 { set r $v  ; set g $p1 ; set b $p2 }
    }
    return [list $r $g $b]
}
