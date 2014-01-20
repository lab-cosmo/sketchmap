package provide progressbar 1.0

namespace eval ::progressbar {
# Routines we export
  namespace export create
  namespace export config
  namespace export increment
# The internal variables
  variable width     ; # The width of the progress bar
  variable height    ; # The height of the progress bar
  variable ntasks    ; # Number of tasks that the bar must perform
  variable itask     ; # The current number of tasks that have been performed
}

proc progressbar::create {window args} {
  variable width
  variable height

  set width($window) [getargs $args "-width" "In call to create progress bar -width was unspecified"]
  set height($window) [getargs $args "-height" "In call to create progress bar -height was unspecified"]

  # Create a canvas
  eval { canvas $window -highlightthickness 0 -borderwidth 0 -background white } \
    -width $width($window)  -height $height($window) 

  return $window
}

proc progressbar::config {window args} {
  variable ntasks
  variable itask

  $window delete task
  update idletasks
  set itask($window) 0
  set ntasks($window) [getargs $args "-ntasks" "In call to configure progress bar -ntasks was unspecified"]
}

proc progressbar::increment {window args} {
  variable itask
  variable width
  variable height
  variable ntasks

  set x1 [ expr $itask($window) * $width($window) / $ntasks($window) ]
  set x2 [ expr ( $itask($window) + 1 ) * $width($window) / $ntasks($window) ] 

  $window create rectangle $x1 0 $x2 $height($window) -tag "task $itask($window)" -outline black -fill black
  update idletasks
  incr itask($window)
}

proc progressbar::getargs {arglist tag {warning 0} } {
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
