var pageQuery=function () { 
  var query_string = {};
  var query = window.location.search.substring(1);
  var vars = query.split("&");
  for (var i=0;i<vars.length;i++) {
    var pair = vars[i].split("=");
    	// If first entry with this name
    if (typeof query_string[pair[0]] === "undefined") {
      query_string[pair[0]] = pair[1];
    	// If second entry with this name
    } else if (typeof query_string[pair[0]] === "string") {
      var arr = [ query_string[pair[0]], pair[1] ];
      query_string[pair[0]] = arr;
    	// If third or later entry with this name
    } else {
      query_string[pair[0]].push(pair[1]);
    }
  } 
    return query_string;
} ();

if (!pageQuery.page) pageQuery.page="main"; //default

function httpGet(theUrl, entity, append, silent)
{
   entity = typeof entity !== "undefined" ? entity : "";
   append = typeof append !== "undefined" ? append: false;
   silent = typeof silent !== "undefined" ? silent: true;

   // Reads a file from URL and returns its content as a string
   var xmlhttp;
   if (window.XMLHttpRequest)
   {// code for IE7+, Firefox, Chrome, Opera, Safari
      xmlhttp=new XMLHttpRequest();
   }
   else
   {// code for IE6, IE5
      xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
   }

   if (entity!="") {
      xmlhttp.onreadystatechange=function() {
         if (append) prefix=document.getElementById(entity).innerHTML+" "; else prefix=""; 
         if (xmlhttp.readyState == 4 && (xmlhttp.status == 200 || (xmlhttp.status == 0 && xmlhttp.responseText))) {
            document.getElementById(entity).innerHTML = prefix + xmlhttp.responseText;      
         } else {
            if (!silent) document.getElementById(entity).innerHTML = prefix + "<b> Error. HTTP " + xmlhttp.status + " </b>";           
         }      
      }
   }      
   try {
      xmlhttp.overrideMimeType("text/plain; charset=utf-8");
      xmlhttp.open("GET", theUrl, entity!="" );    
      xmlhttp.send();    
   } catch(err) {  if (!silent) { document.write("Could not load resources from "+theUrl); } throw(err); }

   return xmlhttp.responseText;   
}

function httpExists(theUrl)
{
   var xmlhttp;
   if (window.XMLHttpRequest)
   {// code for IE7+, Firefox, Chrome, Opera, Safari
      xmlhttp=new XMLHttpRequest();
   }
   else
   {// code for IE6, IE5
      xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
   }
   try{
      xmlhttp.overrideMimeType("text/plain; charset=utf-8");
      xmlhttp.open("GET", theUrl, false );    
      xmlhttp.send();    
      if (xmlhttp.readyState == 4 && (xmlhttp.status == 200 || (xmlhttp.status == 0 && xmlhttp.responseText))) return true;
      else return false;
   } catch(err) {return false;}
}

function httpInclude(theUrl)
{  // Reads from URL and writes to document in place 
   document.write(httpGet(theUrl));
}

function composeHeader(pageName)
{
   try { httpGet("heads/base.html", "header", true); }
   catch(err) { document.write("Could not load page header."); }

   try { httpGet("heads/"+pageName+".html", "header", true); }
   catch(err) { httpGet("heads/default.html", "header", true) }   
}

function mkMenu(item,label)
{
   var menu=document.getElementById("menu");
   menu.innerHTML = menu.innerHTML + " <a class='menu " + (item==pageQuery.page?"sel":"")+ "' href='index.html?page="+item+"'>"+label+"</a>";
}

function mkSubmenu(item,label)
{
   var submenu=document.getElementById("submenu");
   submenu.innerHTML = submenu.innerHTML + "<a class='submenu "+(item==pageQuery.psub?"subsel":"submenu")+"' href='index.html?page="+pageQuery.page+"&amp;psub="+item+"'>"+label+"</a>";
}


function toggle2(showHideDiv, switchTextDiv) {
	var ele = document.getElementById(showHideDiv);
	var text = document.getElementById(switchTextDiv);
	if(ele.style.display == "block") {
      ele.style.display = "none";
		text.innerHTML = "expand +";
  	}
	else {
		ele.style.display = "block";
		text.innerHTML = "collapse -";
	}
}

function toggle(showHideDiv, switchTextDiv) {
	var ele = document.getElementById(showHideDiv);
	var text = document.getElementById(switchTextDiv);
	if(ele.style.display == "block") {
      ele.style.display = "none";
		text.innerHTML = text.innerHTML.replace("-&nbsp;","+&nbsp;");
  	}
	else {
		ele.style.display = "block";
		text.innerHTML = text.innerHTML.replace("+&nbsp;","-&nbsp;");
	}
}   

//<![CDATA[
var txtControl='[+-]';
var txtCtrlExpanded='<span style="float:right; display:block; clear:none;">&nbsp;collapse&nbsp;-</span>';
var txtCtrlCollapsed='<span style="float:right; display:block; clear:none;">&nbsp;expand&nbsp;+</span>'; 
//]]>   
function safe_toggle(showHideDiv, switchTextDiv) {
	var ele = document.getElementById(showHideDiv);
	var text = document.getElementById(switchTextDiv);
	text.innerHTML = text.innerHTML.replace(txtControl,txtCtrlExpanded)   	
	if(ele.style.display == "block") {
      ele.style.display = "none";
		text.innerHTML = text.innerHTML.replace(txtCtrlExpanded,txtCtrlCollapsed);
  	}
	else {
		ele.style.display = "block";
		text.innerHTML = text.innerHTML.replace(txtCtrlCollapsed,txtCtrlExpanded);
	}
}    
