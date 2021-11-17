/********************Stuff for the Language Filter**************************/
function processClicks() {
//  alert("You released the mouse button.");
	makeInvisible( getElt( "FILTER_POPUP" ) );
	return true;
}

function languageFilterPopup()
{
	divObj = getElt( "FILTER_POPUP" );
	forWhat = getElt( "FILTER" );

	var forWhatRect = new Rectangle;
	getElementSize( forWhat, forWhatRect );
	var popupRect = new Rectangle;
	divObj.style.top = (forWhatRect.top + forWhatRect.height + 2)+"px";
	divObj.style.left = "10px";
	makeVisible( divObj );

	// Now arrange for it to get popped down when the user clicks
	//   outside the box.
	if (window.Event) // Navigator 4.0x
	  document.captureEvents(Event.MOUSEUP);
	document.onmouseup = processClicks;
}

var filter = "all";

function initFilter()
{
	var allcookies = document.cookie;
	var pos = allcookies.indexOf( "languageFilter=" );
	if ( pos >= 0 ) {
		var start = pos+15;
		var end = allcookies.indexOf(";", start);
		if ( end == -1 ) end = allcookies.length;
		filter = allcookies.substr(start,end-start);
		if ( filter != "all" && filter != "cpp" &&
			 filter != "cs" && filter != "vb" &&
			 filter != "js" )
		  filter = "all";
	}
}

initFilter();

function showlang( name, bShow )
{
	var a = document.getElementsByTagName("div");

	for ( var i = 0; i < a.length; ++i )
	{
		if ( a(i).NAME == name )
		{
			a(i).style.display = (bShow ? "" : "none");
		}
	}
}

function showProperLanguage()
{
	showlang( "Decl_Cpp", ( filter == "all" || filter == "cpp" ) );
	showlang( "Decl_Cs", ( filter == "all" || filter == "cs" ) );
	showlang( "Decl_VisualBasic", ( filter == "all" || filter == "vb" ) );
	showlang( "Decl_JScript", ( filter == "all" || filter == "js" ) );
}

function setFilter( f ) {
	filter = f;
	document.cookie = "languageFilter="+f;

	showProperLanguage();
}


/*********************************** Rectangle *****************************/

function Rectangle( left, top, width, height )
{
	this.left = left;
	this.top = top;
	this.width = width;
	this.height = height;
}

/*************************** Browser-dependent stuff ***********************/

function makeVisible( obj )
{
	obj.style.visibility = "visible";
}

function makeInvisible( obj )
{
	obj.style.visibility = "hidden";
}

// This function populates r (a rectangle) with the visible area of
//  the window or frame
function getVisibleArea( r )
{
	var o3_frame = window;

	r.top = (navigator.family == 'ie4') ? document.body.scrollTop : window.pageYOffset;
	r.left = (navigator.family == 'ie4') ? document.body.scrollLeft : window.pageXOffset;

	r.width = o3_frame.document.body.clientWidth;
	if (!r.width) r.width = o3_frame.innerWidth; // was screwed in mozilla, fixed now?
	if (!r.width) r.width = o3_frame.outerWidth;

	r.height = o3_frame.document.body.clientHeight;
	if (!r.height) r.height = o3_frame.innerHeight; // was screwed in mozilla, fixed now?
	if (!r.height) r.height = o3_frame.outerHeight;
}

function getRealLeft(el) {
    xPos = el.offsetLeft;
    tempEl = el.offsetParent;
    while (tempEl != null) {
        xPos += tempEl.offsetLeft;
        tempEl = tempEl.offsetParent;
    }
    return xPos;
}

function getRealTop(el) {
    yPos = el.offsetTop;
    tempEl = el.offsetParent;
    while (tempEl != null) {
        yPos += tempEl.offsetTop;
        tempEl = tempEl.offsetParent;
    }
    return yPos;
}

// This gets the extent of the given object
function getElementSize( obj, r )
{
	if ( typeof obj.offsetLeft != 'undefined' ) {
		r.left = getRealLeft( obj );
		r.top = getRealTop( obj );
		r.width = obj.offsetWidth;
		r.height = obj.offsetHeight;
	}
	else {
		alert( "Got no idea how to figure out the element size." );
	}
}

function setElementPos( obj, r )
{
	obj.style.left = r.left + "px";
	obj.style.top = r.top + "px";
	obj.style.width = r.width + "px";
	obj.style.height = r.height + "px";
}


/********************************* Popup Support ***************************/

var timerHandle = null;
var verticalPad = 5;
var divObj = null;

function showit()
{
	makeVisible( divObj );
}

function getElt( o )
{
	if ( typeof(document.getElementById) == 'undefined' )
	  return document.all( o );
	else
	  return document.getElementById( o );
}

function popup( forWhat, html )
{
	if ( document.body.insertAdjacentHTML ) {
		document.body.insertAdjacentHTML( 'AfterBegin', '<DIV ID=popup CLASS=popup>xyz</DIV>' );
		divObj = getElt( "popup" );
	}
	else {
		divObj = document.createElement( "DIV" );
		divObj.className = "popup";
		document.body.appendChild( divObj );
	}

	divObj.innerHTML = html.replace(/\|q/g, "'").replace(/\|d/g, '"').replace(/\|n/g, '\n').replace(/\|\|/g, "|");

	clearTimeout( timerHandle );

	var forWhatRect = new Rectangle;
	getElementSize( forWhat, forWhatRect );
	var popupRect = new Rectangle;
	getElementSize( divObj, popupRect );
	var screenSize = new Rectangle;
	getVisibleArea( screenSize );

	// Now what we need to do is decide where to put the
	// popup.  Ideally, we center it right above forWhat.
	// If there's not enough space above it, we then try
	// for below it.  If there's not enough space either
	// way, we go for the one that offers the most room.

	// We never set the width; we just assume that it will
	// never be made wider than the screen.

	var idealLeft = forWhatRect.left - ( popupRect.width - forWhatRect.width ) / 2;
	if ( idealLeft < screenSize.left )
	  idealLeft = screenSize.left;
	else if ( idealLeft + popupRect.width > screenSize.left + screenSize.width )
	  idealLeft = screenSize.left + screenSize.width - popupRect.width;

	var idealTop = forWhatRect.top - popupRect.height - verticalPad;
	if ( idealTop < screenSize.top )
	  idealTop = forWhatRect.top + forWhatRect.height + verticalPad;


	popupRect.left = idealLeft;
	popupRect.top = idealTop;

	setElementPos( divObj, popupRect );

	timerHandle = setTimeout( "showit();", 200 );
}

function clearpopup()
{
	clearTimeout( timerHandle );
	timerHandle = null;

	if ( divObj != null ) {
		makeInvisible( divObj );

		var screenSize = new Rectangle;
		getVisibleArea( screenSize );
		divObj = null;
	}
}

//------------------------------------------------------------------------------------------------------------------

function InitDerivation()
{
	var allcookies = document.cookie;
	var pos = allcookies.indexOf( "whichDerivation=" );
	if ( pos == -1 )
	  ShowBriefDerivation();
	else {
		var start = pos+16;
		var end = allcookies.indexOf(";", start);
		if ( end == -1 ) end = allcookies.length;
		var value = allcookies.substr(start,end);
		if ( value == "full" )
		  ShowFullDerivation();
		else
		  ShowBriefDerivation();
	}
}

function ShowBriefDerivation()
{
	var innerHTML = getElt('BriefDerivation').innerHTML;
	if ( navigator.org == 'netscape' )
	  innerHTML = unescape( innerHTML );
	getElt('DualDerivation').innerHTML = innerHTML;
	document.cookie = "whichDerivation=brief";
}

function ShowFullDerivation()
{
	var innerHTML = getElt('FullDerivation').innerHTML;
	if ( navigator.org == 'netscape' )
	  innerHTML = unescape( innerHTML );
	getElt('DualDerivation').innerHTML = innerHTML;
	document.cookie = "whichDerivation=full";
}


function openTOC( url )
{
	if ( isFrameView() )
	  window.open( url, "TOC" );
	else
	  window.open( url, "TOC",
		"toolbar=no,width=350,height=600,directories=no,status=no,scrollbars=yes,resize=yes,menubar=no" );
}

function openBody(url)
{
	if ( url != "" ) {
		window.open( url, "Body" );
	}
}

function isFrameView()
{
	return window.parent != null && window.parent != window.self;
}

function toggleFrames()
{
	if ( isFrameView() ) {
		window.parent.open( unescape(document.URL), "_self" );
	}
	else {
		// open the frameset with the current page
		//  in the body window.
		var querystring = "bodyURL=" + document.URL;
		if ( window.contentsWindowURL != null )
		  querystring += ",contentsURL=" + window.contentsWindowURL;
		window.parent.open( "contentsf.html" + "?" + querystring, "_self" );
	}
}

function accessError(message,url,line) {
	window.close();
	return true;
}

function toggleTOC(menuName)
{
	var menu = getElt('sub_'+menuName);
	var vSrc = getElt(menuName);

	if ( menu.style.display == "none" ) {
		// Change from + to -
		vSrc.src = "true.gif";
		menu.style.display = "";
	}
	else {
		vSrc.src = "false.gif";
		menu.style.display = "none";
	}
}

function getArgs() {
	var args = new Object();
	var query = location.search.substring(1);
	var pairs = query.split(",");
	for ( var i = 0; i < pairs.length; i++ ) {
		var pos = pairs[i].indexOf('=');
		if ( pos == -1 ) continue;
		var argname = pairs[i].substring(0,pos);
		var value = pairs[i].substring(pos+1);
		args[argname] = unescape(value);
	}
	return args;
}


function ViewInVB( objname )
{
	if ( VSControl.IsVBRunning() ) {
		wentok = VSControl.ShowInVB( objname );
		if ( !wentok ) {
			alert( "Something went wrong -- probably VB is not open to the right project" );
		}
	}
	else {
	    alert( "VB or the VB plug-in is not running" );
	}
}

function ViewInDevStudio( project, file, line, defaultPath )
{
	var helpURL = unescape( document.URL.replace(/[^\/\\]*$/, "" ) ) + "GoToSourceHelp.html";

	if ( !VSControl.ShowInDevStudio( project, file, line, defaultPath, helpURL ) ) {
		// Could just be the user canceled the action.
	}
}

function TOCHeight(depth) {
  if ( depth < 2 )
    return 22;
  else if ( depth < 4 )
    return 20;
  else
    return 18;
}


args = getArgs();

window.focus();





////////////////////////////////////////////////////////////////////////////////////
//  Begin UA.JS

// ua.js - Detect Browser
// Requires JavaScript 1.1
/*
The contents of this file are subject to the Netscape Public
License Version 1.1 (the "License"); you may not use this file
except in compliance with the License. You may obtain a copy of
the License at http://www.mozilla.org/NPL/

Software distributed under the License is distributed on an "AS
IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or
implied. See the License for the specific language governing
rights and limitations under the License.

The Initial Developer of the Original Code is Bob Clary.

Contributor(s): Bob Clary, Original Work, Copyright 1999-2000
                Bob Clary, Netscape Communications, Copyright 2001

Alternatively, the contents of this file may be used under the
terms of the GNU Public License (the "GPL"), in which case the
provisions of the GPL are applicable instead of those above.
If you wish to allow use of your version of this file only
under the terms of the GPL and not to allow others to use your
version of this file under the NPL, indicate your decision by
deleting the provisions above and replace them with the notice
and other provisions required by the GPL.  If you do not delete
the provisions above, a recipient may use your version of this
file under either the NPL or the GPL.
*/

// work around bug in xpcdom Mozilla 0.9.1
window.saveNavigator = window.navigator;

// Handy functions
function noop() {}
function noerror() { return true; }

function defaultOnError(msg, url, line)
{
	// customize this for your site
	if (top.location.href.indexOf('_files/errors/') == -1)
		top.location = '/evangelism/xbProjects/_files/errors/index.html?msg=' + escape(msg) + '&url=' + escape(url) + '&line=' + escape(line);
}

// Display Error page... 
// XXX: more work to be done here
//
function reportError(message)
{
	// customize this for your site
	if (top.location.href.indexOf('_files/errors/') == -1)
		top.location = '/evangelism/xbProjects/_files/errors/index.html?msg=' + escape(message);
}

function pageRequires(cond, msg, redirectTo)
{
	if (!cond)
	{
		msg = 'This page requires ' + msg;
		top.location = redirectTo + '?msg=' + escape(msg);
	}
	// return cond so can use in <A> onclick handlers to exclude browsers
	// from pages they do not support.
	return cond;
}

function detectBrowser()
{
	var oldOnError = window.onerror;
	var element = null;
	
	window.onerror = defaultOnError;

	navigator.OS		= '';
	navigator.version	= 0;
	navigator.org		= '';
	navigator.family	= '';

	var platform;
	if (typeof(window.navigator.platform) != 'undefined')
	{
		platform = window.navigator.platform.toLowerCase();
		if (platform.indexOf('win') != -1)
			navigator.OS = 'win';
		else if (platform.indexOf('mac') != -1)
			navigator.OS = 'mac';
		else if (platform.indexOf('unix') != -1 || platform.indexOf('linux') != -1 || platform.indexOf('sun') != -1)
			navigator.OS = 'nix';
	}

	var i = 0;
	var ua = window.navigator.userAgent.toLowerCase();
	
	if (ua.indexOf('opera') != -1)
	{
		i = ua.indexOf('opera');
		navigator.family	= 'opera';
		navigator.org		= 'opera';
		navigator.version	= parseFloat('0' + ua.substr(i+6), 10);
	}
	else if ((i = ua.indexOf('msie')) != -1)
	{
		navigator.org		= 'microsoft';
		navigator.version	= parseFloat('0' + ua.substr(i+5), 10);
		
		if (navigator.version < 4)
			navigator.family = 'ie3';
		else
			navigator.family = 'ie4'
	}
	else if (typeof(window.controllers) != 'undefined' && typeof(window.locationbar) != 'undefined')
	{
		i = ua.lastIndexOf('/')
		navigator.version = parseFloat('0' + ua.substr(i+1), 10);
		navigator.family = 'gecko';

		if (ua.indexOf('netscape') != -1)
			navigator.org = 'netscape';
		else if (ua.indexOf('compuserve') != -1)
			navigator.org = 'compuserve';
		else
			navigator.org = 'mozilla';
	}
	else if ((ua.indexOf('mozilla') !=-1) && (ua.indexOf('spoofer')==-1) && (ua.indexOf('compatible') == -1) && (ua.indexOf('opera')==-1)&& (ua.indexOf('webtv')==-1) && (ua.indexOf('hotjava')==-1))
	{
	    var is_major = parseFloat(navigator.appVersion);
    
		if (is_major < 4)
			navigator.version = is_major;
		else
		{
			i = ua.lastIndexOf('/')
			navigator.version = parseFloat('0' + ua.substr(i+1), 10);
		}
		navigator.org = 'netscape';
		navigator.family = 'nn' + parseInt(navigator.appVersion);
	}
	else if ((i = ua.indexOf('aol')) != -1 )
	{
		// aol
		navigator.family	= 'aol';
		navigator.org		= 'aol';
		navigator.version	= parseFloat('0' + ua.substr(i+4), 10);
	}

	navigator.DOMCORE1	= (typeof(document.getElementsByTagName) != 'undefined' && typeof(document.createElement) != 'undefined');
	navigator.DOMCORE2	= (navigator.DOMCORE1 && typeof(document.getElementById) != 'undefined' && typeof(document.createElementNS) != 'undefined');
	navigator.DOMHTML	= (navigator.DOMCORE1 && typeof(document.getElementById) != 'undefined');
	navigator.DOMCSS1	= ( (navigator.family == 'gecko') || (navigator.family == 'ie4') );

	navigator.DOMCSS2   = false;
	if (navigator.DOMCORE1)
	{
		element = document.createElement('p');
		navigator.DOMCSS2 = (typeof(element.style) == 'object');
	}

	navigator.DOMEVENTS	= (typeof(document.createEvent) != 'undefined');
	navigator.CollapseTOC = ( navigator.family == 'ie4' || navigator.family == 'gecko' );

	window.onerror = oldOnError;
}

detectBrowser();

/************************ overloadExpander *******************************/
function hasAttribute( o, attribute, value )
{
	return ( o.attributes != null &&
		     o.attributes[attribute] != null &&
			 o.attributes[attribute].nodeValue == value );
}

function isOverloadDetails( o )
{
	return hasAttribute( o, "class", "OverloadDetails" )
}

function findChild( p, testFun )
{
	if ( testFun( p ) )
		return p;
	for ( var o = p.firstChild; o != null; o = o.nextSibling )
	{
		var c = findChild( o, testFun );
		if ( c != null )
		  return c;
	}
	return null;
}

function toggleShowOverload(img)
{
	var p = img;
	while ( img != null && !hasAttribute( p, "class", "OverloadContainer" ) )
		p = p.parentNode;

	o = findChild( p, isOverloadDetails );
	if ( o.style.display == "" )
	{
		o.style.display = "none";
		img.src = "overload_plus.gif";
	}
	else
	{
		o.style.display = "";
		img.src = "overload_minus.gif";
	}
}


function getElementsByClassName(oElm, strTagName, strClassName)
{
/*
    Written by Jonathan Snook, http://www.snook.ca/jonathan
    Add-ons by Robert Nyman, http://www.robertnyman.com
*/
    var arrElements = (strTagName == "*" && document.all)? document.all : oElm.getElementsByTagName(strTagName);
    var arrReturnElements = new Array();
    strClassName = strClassName.replace(/\-/g, "\\-");
    var oRegExp = new RegExp("(^|\\s)" + strClassName + "(\\s|$)");
    var oElement;
    for(var i=0; i<arrElements.length; i++){
        oElement = arrElements[i];      
        if(oRegExp.test(oElement.className)){
            arrReturnElements.push(oElement);
        }   
    }
    return (arrReturnElements)
}

var showAllOverloads = false;

function showAllOverloadDetails( setting )
{
	var aElt = getElt( "OverloadShowAll" );
	if ( aElt != null )
	{
		aElt.innerHTML = (showAllOverloads ? "[Hide Overload Details]" : "[Show All Overload Details]");
		// If aElt is null, this isn't an overloaded function

		var allExpanders = getElementsByClassName( document, "img", "OverloadExpander" );
		for ( var i = 0; i < allExpanders.length; ++i )
		{
			var img = allExpanders[i];
			img.src = (showAllOverloads ? "overload_minus.gif" : "overload_plus.gif");
		}

		var allDetails = getElementsByClassName( document, "div", "OverloadDetails" );
		for ( var i = 0; i < allDetails.length; ++i )
		{
			var d = allDetails[i];
			d.style["display"] = (showAllOverloads ? "" : "none");
		}
	}
}

function toggleShowAllOverloads()
{
	showAllOverloads = !showAllOverloads;
	showAllOverloadDetails();
}

function initShowAllOverloads()
{
	var allcookies = document.cookie;
	var pos = allcookies.indexOf( "showAllOverloads=" );
	if ( pos >= 0 ) {
		var start = pos+17;
		var end = allcookies.indexOf(";", start);
		if ( end == -1 ) end = allcookies.length;
		showAllOverloads = (allcookies.substr(start,end) == "true");
	}
	showAllOverloadDetails(showAllOverloads);
}

initShowAllOverloads();

function onLoad()
{
	showAllOverloadDetails( showAllOverloads );
	showProperLanguage();
}
