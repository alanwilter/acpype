// Core javascript helper functions

// basic browser identification & version
var isOpera = (navigator.userAgent.indexOf("Opera")>=0) && parseFloat(navigator.appVersion);
var isIE = ((document.all) && (!isOpera)) && parseFloat(navigator.appVersion.split("MSIE ")[1].split(";")[0]);

// Cross-browser event handlers.
function addEvent(obj, evType, fn) {
    if (obj.addEventListener) {
        obj.addEventListener(evType, fn, false);
        return true;
    } else if (obj.attachEvent) {
        var r = obj.attachEvent("on" + evType, fn);
        return r;
    } else {
        return false;
    }
}

function removeEvent(obj, evType, fn) {
    if (obj.removeEventListener) {
        obj.removeEventListener(evType, fn, false);
        return true;
    } else if (obj.detachEvent) {
        obj.detachEvent("on" + evType, fn);
        return true;
    } else {
        return false;
    }
}

// quickElement(tagType, parentReference, textInChildNode, [, attribute, attributeValue ...]);
function quickElement() {
    var obj = document.createElement(arguments[0]);
    if (arguments[2] != '' && arguments[2] != null) {
        var textNode = document.createTextNode(arguments[2]);
        obj.appendChild(textNode);
    }
    var len = arguments.length;
    for (var i = 3; i < len; i += 2) {
        obj.setAttribute(arguments[i], arguments[i+1]);
    }
    arguments[1].appendChild(obj);
    return obj;
}

// ----------------------------------------------------------------------------
// Cross-browser xmlhttp object
// from http://jibbering.com/2002/4/httprequest.html
// ----------------------------------------------------------------------------
var xmlhttp;
/*@cc_on @*/
/*@if (@_jscript_version >= 5)
    try {
        xmlhttp = new ActiveXObject("Msxml2.XMLHTTP");
    } catch (e) {
        try {
            xmlhttp = new ActiveXObject("Microsoft.XMLHTTP");
        } catch (E) {
            xmlhttp = false;
        }
    }
@else
    xmlhttp = false;
@end @*/
if (!xmlhttp && typeof XMLHttpRequest != 'undefined') {
  xmlhttp = new XMLHttpRequest();
}

// ----------------------------------------------------------------------------
// Find-position functions by PPK
// See http://www.quirksmode.org/js/findpos.html
// ----------------------------------------------------------------------------
function findPosX(obj) {
    var curleft = 0;
    if (obj.offsetParent) {
        while (obj.offsetParent) {
            curleft += obj.offsetLeft - ((isOpera) ? 0 : obj.scrollLeft);
            obj = obj.offsetParent;
        }
        // IE offsetParent does not include the top-level
        if (isIE && obj.parentElement){
            curleft += obj.offsetLeft - obj.scrollLeft;
        }
    } else if (obj.x) {
        curleft += obj.x;
    }
    return curleft;
}

function findPosY(obj) {
    var curtop = 0;
    if (obj.offsetParent) {
        while (obj.offsetParent) {
            curtop += obj.offsetTop - ((isOpera) ? 0 : obj.scrollTop);
            obj = obj.offsetParent;
        }
        // IE offsetParent does not include the top-level
        if (isIE && obj.parentElement){
            curtop += obj.offsetTop - obj.scrollTop;
        }
    } else if (obj.y) {
        curtop += obj.y;
    }
    return curtop;
}

//-----------------------------------------------------------------------------
// Date object extensions
// ----------------------------------------------------------------------------
Date.prototype.getCorrectYear = function() {
    // Date.getYear() is unreliable --
    // see http://www.quirksmode.org/js/introdate.html#year
    var y = this.getYear() % 100;
    return (y < 38) ? y + 2000 : y + 1900;
}

Date.prototype.getTwoDigitMonth = function() {
    return (this.getMonth() < 9) ? '0' + (this.getMonth()+1) : (this.getMonth()+1);
}

Date.prototype.getTwoDigitDate = function() {
    return (this.getDate() < 10) ? '0' + this.getDate() : this.getDate();
}

Date.prototype.getTwoDigitHour = function() {
    return (this.getHours() < 10) ? '0' + this.getHours() : this.getHours();
}

Date.prototype.getTwoDigitMinute = function() {
    return (this.getMinutes() < 10) ? '0' + this.getMinutes() : this.getMinutes();
}

Date.prototype.getTwoDigitSecond = function() {
    return (this.getSeconds() < 10) ? '0' + this.getSeconds() : this.getSeconds();
}

Date.prototype.getISODate = function() {
    return this.getCorrectYear() + '-' + this.getTwoDigitMonth() + '-' + this.getTwoDigitDate();
}

Date.prototype.getHourMinute = function() {
    return this.getTwoDigitHour() + ':' + this.getTwoDigitMinute();
}

Date.prototype.getHourMinuteSecond = function() {
    return this.getTwoDigitHour() + ':' + this.getTwoDigitMinute() + ':' + this.getTwoDigitSecond();
}

// ----------------------------------------------------------------------------
// String object extensions
// ----------------------------------------------------------------------------
String.prototype.pad_left = function(pad_length, pad_string) {
    var new_string = this;
    for (var i = 0; new_string.length < pad_length; i++) {
        new_string = pad_string + new_string;
    }
    return new_string;
}

// ----------------------------------------------------------------------------
// Get the computed style for and element
// ----------------------------------------------------------------------------
function getStyle(oElm, strCssRule){
    var strValue = "";
    if(document.defaultView && document.defaultView.getComputedStyle){
        strValue = document.defaultView.getComputedStyle(oElm, "").getPropertyValue(strCssRule);
    }
    else if(oElm.currentStyle){
        strCssRule = strCssRule.replace(/\-(\w)/g, function (strMatch, p1){
            return p1.toUpperCase();
        });
        strValue = oElm.currentStyle[strCssRule];
    }
    return strValue;
}

/* gettext identity library */
function gettext(msgid) { return msgid; }
function ngettext(singular, plural, count) { return (count == 1) ? singular : plural; }
function gettext_noop(msgid) { return msgid; }
function interpolate(fmt, obj, named) {
  if (named) {
    return fmt.replace(/%\(\w+\)s/g, function(match){return String(obj[match.slice(2,-2)])});
  } else {
    return fmt.replace(/%s/g, function(match){return String(obj.shift())});
  }
}

// Finds all fieldsets with class="collapse", collapses them, and gives each
// one a "Show" link that uncollapses it. The "Show" link becomes a "Hide"
// link when the fieldset is visible.

function findForm(node) {
    // returns the node of the form containing the given node
    if (node.tagName.toLowerCase() != 'form') {
        return findForm(node.parentNode);
    }
    return node;
}

var CollapsedFieldsets = {
    collapse_re: /\bcollapse\b/,   // Class of fieldsets that should be dealt with.
    collapsed_re: /\bcollapsed\b/, // Class that fieldsets get when they're hidden.
    collapsed_class: 'collapsed',
    init: function() {
        var fieldsets = document.getElementsByTagName('fieldset');
        var collapsed_seen = false;
        for (var i = 0, fs; fs = fieldsets[i]; i++) {
            // Collapse this fieldset if it has the correct class, and if it
            // doesn't have any errors. (Collapsing shouldn't apply in the case
            // of error messages.)
            // wb104: changed this 1 Nov 2011 so that can see options in forms
            //if (fs.className.match(CollapsedFieldsets.collapse_re) && !CollapsedFieldsets.fieldset_has_errors(fs)) {
            if (false) {
                collapsed_seen = true;
                // Give it an additional class, used by CSS to hide it.
                fs.className += ' ' + CollapsedFieldsets.collapsed_class;
                // (<a id="fieldsetcollapser3" class="collapse-toggle" href="#">Show</a>)
                var collapse_link = document.createElement('a');
                collapse_link.className = 'collapse-toggle';
                collapse_link.id = 'fieldsetcollapser' + i;
                collapse_link.onclick = new Function('CollapsedFieldsets.show('+i+'); return false;');
                collapse_link.href = '#';
                collapse_link.title = 'More or Less Options'
                collapse_link.innerHTML = gettext('&raquo;');
                var h2 = fs.getElementsByTagName('a')[0];
                //h2.appendChild(document.createTextNode(' ('));
                //h2.appendChild(collapse_link);
                h2.parentNode.replaceChild(collapse_link, h2);
                //h2.appendChild(document.createTextNode(')'));
            }
        }
        if (collapsed_seen) {
            // Expand all collapsed fieldsets when form is submitted.
            addEvent(findForm(document.getElementsByTagName('fieldset')[0]), 'submit', function() { CollapsedFieldsets.uncollapse_all(); });
        }
    },
    fieldset_has_errors: function(fs) {
        // Returns true if any fields in the fieldset have validation errors.
        var divs = fs.getElementsByTagName('span');
        for (var i=0; i<divs.length; i++) {
            if (divs[i].className.match(/\berrors\b/)) {
                return true;
            }
        }
        return false;
    },
    show: function(fieldset_index) {
        var fs = document.getElementsByTagName('fieldset')[fieldset_index];
        // Remove the class name that causes the "display: none".
        fs.className = fs.className.replace(CollapsedFieldsets.collapsed_re, '');
        // Toggle the "Show" link to a "Hide" link
        var collapse_link = document.getElementById('fieldsetcollapser' + fieldset_index);
        collapse_link.onclick = new Function('CollapsedFieldsets.hide('+fieldset_index+'); return false;');
        collapse_link.innerHTML = gettext('&laquo;');
    },
    hide: function(fieldset_index) {
        var fs = document.getElementsByTagName('fieldset')[fieldset_index];
        // Add the class name that causes the "display: none".
        fs.className += ' ' + CollapsedFieldsets.collapsed_class;
        // Toggle the "Hide" link to a "Show" link
        var collapse_link = document.getElementById('fieldsetcollapser' + fieldset_index);
        collapse_link.onclick = new Function('CollapsedFieldsets.show('+fieldset_index+'); return false;');
        collapse_link.innerHTML = gettext('&raquo;');
    },

    uncollapse_all: function() {
        var fieldsets = document.getElementsByTagName('fieldset');
        for (var i=0; i<fieldsets.length; i++) {
            if (fieldsets[i].className.match(CollapsedFieldsets.collapsed_re)) {
                CollapsedFieldsets.show(i);
            }
        }
    }
}

addEvent(window, 'load', CollapsedFieldsets.init);
