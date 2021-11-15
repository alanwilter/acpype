/**************************************************************
*
* Progress Circle 1.1.0
* Updated support for bootstrap
*
**************************************************************/

( function( $ ){
	var ProgressCircle = function( element, options ){

		var settings          = $.extend( {}, $.fn.progressCircle.defaults, options );
		var thicknessConstant = 0.02;
		var nRadian           = 0;

		computePercent();
		setThickness();

		var border      = ( settings.thickness * thicknessConstant ) + 'em';
		var offset      = ( 1 - thicknessConstant * settings.thickness * 2 ) + 'em';
		var circle      = $( element );
		var progCirc    = circle.find( '.prog-circle' );
		var circleDiv   = progCirc.find( '.bar' );
		var circleSpan  = progCirc.children( '.percenttext' );
		var circleFill  = progCirc.find( '.fill' );
		var circleSlice = progCirc.find( '.slice' );

		if ( settings.nPercent == 0 ) {
			circleSlice.hide();
		} else {
			resetCircle();
			transformCircle( nRadians, circleDiv );
		}
		setBorderThickness();
		updatePercentage();
		setCircleSize();

		function computePercent () {
			settings.nPercent > 100 || settings.nPercent < 0 ? settings.nPercent = 0 : settings.nPercent;
			nRadians = ( 360 * settings.nPercent ) / 100;
		}

		function setThickness () {
			if ( settings.thickness > 10 ) {
				settings.thickness = 10;
			} else if ( settings.thickness < 1 ) {
				settings.thickness = 1;
			} else {
				settings.thickness = Math.round( settings.thickness );
			}
		}

		function setCircleSize ( ) {
			progCirc.css( 'font-size', settings.circleSize + 'px' );
		}

		function transformCircle ( nRadians, cDiv ) {
			var rotate = "rotate(" + nRadians + "deg)";
	    cDiv.css({
	      "-webkit-transform" : rotate,
	      "-moz-transform"    : rotate,
	      "-ms-transform"     : rotate,
	      "-o-transform"      : rotate,
	      "transform"         : rotate
	    });
	    if( nRadians > 180 ) {
		    transformCircle( 180, circleFill );
		    circleSlice.addClass( ' clipauto ');
			}
		}

		function setBorderThickness () {
			progCirc.find(' .slice > div ').css({
				'border-width' : border,
				'width'        : offset,
				'height'       : offset
			})
			if( settings.barOverride ) {
				progCirc.find(' .slice > div.bar ').css({
					'width'        : '1em',
					'height'       : '1em'
				})
			}
			progCirc.find('.after').css({
				'top'    : border,
				'left'   : border,
				'width'  : offset,
				'height' : offset
			})
		}

		function resetCircle () {
			circleSlice.show();
			circleSpan.text( '' );
			circleSlice.removeClass( 'clipauto' )
			transformCircle( 20, circleDiv );
			transformCircle( 20, circleFill );
			return this;
		}

		function updatePercentage () {
			settings.showPercentText && circleSpan.text( settings.nPercent + '%' );
		}
	};

	$.fn.progressCircle = function( options ) {
		return this.each( function( key, value ){
      var element = $( this );
      if ( element.data( 'progressCircle' ) ) {
      	var progressCircle = new ProgressCircle( this, options );
      	return element.data( 'progressCircle' );
      }
    	$( this ).append( '<div class="prog-circle">' +
	                 			'	<div class="percenttext"> </div>' +
	                 			'	<div class="slice">' +
	                 			'		<div class="bar"> </div>' +
	                 			'		<div class="fill"> </div>' +
	                 			'	</div>' +
	                 			'	<div class="after"> </div>' +
	                 			'</div>');
      var progressCircle = new ProgressCircle( this, options );
      element.data( 'progressCircle', progressCircle );
  	});
	};

	$.fn.progressCircle.defaults = {
		nPercent        : 50,
		showPercentText : true,
		circleSize      : 100,
		thickness       : 3,
		barOverride     : true
	};

})( jQuery );