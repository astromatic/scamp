<!DOCTYPE html>
<html>
	<head>
		<title>Processing Summary</title>
		<link 
			rel="stylesheet" 
			href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" 
			integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" 
			crossorigin="anonymous">
		<link 
			rel="stylesheet" 
			href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css" 
			integrity="sha384-rHyoN1iRsVXV4nD0JutlnGaslCJuC7uwjduW9SVrLvRYooPp2bWYgmgJQIXwl/Sp" 
			crossorigin="anonymous">
		<link 
			rel="stylesheet" 
			href="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css" />
		<link 
			rel="stylesheet" 
			href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.12.1/bootstrap-table.min.css">
		<link 
			rel="stylesheet" 
			type="text/css" 
			href="https://cdn.datatables.net/1.10.18/css/dataTables.bootstrap.min.css">
		<style type="text/css">
			@-ms-viewport	 { width: device-width; }
			@-o-viewport	  { width: device-width; }
			@viewport		 { width: device-width; }
			canvas {
				-moz-user-select: none;
				-webkit-user-select: none;
				-ms-user-select: none;
			}
			.chart-container {
				width: 350px;
				margin-left: 40px;
				margin-right: 40px;
				margin-bottom: 40px;
			}
			.container {
				display: flex;
				flex-direction: row;
				flex-wrap: wrap;
				justify-content: center;
			}
		</style>
		<!--[if lt IE 9]>
			<script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
			<script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
		<![endif]-->
	</head>
	<body>
		<nav class="navbar navbar-inverse " role="navigation">
			<div class="container-fluid">
				<div class="navbar-header">
					<a 
						style="padding-top:4px;" 
						class="navbar-brand" 
						target="_blank" href="http://astromatic.net">
						<img 
							style="height: 50px; margin-top:0px;" 
							alt="astromatic" 
							src="http://astromatic.net/xsl/astromatic.png">
					</a>
				</div>
				<div class="nav navbar-nav navbar-right">
					<p class="navbar-text pull-right">
						<strong>
							<a target="_blank" 
								href="https://github.com/astromatic/scamp">
								<span id="soft" 
									class="label label-primary" 
									style="font-size:50;">
								</span>
							</a>
						</strong> 
						completed on <strong><span id="date"></span></strong>
						at <strong><span id="time"></span></strong> 
						using <strong><span id="nthreads"></span> </strong>
						threads (run time: <strong><span id="runtime"></span></strong>)
						started by user <strong><span id="username"></span></strong> 
						in <strong><span id="rundir"></span></strong>.
					</p>
				</div> <!-- navright -->
			</div> <!--container-fluid-->
		</nav> <!-- nav -->

		<div class="container-fluid role="main"">
			<!-- WARNINGS TABLE -->
			<div class="row">
				<div class="col-xs-12 col-sm-12 col-md-12">
					<div id ="warningDiv" class="alert alert-warning">
						<button 
							class="btn btn-warning" 
							type="button" data-toggle="collapse" 
							data-target="#warningsTableCollapse" 
							aria-expanded="false" 
							aria-controls="warningsTableCollapse">
							<strong>Warnings</strong><span id="warningsShort"></span>
						</button>
						<div id="warningsTableCollapse" class="collapse panel-body table-responsive">
							<table id="warningsTable" class="table table-responsive table-striped">
								<thead>
									<tr>
										<th>Date</th>
										<th>Time</th>
										<th>Message</th>
									</tr>
								</thead>
								<tbody>
								</tbody>
							</table>
						</div>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->



			<div class="row">
				<!--  HISTOGRAMS -->
				<div class="col-xs-12 col-sm-12 col-md-12 col-lg-7">
					<div class="row">
						<div class="col-xs-12 col-sm-12 col-md-12">
							<div class="btn-group" role="group" aria-label="Basic example">
								<div class="btn-group" role="group" aria-label="Basic example">
									<div class="dropdown">
										<button 
											type="button" 
											class="btn btn-primary dropdown-toggle" 
											data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
											<span id="plot1_btn1">Table</span>
											<span class='caret'></span>
										</button>
										<ul id="plot1_tbl" class="dropdown-menu table-dropdown">
										</ul>
									</div>
								</div>
								<div class="btn-group" role="group" aria-label="Basic example">
									<div class="dropdown">
										<button 
											type="button" 
											class="btn btn-primary dropdown-toggle" 
											data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
											<span id ="plot1_btn2">Column</span>
											<span class='caret'></span>
										</button>
										<ul id="plot1_cols" class="dropdown-menu column-dropdown">
										</ul>
									</div>
								</div>
							</div>
							<div id="plot1">
							</div>
						</div>
					</div>
					<div class="row">
						<div class="col-xs-12 col-sm-12 col-md-12">
							<div class="btn-group" role="group" aria-label="Basic example">
								<div class="btn-group" role="group" aria-label="Basic example">
									<div class="dropdown">
										<button 
											type="button" 
											class="btn btn-primary dropdown-toggle" 
											data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
											<span id ="plot2_btn1">Table</span>
											<span class='caret'></span>
										</button>
										<ul id="plot2_tbl" class="dropdown-menu table-dropdown">
										</ul>
									</div>
								</div>
								<div class="btn-group" role="group" aria-label="Basic example">
									<div class="dropdown">
										<button 
											type="button" 
											class="btn btn-primary dropdown-toggle" 
											data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
											<span id ="plot2_btn2">Column</span>
											<span class='caret'></span>
										</button>
										<ul id="plot2_cols" class="dropdown-menu column-dropdown">
										</ul>
									</div>
								</div>
							</div>
							<div id="plot2">
							</div>
						</div>
					</div>
				</div> <!-- end col -->
				<!-- ALADIN -->
				<div class="col-xs-12 col-sm-12 col-md-12 col-lg-5">
					<div class="panel panel-default">
						<div class="panel-heading text-center">Catalogs Footprint</div>
						<div id="aladin-lite-div" class="panel-body" style="width: 100%; height:600px;"></div>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->



			<!-- FIELD TABLE -->
			<div class="row">
				<div class="col-md-12">
					<div id="inputFiles" class="panel panel-default">
						<div class="panel-heading text-left">
							<button 
								class="btn btn-info" 
								type="button" 
								data-toggle="collapse" 
								data-target="#fieldsTableCollapse" 
								aria-expanded="false" 
								aria-controls="fieldsTableCollapse">
								<strong>Summary Table on Input Files</strong> (Fields)
							</button>
						</div>
						<div id="fieldsTableCollapse" class="collapse panel-body table-responsive">
							<table 
								tabindex='1'
								id="fieldsTable" 
								class="table table-hover table-bordered table-striped table-sm" 
								cellspacing="0">
								<thead>
									<tr>
										<th data-sortable="true">#</th>
										<th data-sortable="true" data-field="Filename">Filename</th>
										<th>Identifier</th>
										<th>Next</th>
										<th>Ndet</th>
										<th>Flags</th>
										<th>G</th>
										<th>A</th>
										<th>P</th>
										<th>Date</th>
										<th>Exposure Time</th>
										<th>Air Mass</th>
										<th>Right Ascension</th>
										<th>Declination</th>
										<th>Radius</th>
										<th>Pixel scale</th>
										<th class="showmatch">&#916;Pixel scale</th>
										<th class="showmatch">&#916;Position Angle</th>
										<th class="showmatch">A/S contrast</th>
										<th class="showmatch">&#916;X</th>
										<th class="showmatch">&#916;Y</th>
										<th class="showmatch">X/Y contrast</th>
										<th>&#967;2int</th>
										<th>&#967;2int High S/N</th>
										<th>&#967;2ref</th>
										<th>&#967;2ref High S/N</th>
										<th>Mag &#916;SP</th>
									</tr>
								</thead>
								<tbody>
								</tbody>
							</table>
						</div>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->



			<!-- GROUPS TABLE -->
			<div class="row">
				<div class="col-md-12">
					<div id="groupProperties" class="panel panel-default">
						<div class="panel-heading text-left">
							<button 
								class="btn btn-info" 
								type="button" 
								data-toggle="collapse" 
								data-target="#groupPropertiesCollapse" 
								aria-expanded="false" 
								aria-controls="groupPropertiesCollapse">
								<strong>Group Properties</strong> (Fgroups)
							</button>
						</div>
					<div id="groupPropertiesCollapse" class="collapse panel-body table-responsive">
						<table id="groupsTable" class="table table-bordered table-striped">
							<thead>
								<tr>
									<th>Group name</th>
									<th class="showplot">Group Plot</th>
									<th>Index</th>
									<th>NFields</th>
									<th>Right Ascension</th>
									<th>Declination</th>
									<th>Pixel scale</th>
									<th>Maximum radius</th>
									<th>Astrom. Ref. Catalog</th>
									<th>Astrom. Ref. Band</th>
									<th class="showplot">&#967;2 Plot</th>
									<th>Astrom. &#963;int</th>
									<th>Astrom. &#961;int</th>
									<th>Astrom. &#967;2int</th>
									<th>Astrom. N2int</th>
									<th>Astrom. &#963;int High S/N</th>
									<th>Astrom. &#961;int High S/N</th>
									<th>Astrom. &#947;2int High S/N</th>
									<th>Astrom. N2int High S/N</th>
									<th class="showplot">Astrom. 1-D Int. Error Plot</th>
									<th class="showplot">Astrom. 2-D Int. Error Plot</th>
									<th>Astrom. &#916;RA ref, &#916; DEC ref</th>
									<th>Astrom. &#963;ref</th>
									<th>Astrom. &#961;ref</th>
									<th>Astrom. &#967;2ref</th>
									<th>Astrom. Nref</th>
									<th>Astrom. &#916;RA ref, &#916; DEC ref High S/N</th>
									<th>Astrom. &#963;ref High S/N</th>
									<th>Astrom. &#961;ref High S/N</th>
									<th>Astrom. &#947;2ref High S/N</th>
									<th>Astrom. Nref High S/N</th>
									<th class="showplot">Astrom. 1-D Ref. Error Plot</th>
									<th class="showplot">Astrom. 2-D Ref. Error Plot</th>
									<th>Photom. instruments</th>
									<th>Photom. &#963;int</th>
									<th>Photom. &#967;2int</th>
									<th>Photom. Nint</th>
									<th>Photom. &#963; High S/N</th>
									<th>Astrom. &#947;2int High S/N</th>
									<th>Photom. Nint High S/N</th>
									<th>Photom. &#963;ref</th>
									<th>Photom. &#967;2ref</th>
									<th>Photom. Nref</th>
									<th>Photom. &#963;ref High S/N</th>
									<th>Photom. &#947;2ref High S/N</th>
									<th>Photom. Nref High S/N</th>
									<th class="showplot">Photom. Internal Error Plot</th>
								</tr>
							</thead>
							<tbody>
							</tbody>
						</table>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->



			<div class="row">
				<!-- ASTROINSTRUMENTS TABLE -->
				<div class="col-md-6"> 
					<div id="astroInstruments" class="panel panel-default">
						<div class="panel-heading text-left">
							<button 
								class="btn btn-info" 
								type="button" 
								data-toggle="collapse" 
								data-target="#astrometricInstrumentsCollapse" 
								aria-expanded="false" 
								aria-controls="astrometricInstrumentsCollapse">
								<strong>Astrometric Instruments</strong> (AstroInstruments)
							</button>
						</div>
						<div id="astrometricInstrumentsCollapse" class="collapse panel-body table-responsive">
							<table id="astrometricInstrumentsTable" class="table table-striped table-bordered" >
								<thead>
									<tr>
										<th>Name</th>
										<th>Index</th>
										<th>NFields</th>
										<th>Number of Keywords</th>
										<th>Keywords</th>
										<th class="showplot">Distortion Plot</th>
									</tr>
								</thead>
								<tbody>
								</tbody>
							</table>
						</div>
					</div>
				</div> <!-- end col -->

				<!-- PHOTOINSTRUMENTS TABLE -->
				<div class="col-md-6"> 
					<div id="photoInstruments" class="panel panel-default">
						<div class="panel-heading text-left">
							<button 
								class="btn btn-info" 
								type="button" 
								data-toggle="collapse" 
								data-target="#photoInstrumentsCollapse" 
								aria-expanded="false" 
								aria-controls="photoInstrumentsCollapse">
								<strong>Photometric Instruments</strong> (PhotInstruments)
							</button>
						</div>
						<div id="photoInstrumentsCollapse" class="collapse panel-body table-responsive">
							<table id="photometricInstrumentsTable" class="table table-striped table-bordered">
								<thead>
									<tr>
										<th>Name</th>
										<th>Index</th>
										<th>NFields</th>
										<th>Output ZP</th>
										<th>Number of Keywords</th>
										<th>Keywords</th>
									</tr>
								</thead>
								<tbody>
								</tbody>
							</table>
						</div>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->



			<!-- SCAMP CONFIG -->
			<div class="row">
				<div class="col-md-12"> 
					<div id="configFile" class="panel panel-default">
						<div class="panel-heading text-left">
							<button 
								class="btn btn-info" 
								type="button" 
								data-toggle="collapse" 
								data-target="#configCollapse" 
								aria-expanded="false" 
								aria-controls="configCollapse">
								<strong>Scamp Configuration</strong>
							</button> 
							<span> 
								<a  href="javascript:void(0);" onclick="downloadScampConf();"> scamp.conf </a> 
							</span>
						</div>
						<div id="configCollapse" class="collapse panel-body table-responsive">
							<table id="configTable" class="table t table-responsive table-striped table-bordered table-sm">
								<thead>
									<tr>
										<th>Config Parameter</th>
										<th>Value</th>
									</tr>
								</thead>
								<tbody>
								</tbody>
							</table>
						</div>
					</div>
				</div> <!-- end col -->
			</div> <!-- end row -->
		</div> <!-- end main container -->



		<!-- MODAL DIALOG image viewer -->
		<div 
			class="modal fade" 
			id="imageModal" 
			tabindex="-1" 
			role="dialog" 
			aria-labelledby="imageModalLabel" 
			aria-hidden="true">
			<div class="modal-dialog modal-dialog-centered" role="document" style="width:1000px;">
				<div class="modal-content">
					<div class="modal-header">
						<h5 class="modal-title" id="imageModalLabel">New message</h5>
						<button type="button" class="close" data-dismiss="modal" aria-label="Close">
							<span aria-hidden="true">&times;</span>
						</button>
					</div>
					<div class="modal-body">
						<img class="modal-img" style="margin-left:auto;margin-right:auto;display:block;"></img>
					</div>
					<div class="modal-footer">
						<button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
					</div>
				</div>
			</div>
		</div>




		<!-- SCRIPTS INCLUDES -->
		<script
			src="https://code.jquery.com/jquery-3.3.1.min.js"
			integrity="sha256-FgpCb/KJQlLNfOu91ta32o/NMZxltwRo8QtmkMRdAu8="
			crossorigin="anonymous"></script>
		<script 
			src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" 
			integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" 
			crossorigin="anonymous"></script>
		<script 
			type="text/javascript" 
			charset="utf8" 
			src="https://cdn.datatables.net/1.10.18/js/jquery.dataTables.min.js"></script>
		<script 
			type="text/javascript" 
			charset="utf8" 
			src="https://cdn.datatables.net/1.10.18/js/dataTables.bootstrap.min.js"></script>
		<script 
			type="text/javascript"
			charset="utf-8"
			src="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js"></script>
		<script 
			type="text/javascript"
			charset="utf-8"
			src="https://cdn.plot.ly/plotly-1.41.0.min.js"></script>
		<script>
			/*! 
			 * IE10 viewport hack for Surface/desktop Windows 8 bug
			 * Copyright 2014-2015 Twitter, Inc.
			 * Licensed under MIT (https://github.com/twbs/bootstrap/blob/master/LICENSE)
			 */
			(
				function () {
					'use strict';
					if (navigator.userAgent.match(/IEMobile\/10\.0/)) {
						var msViewportStyle = document.createElement('style')
						msViewportStyle.appendChild(
							document.createTextNode(
								'@-ms-viewport{width:auto!important}'
							)
						)
						document.querySelector('head').appendChild(msViewportStyle)
					}
				} 
			)();
		</script>
	</body>
</html>

<script>
/*
 * Data contained in "scamp_data" object, is loaded at the very end 
 * of this html file is his own <script> tag.
 */

class TableUtils {

	/* from an array of object (data), return the "value" property of   
	 * object havinig "str" as "name" property. */   
	static getElemVal(str, data) {
		var value = "";
		$.each(data, function (i, elem) {
			if (elem.name == str) {
				value = elem.value;
				return;
			}
		});
		return value;
	}

	static getFlagValHelper(b, t) {
		if (b)
			return t;
		return "-";
	}

	static getElemAverageValHelper(a) {
		var value = 0.0;
		var n = 0;
		$.each(a, function(i, elem) {
			value += elem;
			n++;
		});
		return value / n;
	}

	static getElemListValHelperFixedHelper(arr, unit, fix) {
		var value = "";
		$.each(arr, function(i, elem) {
			if (fix >= 0) {
				elem = elem.toFixed(fix);
			}
			value += elem + unit + " ";
		});
		return value;
	}

	static getElemListValHelper(arr, unit) {
		return TableUtils.getElemListValHelperFixedHelper(arr, unit, -1);
	}


	static getRaValHelper(value) {
		var a = Math.floor(value[0] / 15.0);
		var b = Math.floor((value[0] * 4) % 60);
		var c = Math.floor((value[0] *240) % 60);
		return a + ":" + b + ":" + c.toFixed(2);
	}


	static getDecValHelper(value) {
		var sign = "";
		if (value[1] < 0) {
			sign = "-";
			value[1] = 0 - value[1];
		}
		var a = Math.floor(value[1]);
		var b = Math.floor((value[1] * 60) % 60);
		var c = Math.floor((value[1] * 3600) % 60);
		return sign + a + ":" + b + ":" + c.toFixed(2);
	}


	static generateImageColHelper(imageUrl) {
		if (imageUrl == null)
			return "<td></td>";

		var value = "";
		/* value += "<td><a type='button' rel='popover' data-img='"+imageUrl+"'>"; */
		value += "<td><a type='button' data-toggle='modal' data-target='#imageModal' data-imgurl='" + imageUrl.value + "'>";
		value += "<img width='100' class='img-fluid' src='" + imageUrl.value + "' />";
		value += "</a></td>";
		return value;
	}

}


aladin = null;
class AladinScampUtils {
	static init() {
		/* first initialize aladin */
		aladin = A.aladin('#aladin-lite-div', {
			survey: "P/DSS2/color", 
			fov: 60,
			showReticle: true,
			showZoomControl: true,
			showFullScreenControl: true,
			showLayersControl: true,
			showGotoControl: true,
			showShareControl: true,
			fulscreen: true
		});

		var left_max	= 1000;
		var top_max	 = -1;
		var right_max   = -1;
		var bottom_max  = 1000;
		$.each(scamp_data.Fields, function(i, field) {
			for (var i=0; i<field.Set_Polygon.value.length; i++) {
				var poly = field.Set_Polygon.value[i];
				if (poly[0][0] < left_max)
					left_max = poly[0][0];
				if (poly[0][1] > top_max)
					top_max = poly[0][1];
				if (poly[1][0] > right_max)
					right_max = poly[1][0];
				if (poly[2][1] < bottom_max)
					bottom_max = poly[2][1];
			}
		})

		var center_lng = (left_max + right_max) / 2;
		var center_lat = (top_max + bottom_max) / 2;
		aladin.gotoRaDec(center_lng, center_lat);

		var height_max = Math.abs(top_max - bottom_max);
		aladin.setFov(height_max * 4);

		AladinScampUtils.draw(-1);
	}


	static draw(higlight) {
		aladin.removeLayers();
		var lineWidth = 1;
		var color = '#bbbbbbff';
		var overlay;
		var poly;
		var drawline;

		$.each(scamp_data.Fields, function(i, field) {
			if (higlight == i)
				return;

			overlay = A.graphicOverlay({color: color, lineWidth: lineWidth});
			aladin.addOverlay(overlay);

			for (var i=0; i<field.Set_Polygon.value.length; i++) {
				poly = field.Set_Polygon.value[i];
				drawline = poly;
				drawline.push(poly[0]);
				overlay.add(A.polyline(drawline));
			}
		})

		if (higlight < 0)
			return;

		/* remain the hiligted one */
		lineWidth = 10;
		color = '#ff0000ff';
		var field = scamp_data.Fields[higlight];
		var overlay = A.graphicOverlay({color: color, lineWidth: lineWidth});
		aladin.addOverlay(overlay);
		for (var i=0; i<field.Set_Polygon.value.length; i++) {
			var poly = field.Set_Polygon.value[i];
			var drawline = poly;
			drawline.push(poly[0]);
			overlay.add(A.polyline(drawline));
		}
	}
}

class Histograms {
	static getPlotColors() {
		return {
			'#plot1': '#fce84f',
			'#plot2': '#4e9a06'
		};
	}

	static getFirstColumn() {
		return {
			'Fields': 'Observation_Date',
			'Fgroups': 'Max_Radius',
			'AstroInstruments': 'NFields',
			'PhotInstruments': 'NFields'
		};
	}


	static setDropDownMenuCallback() {
		$(".dropdown-menu li a").click(function() {
			var value   = $(this).text();
			var listId  = $(this).parent().parent().attr('id');
			var isColumn = listId.endsWith('_cols');
			var divId = '';
			if (isColumn) {
				divId = listId.substring(0, listId.length - 5)
				var tableId = $('#' + divId + '_btn1').text();
				Histograms.generate(tableId, value, '#' + divId, Histograms.getPlotColors()['#'+divId]);
			} else {
				divId = listId.substring(0, listId.length - 4)
				Histograms.generate(value, Histograms.getFirstColumn()[value], '#' + divId, Histograms.getPlotColors()['#'+divId]);
			}
		});
	}

	static generate(table, column, divId, color) {
		var fieldsCols =  `<li class="dropdown-header">Column</li>
			<li role="separator" class="divider"></li>
			<li><a href="javascript:void(0);">Observation_Date</a></li>
			<li><a href="javascript:void(0);">Exposure_Time</a></li>
			<li><a href="javascript:void(0);">Air_Mass</a></li>
			<li><a href="javascript:void(0);">Pixel_Scale</a></li>
			<li><a href="javascript:void(0);">AS_Contrast</a></li>
			<li><a href="javascript:void(0);">XY_Contrast</a></li>`;
		var groupsCols =  `<li class="dropdown-header">Column</li>
			<li role="separator" class="divider"></li>
			<li><a href="javascript:void(0);">Max_Radius</a></li>
			<li><a href="javascript:void(0);">AstromNDets_Internal</a></li>`;
		var astroCols = `<li class="dropdown-header">Column</li>
						 <li role="separator" class="divider"></li>
						 <li><a href="javascript:void(0);">NFields</a></li>`;
		var photoCols = `<li class="dropdown-header">Column</li>
						 <li role="separator" class="divider"></li>
						 <li><a href="javascript:void(0);">NFields</a></li>`;

		$(divId + '_btn1').html(table);
		$(divId + '_btn2').html(column);

		if (table == 'Fields') {
			$(divId + '_cols').html(fieldsCols);
		} else if (table == 'Fgroups') {
			$(divId + '_cols').html(groupsCols);
		} else if (table == 'PhotInstruments') {
			$(divId + '_cols').html(astroCols);
		} else if (table == 'AstroInstruments') {
			$(divId + '_cols').html(photoCols);
		}
		Histograms.setDropDownMenuCallback();

		var d3 = Plotly.d3;
		var gd = d3.select(divId).node();
		var x = [];
		$.each(scamp_data[table], function(i, row) {
			x[i] = row[column].value;
		});

		var layout = {
			showlegend:  false,
			xaxis: {title: column},
			margin: {l:50, r:50, b:80, t:40, pad:4},
			height: 300
		};

		Plotly.newPlot(gd, [{
			type: 'histogram',
			autobinx: false,
			histfunc: "density",
			x: x,
			marker: {color: color, line: {color: '#555',width:1}},
			opacity: 0.8,
		}], layout, {displayModeBar:true,showLink:false,displayLogo:false});

		return gd;
	}

}


function downloadScampConf() {
	/* generate scamp conf string */
	var text = "";
	for (var i = 0; i < scamp_data.Configuration.length; i++) {
		var configObj = scamp_data.Configuration[i];
		text += configObj.name + "\t" + configObj.value + '\n';
	}
	var element = document.createElement('a');
	element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
	element.setAttribute('download', 'scamp.conf');
	element.style.display = 'none';
	document.body.appendChild(element);
	element.click();
	document.body.removeChild(element);
}

$(document).ready(function() {

	console.log(scamp_data);
	AladinScampUtils.init();


	/* build status string */
	$('#soft').text(scamp_data.Software.Name.value +" "+scamp_data.Software.Version.value);
	$('#date').text(scamp_data.Software.Date.value);
	$('#time').text(scamp_data.Software.Time.value);
	$('#nthreads').text(scamp_data.Software.NThreads.value);
	$('#runtime').text(scamp_data.Software.Duration.value + " s");
	$('#username').text(scamp_data.Software.User.value);
	$('#rundir').text(scamp_data.Software.Path.value);



	/* show/hide match option and plots */
	var showmatch = TableUtils.getElemVal("MATCH", scamp_data.Configuration);
	var showplot  = (TableUtils.getElemVal("CHECKPLOT_DEV", scamp_data.Configuration)[0] == "PNG") ? true: false;



	/* 
	 * build fields table 
	 */
	$.each(scamp_data.Fields, function(i, field) {
		var tdClass = "<td>";
		var magLimit = 2.0;
		var fieldNum = i + 1;
		if (field.XY_Contrast.value.toFixed(1) < magLimit) {
			tdClass = '<td class="danger">';
			scamp_data.Warnings.push({
				Date: {value: scamp_data.Software.Date.value},
				Text: {value: "Catalog number " + fieldNum + " failled the crossmatch with the external catalog. (catalog file: " + field.Catalog_Name.value + ")"},
				Time: {value: scamp_data.Software.Time.value}
			});
		} 
		var table_row = "";
		table_row += '<tr>';
		table_row += tdClass +  field.Catalog_Number.value + "</td>";
		table_row += tdClass +  field.Catalog_Name.value + "</td>";
		table_row += tdClass +  field.Image_Ident.value + "</td>";
		table_row += tdClass +  field.NExtensions.value + "</td>";
		table_row += tdClass +  field.NDetect.value + "</td>";
		table_row += tdClass +  TableUtils.getFlagValHelper(field.Ext_Header.value, "H")  + TableUtils.getFlagValHelper(field.Photom_Flag, "P") +  "</td>";
		table_row += tdClass +  field.Group.value + "</td>";
		table_row += tdClass +  field.Astr_Instrum.value + "</td>";
		table_row += tdClass +  field.Phot_Instrum.value + "</td>";
		table_row += tdClass +  field.Observation_Date.value + "</td>";
		table_row += tdClass +  field.Exposure_Time.value.toFixed(3) + "</td>";
		table_row += tdClass +  field.Air_Mass.value.toFixed(2) + "</td>";
		table_row += tdClass +  TableUtils.getRaValHelper(field.Field_Coordinates.value) + "</td>";
		table_row += tdClass +  TableUtils.getDecValHelper(field.Field_Coordinates.value) + "</td>";
		table_row += tdClass +  field.Max_Radius.value.toFixed(3) + "'" + "</td>";
		table_row += tdClass +  TableUtils.getElemAverageValHelper(field.Pixel_Scale.value).toFixed(4) + "''" + "</td>";
		if (showmatch) {
			table_row += tdClass +  field.DPixel_Scale.value.toFixed(4) + "</td>";
			table_row += tdClass +  field.DPos_Angle.value + "°" + "</td>";
			table_row += tdClass +  field.AS_Contrast.value.toFixed(1) + "</td>";
			table_row += tdClass +  field.DX.value.toExponential() + "°" + "</td>";
			table_row += tdClass +  field.DY.value.toExponential() + "°" + "</td>";
			table_row += tdClass +  field.XY_Contrast.value.toFixed(1) + "</td>";
		} else {
			table_row += tdClass + "</td>";
			table_row += tdClass + "</td>";
			table_row += tdClass + "</td>";
			table_row += tdClass + "</td>";
			table_row += tdClass + "</td>";
			table_row += tdClass + "</td>";
		}
		table_row += tdClass +  field.Chi2_Internal.value.toFixed(2) + "</td>";
		table_row += tdClass +  field.Chi2_Internal_HighSN.value.toFixed(2) + "</td>";
		table_row += tdClass +  field.Chi2_Reference.value.toFixed(2) + "</td>";
		table_row += tdClass +  field.Chi2_Reference_HighSN.value.toFixed(2) + "</td>";
		table_row += tdClass +  field.ZeroPoint_Corr.value.toFixed(3) + "</td>";
		table_row += "</tr>";
		$(table_row).appendTo("#fieldsTable tbody");
	});

	/* activate datatable for sorting */
	$('#fieldsTable').DataTable({
		paging: false,
		info: false,
		searching: true
	});

	if (!showmatch) {
		$('#fieldsTable th:nth-child(17)').hide();
		$('#fieldsTable td:nth-child(17)').hide();
		$('#fieldsTable th:nth-child(18)').hide();
		$('#fieldsTable td:nth-child(18)').hide();
		$('#fieldsTable th:nth-child(19)').hide();
		$('#fieldsTable td:nth-child(19)').hide();
		$('#fieldsTable th:nth-child(20)').hide();
		$('#fieldsTable td:nth-child(20)').hide();
		$('#fieldsTable th:nth-child(21)').hide();
		$('#fieldsTable td:nth-child(21)').hide();
		$('#fieldsTable th:nth-child(22)').hide();
		$('#fieldsTable td:nth-child(22)').hide();
	}

	/* configure highlighting of footprint on aladin on click */
	$('#fieldsTable tbody tr').on('click', function(event) {
		if ($(this).hasClass('info'))  {
			$(this).removeClass('info');
			AladinScampUtils.draw(-1);
			return;
		}
		$(this).addClass('info').siblings().removeClass('info');
		AladinScampUtils.draw($(this).children('td:first').text() - 1);
	});





	/* 
	 * build fields groups table 
	 */
	$.each(scamp_data.Fgroups, function(i, group) {
		var table_row = "";
		table_row += "<tr>";
		table_row += "<td>" +  group.Name.value + "</td>";
		if (showplot) {
			table_row += TableUtils.generateImageColHelper(group.FgroupsPlot);
		} else {
			table_row += "<td></td>";
		}
		table_row += "<td>" +  group.Index.value + "</td>";
		table_row += "<td>" +  group.NFields.value + "</td>";
		table_row += "<td>" +  TableUtils.getRaValHelper(group.Field_Coordinates.value) + "</td>";
		table_row += "<td>" +  TableUtils.getDecValHelper(group.Field_Coordinates.value) + "</td>";
		table_row += "<td>" +  TableUtils.getElemAverageValHelper(group.Pixel_Scale.value).toFixed(4) + "''" + "</td>";
		table_row += "<td>" +  group.Max_Radius.value.toFixed(3) + "'" + "</td>";
		table_row += "<td>" +  group.AstRef_Catalog.value + "</td>";
		table_row += "<td>" +  group.AstRef_Band.value + "</td>";
		if (showplot) {
			table_row += TableUtils.generateImageColHelper(group.Chi2Plot);
		} else {
			table_row += "<td></td>";
		}
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.AstromSigma_Internal.value, "'' ", 4) + "</td>";
		table_row += "<td>" +  group.AstromCorr_Internal.value.toFixed(5) + "</td>";
		table_row += "<td>" +  group.AstromChi2_Internal.value.toFixed(1) + "</td>";
		table_row += "<td>" +  group.AstromNDets_Internal.value + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.AstromSigma_Internal_HighSN.value, "'' ", 4) + "</td>";
		table_row += "<td>" +  group.AstromCorr_Internal_HighSN.value.toFixed(5)  + "</td>";
		table_row += "<td>" +  group.AstromChi2_Internal_HighSN.value.toFixed(1) + "</td>";
		table_row += "<td>" +  group.AstromNDets_Internal_HighSN.value + "</td>";
		if (showplot) {
			table_row += TableUtils.generateImageColHelper(group.IntErr1DimPlot);
			table_row += TableUtils.generateImageColHelper(group.IntErr2DimPlot);
		} else {
			table_row += "<td></td>";
			table_row += "<td></td>";
		}
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.AstromOffset_Reference.value, "'' ", 4) + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.AstromSigma_Reference.value, "'' ", 3) + "</td>";
		table_row += "<td>" +  group.AstromCorr_Reference.value.toFixed(4) + "</td>";
		table_row += "<td>" +  group.AstromChi2_Reference.value.toFixed(1) + "</td>";
		table_row += "<td>" +  group.AstromNDets_Reference.value + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.AstromOffset_Reference_HighSN.value, "'' ", 4) + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.AstromSigma_Reference_HighSN.value, "'' ", 3) + "</td>";
		table_row += "<td>" +  group.AstromCorr_Reference_HighSN.value.toFixed(4) + "</td>";
		table_row += "<td>" +  group.AstromChi2_Reference_HighSN.value.toFixed(1) + "</td>";
		table_row += "<td>" +  group.AstromNDets_Reference_HighSN.value + "</td>";
		if (showplot) {
			table_row += TableUtils.generateImageColHelper(group.RefErr1DimPlot);
			table_row += TableUtils.generateImageColHelper(group.RefErr2DimPlot);
		} else {
			table_row += "<td></td>";
			table_row += "<td></td>";
		}
		table_row += "<td>" +  TableUtils.getElemListValHelper(group.PhotInstru_Name.value,", ") + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.PhotSigma_Internal.value, " ", 6) + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.PhotChi2_Internal.value, " ", 4) + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelper(group.PhotNDets_Internal.value, " ") + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.PhotSigma_Internal_HighSN.value, " ", 6) + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.PhotChi2_Internal_HighSN.value, " ", 2) + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelper(group.PhotNDets_Internal_HighSN.value, " ") + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.PhotSigma_Reference.value, " ", 6) + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.PhotChi2_Reference.value, " ", 6) + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelper(group.PhotNDets_Reference.value, " ") + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.PhotSigma_Reference_HighSN.value, " ", 6) + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelperFixedHelper(group.PhotChi2_Reference_HighSN.value, " ", 6) + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelper(group.PhotNDets_Reference_HighSN.value, " ") + "</td>";
		if (showplot) {
			table_row += TableUtils.generateImageColHelper(group.PhotErrPlot);
		} else {
			table_row += "<td></td>";
		}
		table_row += "</tr>";
		$(table_row).appendTo("#groupsTable tbody");
	});

	if (!showplot) {
		$('#groupsTable th:nth-child(2)').hide();
		$('#groupsTable td:nth-child(2)').hide();
		$('#groupsTable th:nth-child(11)').hide();
		$('#groupsTable td:nth-child(11)').hide();
		$('#groupsTable th:nth-child(20)').hide();
		$('#groupsTable td:nth-child(20)').hide();
		$('#groupsTable th:nth-child(21)').hide();
		$('#groupsTable td:nth-child(21)').hide();
		$('#groupsTable th:nth-child(32)').hide();
		$('#groupsTable td:nth-child(32)').hide();
		$('#groupsTable th:nth-child(33)').hide();
		$('#groupsTable td:nth-child(33)').hide();
		$('#groupsTable th:nth-child(47)').hide();
		$('#groupsTable td:nth-child(47)').hide();
		$('#astrometricInstrumentsTable th:nth-child(6)').hide();   
		$('#astrometricInstrumentsTable td:nth-child(6)').hide();  
	}

	$('#groupsTable').DataTable({
		paging: false,
		info: false,
		searching: false
	});




	/* 
	 * build astrometric instruments table 
	 */
	$.each(scamp_data.AstroInstruments, function(i, astroinstru) {
		var table_row = "";
		table_row += "<tr>";
		table_row += "<td>" +  astroinstru.Name.value + "</td>";
		table_row += "<td>" +  astroinstru.Index.value + "</td>";
		table_row += "<td>" +  astroinstru.NFields.value + "</td>";
		table_row += "<td>" +  astroinstru.NKeys.value + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelper(astroinstru.Keys.value, " ") + "</td>";
		if (showplot) {
			table_row += TableUtils.generateImageColHelper(astroinstru.DistPlot);
		} else {
			table_row += "<td></td>";
		}
		table_row += "</tr>";
		$(table_row).appendTo("#astrometricInstrumentsTable tbody");
	});

	$('#astrometricInstrumentsTable').DataTable({
		paging: false,
		info: false,
		searching: false
	});




	/* 
	 * build photometric instruments table 
	 */
	$.each(scamp_data.PhotInstruments, function(i, photoinstru) {
		var table_row = "";
		table_row += "<tr>";
		table_row += "<td>" +  photoinstru.Name.value + "</td>";
		table_row += "<td>" +  photoinstru.Index.value + "</td>";
		table_row += "<td>" +  photoinstru.NFields.value + "</td>";
		table_row += "<td>" +  photoinstru.MagZeroPoint_Output.value + "</td>";
		table_row += "<td>" +  photoinstru.NKeys.value + "</td>";
		table_row += "<td>" +  TableUtils.getElemListValHelper(photoinstru.Keys.value, " ") + "</td>";
		table_row += "</tr>";
		$(table_row).appendTo("#photometricInstrumentsTable tbody");
	});

	$('#photometricInstrumentsTable').DataTable({
		paging: false,
		info: false,
		searching: false
	});



	/* 
	 * build configuration table 
	 */
	$.each(scamp_data.Configuration, function(i, config) {
		var table_row = "";
		table_row += "<tr>";
		table_row += "<td>" +  config.name + "</td>";

		var value = "";
		if (config.datatype.includes("array")) {
			for (var i = 0; i < config.value.length; i++) {
				value += config.value[i] + ", ";
			}
		} else {
			value = config.value;
		}
		table_row += "<td>" +  value + "</td>";
		table_row += "</tr>";
		$(table_row).appendTo("#configTable tbody");
	});

	$('#configTable').DataTable({
		scrollY: "400px",
		scrollCollapse: true,
		paging: false,
		info: false,
		searching: false
	});




	/* 
	 * build warnings table 
	 */
	$.each(scamp_data.Warnings, function(i, warn) {
		var table_row = "";
		table_row += "<tr>";
		table_row += "<td>" +  warn.Date.value + "</td>";
		table_row += "<td>" +  warn.Time.value + "</td>";
		table_row += "<td>" +  warn.Text.value + "</td>";
		$(table_row).appendTo("#warningsTable tbody");
	});

	$('#warningsTable').DataTable({
		paging: false,
		info: false,
		searching: false
	});

	if (scamp_data.Warnings.length == 0) {
		$("#warningDiv").hide();
	} else {
		var nwarn = $('#warningsTable tbody tr').length;
		$('#warningsShort').html(' (' + nwarn + ')');
	}


	$('#imageModal').on('show.bs.modal', function (event) {
		var button = $(event.relatedTarget); // Button that triggered the modal
		var img = button.data('imgurl'); // Extract info from data-* attributes
		var modal = $(this);
		modal.find('.modal-title').text("Image: " + img);
		imgtag = modal.find('.modal-img');
		imgtag.attr("src", img);
		imgtag.attr("title", img);
		imgtag.attr("alt", img);
	})

	$('a[rel=popover]').popover({
		animation: true, 
		container: 'body', 
		html: true, 
		placement: 'bottom', 
		content: function() {
			return "<img src='"+$(this).data('img') + "' />";
		}
	});


	var gd1 = Histograms.generate('Fields', 'Observation_Date', '#plot1', Histograms.getPlotColors()['#plot1']);
	var gd2 = Histograms.generate('Fields', 'XY_Contrast', '#plot2', Histograms.getPlotColors()['#plot2']);

	window.onresize = function() {
		Plotly.Plots.resize(gd1);
		Plotly.Plots.resize(gd2);
	}


	$(".table-dropdown").each(function(index) {
		$(this).append(`
			<li class="dropdown-header">Tables</li>
			<li role="separator" class="divider"></li>
			<li><a href="javascript:void(0);">Fields</a></li>
			<li><a href="javascript:void(0);">Fgroups</a></li>
			<li><a href="javascript:void(0);">PhotInstruments</a></li>
			<li><a href="javascript:void(0);">AstroInstruments</a></li>
		`);
	});

	Histograms.setDropDownMenuCallback();

	$("#configCollapse").on('shown.bs.collapse', function() {
		$('#configTable').DataTable().columns.adjust().draw();
	});


	document.onkeydown = checkKey;

	function checkKey(e) {
		e = e || window.event;
		console.log('has focus: ' + $(':focus').attr('id'));
		if ($(':focus').attr('id') == 'fieldsTable') {
			if (e.keyCode == '38') {
				// up arrow
				e.preventDefault();
				$('#fieldsTable tbody tr.info').prev().click();
			} else if (e.keyCode == '40') {
				// down arrow
				e.preventDefault();
				$('#fieldsTable tbody tr.info').next().click();
			}
		}
	}
});
</script>
<script>
	var scamp_data = $.parseJSON('{ "Software": { "Name": { "name": "Name", "datatype": "string", "ucd": "meta.title;meta.software", "value": "SCAMP" }, "Version": { "name": "Version", "datatype": "string", "ucd": "meta.version;meta.software", "value": "2.7.3" }, "Url": { "name": "Url", "datatype": "string", "ucd": "meta.ref.url;meta.software", "value": "http:\/\/astromatic.net\/software\/scamp" }, "Author": { "name": "Author", "datatype": "string", "ucd": "meta.bib.author;meta.software", "value": "Emmanuel Bertin" }, "Ref": { "name": "Ref", "datatype": "string", "ucd": "meta.bib.bibcode;meta.software", "value": "2006ASPC..351..112B" }, "NThreads": { "name": "NThreads", "datatype": "int", "ucd": "meta.number;meta.software", "value": 4 }, "Date": { "name": "Date", "datatype": "string", "ucd": "time.end;meta.software", "value": "2018-10-08" }, "Time": { "name": "Time", "datatype": "string", "ucd": "time.end;meta.software", "value": "09:28:44" }, "Duration": { "name": "Duration", "datatype": "float", "ucd": "time.duration;meta.software", "value": 56.000000 }, "User": { "name": "User", "datatype": "string", "ucd": "meta.curation", "value": "serre" }, "Path": { "name": "Path", "datatype": "string", "ucd": "meta.dataset", "value": "\/home\/serre\/src\/scamp" } }, "Fields": [ { "Catalog_Number": { "name": "Catalog_Number", "datatype": "int", "ucd": "meta.record;meta.table;meta.file", "value": 1 }, "Catalog_Name": { "name": "Catalog_Name", "datatype": "string", "ucd": "meta.id;meta.table;meta.file", "value": "744331p.cat" }, "Image_Ident": { "name": "Image_Ident", "datatype": "string", "ucd": "meta.id;obs.field", "value": "13h" }, "NExtensions": { "name": "NExtensions", "datatype": "int", "ucd": "meta.record", "value": 36 }, "NAxis": { "name": "NAxis", "datatype": "int", "ucd": "pos.wcs.naxis", "value": 2 }, "Lng_Axis": { "name": "Lng_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.ra", "value": 0 }, "Lat_Axis": { "name": "Lat_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.dec", "value": 1 }, "Ext_Header": { "name": "Ext_Header", "datatype": "boolean", "ucd": "meta.code", "value": false }, "NDetect": { "name": "NDetect", "datatype": "int", "ucd": "meta.number;src", "value": 30934 }, "Group": { "name": "Group", "datatype": "int", "ucd": "meta.id.parent;meta.dataset", "value": 1 }, "Astr_Instrum": { "name": "Astr_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "A1" }, "Phot_Instrum": { "name": "Phot_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "P1" }, "Photom_Flag": { "name": "Photom_Flag", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Photom_Link": { "name": "Photom_Link", "datatype": "boolean", "ucd": "meta.code;phot", "value": true }, "Observation_Date": { "name": "Observation_Date", "datatype": "double", "ucd": "time.epoch;obs.field", "unit": "yr", "value": 2004.353010 }, "Exposure_Time": { "name": "Exposure_Time", "datatype": "float", "ucd": "time.duration;obs.exposure", "value": 512.204000 }, "Air_Mass": { "name": "Air_Mass", "datatype": "float", "ucd": "obs.airMass", "value": 1.132000 }, "Field_Coordinates": { "name": "Field_Coordinates", "datatype": "double array", "ucd": "pos.eq.ra;pos.eq.dec;obs.field", "unit": "%s", "value": [ 203.650675, 37.913065 ] }, "Pixel_Scale": { "name": "Pixel_Scale", "datatype": "float array", "ucd": "instr.scale;instr.pixel;stat.mean", "unit": "%s", "value": [ 0.186001, 0.186001 ] }, "Max_Radius": { "name": "Max_Radius", "datatype": "float", "ucd": "phys.size.radius", "unit": "%s", "value": 42.650690 }, "ZeroPoint_Corr": { "name": "ZeroPoint_Corr", "datatype": "float", "ucd": "phot.mag;phot.calib;arith.zp", "unit": "mag", "value": 0.000054 }, "DPixel_Scale": { "name": "DPixel_Scale", "datatype": "float", "ucd": "instr.scale;instr.pixel;arith.ratio", "value": 0.999902 }, "DPos_Angle": { "name": "DPos_Angle", "datatype": "float", "ucd": "pos.posAng;obs.image;arith.diff", "unit": "deg", "value": 0.002772 }, "AS_Contrast": { "name": "AS_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 26.500492 }, "DX": { "name": "DX", "datatype": "float", "ucd": "pos.eq.ra;arith.diff", "unit": "deg", "value": 0.000198 }, "DY": { "name": "DY", "datatype": "float", "ucd": "pos.eq.dec;arith.diff", "unit": "deg", "value": 0.000018 }, "XY_Contrast": { "name": "XY_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 32.271492 }, "Shear": { "name": "Shear", "datatype": "float", "ucd": "phys.size.axisRatio;obs.image", "value": 0.000032 }, "Shear_PosAngle": { "name": "Shear_PosAngle", "datatype": "float", "ucd": "pos.posAng;obs.image", "unit": "deg", "value": 38.105766 }, "Chi2_Internal": { "name": "Chi2_Internal", "datatype": "float", "ucd": "stat.fit.chi2", "value": 92.241987 }, "NDeg_Internal": { "name": "NDeg_Internal", "datatype": "int", "ucd": "stat.fit.dof", "value": 156344 }, "Chi2_Internal_HighSN": { "name": "Chi2_Internal_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 450.014824 }, "NDeg_Internal_HighSN": { "name": "NDeg_Internal_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 16859 }, "AstromOffset_Reference": { "name": "AstromOffset_Reference", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ 0.000268, -0.009764 ] }, "Astrom_Reference": { "name": "AstromSigma_Reference", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.315247, 0.360762 ] }, "AstromCorr_Reference": { "name": "AstromCorr_Reference", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.077123 }, "Chi2_Reference": { "name": "Chi2_Reference", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.525224 }, "NDeg_Reference": { "name": "NDeg_Reference", "datatype": "int", "ucd": "stat.fit.dof", "value": 759 }, "AstromOffset_Reference_HighSN": { "name": "AstromOffset_Reference_HighSN", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ 0.000769, -0.011128 ] }, "AstromSigma_Reference_HighSN": { "name": "AstromSigma_Reference_HighSN", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.314002, 0.357438 ] }, "AstromCorr_Reference_HighSN": { "name": "AstromCorr_Reference_HighSN", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.070054 }, "Chi2_Reference_HighSN": { "name": "Chi2_Reference_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.447135 }, "NDeg_Reference_HighSN": { "name": "NDeg_Reference_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 754 }, "Set_Polygon": { "name": "Set_Polygon", "datatype": "float array", "ucd": "qsldfjklqksj", "unit": "%s", "value": [ [ [ 204.131761, 38.415582 ], [ 204.269405, 38.415079 ], [ 204.272956, 38.177062 ], [ 204.135347, 38.177283 ] ], [ [ 203.993046, 38.415540 ], [ 204.131579, 38.415393 ], [ 204.134817, 38.177462 ], [ 203.996870, 38.177068 ] ], [ [ 203.853795, 38.415631 ], [ 203.992818, 38.415665 ], [ 203.996287, 38.176990 ], [ 203.856928, 38.176533 ] ], [ [ 203.714166, 38.415107 ], [ 203.853377, 38.415689 ], [ 203.856986, 38.176355 ], [ 203.717460, 38.175257 ] ], [ [ 203.574709, 38.413863 ], [ 203.714014, 38.414990 ], [ 203.717011, 38.175476 ], [ 203.577525, 38.174752 ] ], [ [ 203.434790, 38.412397 ], [ 203.574311, 38.413802 ], [ 203.577284, 38.174728 ], [ 203.437985, 38.173012 ] ], [ [ 203.295957, 38.410038 ], [ 203.434581, 38.412596 ], [ 203.437248, 38.173275 ], [ 203.298501, 38.171509 ] ], [ [ 203.157465, 38.407564 ], [ 203.295101, 38.409959 ], [ 203.298452, 38.171196 ], [ 203.159579, 38.169567 ] ], [ [ 203.019810, 38.404501 ], [ 203.157457, 38.407342 ], [ 203.159427, 38.169540 ], [ 203.021468, 38.167255 ] ], [ [ 204.135568, 38.157182 ], [ 204.272886, 38.156766 ], [ 204.274978, 37.916718 ], [ 204.137807, 37.916730 ] ], [ [ 203.996734, 38.156845 ], [ 204.135224, 38.157075 ], [ 204.137360, 37.916572 ], [ 203.999205, 37.916038 ] ], [ [ 203.857275, 38.156271 ], [ 203.996465, 38.156792 ], [ 203.998948, 37.915937 ], [ 203.860033, 37.915042 ] ], [ [ 203.717572, 38.155496 ], [ 203.857163, 38.156192 ], [ 203.859800, 37.915137 ], [ 203.720522, 37.914291 ] ], [ [ 203.577868, 38.154588 ], [ 203.717341, 38.155519 ], [ 203.720337, 37.914355 ], [ 203.580985, 37.913198 ] ], [ [ 203.438043, 38.152923 ], [ 203.577574, 38.154406 ], [ 203.580676, 37.913029 ], [ 203.441613, 37.911826 ] ], [ [ 203.298822, 38.151413 ], [ 203.437735, 38.153080 ], [ 203.441294, 37.911960 ], [ 203.302567, 37.910602 ] ], [ [ 203.159900, 38.149451 ], [ 203.298849, 38.151351 ], [ 203.302372, 37.910583 ], [ 203.163965, 37.908970 ] ], [ [ 203.022209, 38.147028 ], [ 203.159714, 38.149231 ], [ 203.163705, 37.908710 ], [ 203.026185, 37.907087 ] ], [ [ 204.275143, 37.678200 ], [ 204.138194, 37.676527 ], [ 204.137915, 37.917488 ], [ 204.274907, 37.917427 ] ], [ [ 204.138128, 37.676801 ], [ 204.000527, 37.675962 ], [ 203.999315, 37.916819 ], [ 204.137499, 37.917460 ] ], [ [ 204.000239, 37.675958 ], [ 203.862258, 37.674873 ], [ 203.860034, 37.915901 ], [ 203.998846, 37.916846 ] ], [ [ 203.861643, 37.674743 ], [ 203.723273, 37.673972 ], [ 203.720491, 37.914990 ], [ 203.859743, 37.915879 ] ], [ [ 203.723061, 37.673890 ], [ 203.584390, 37.672708 ], [ 203.580974, 37.913963 ], [ 203.720340, 37.915063 ] ], [ [ 203.584113, 37.672815 ], [ 203.445700, 37.671711 ], [ 203.441305, 37.912719 ], [ 203.580314, 37.914069 ] ], [ [ 203.445463, 37.671790 ], [ 203.307246, 37.670525 ], [ 203.302510, 37.911392 ], [ 203.441264, 37.912765 ] ], [ [ 203.306927, 37.670533 ], [ 203.169609, 37.669328 ], [ 203.163822, 37.909673 ], [ 203.301997, 37.911366 ] ], [ [ 203.169195, 37.669381 ], [ 203.032255, 37.668211 ], [ 203.026358, 37.907859 ], [ 203.163849, 37.909663 ] ], [ [ 204.273820, 37.419981 ], [ 204.137786, 37.418687 ], [ 204.138381, 37.656441 ], [ 204.274934, 37.658032 ] ], [ [ 204.137542, 37.419027 ], [ 204.001200, 37.417210 ], [ 204.000635, 37.655893 ], [ 204.138303, 37.656832 ] ], [ [ 204.000917, 37.417319 ], [ 203.863694, 37.415456 ], [ 203.862284, 37.654750 ], [ 204.000146, 37.655897 ] ], [ [ 203.863753, 37.415560 ], [ 203.725892, 37.414311 ], [ 203.723615, 37.653677 ], [ 203.861939, 37.654564 ] ], [ [ 203.725476, 37.414648 ], [ 203.588393, 37.413481 ], [ 203.584574, 37.652620 ], [ 203.723209, 37.653747 ] ], [ [ 203.588092, 37.413218 ], [ 203.450235, 37.412320 ], [ 203.445677, 37.651197 ], [ 203.584430, 37.652527 ] ], [ [ 203.450170, 37.412597 ], [ 203.313207, 37.411783 ], [ 203.307744, 37.650412 ], [ 203.445605, 37.651732 ] ], [ [ 203.312737, 37.412119 ], [ 203.176068, 37.410822 ], [ 203.169799, 37.649241 ], [ 203.307213, 37.650278 ] ], [ [ 203.175746, 37.410824 ], [ 203.039893, 37.410187 ], [ 203.032700, 37.648071 ], [ 203.169547, 37.649063 ] ] ] } }, { "Catalog_Number": { "name": "Catalog_Number", "datatype": "int", "ucd": "meta.record;meta.table;meta.file", "value": 2 }, "Catalog_Name": { "name": "Catalog_Name", "datatype": "string", "ucd": "meta.id;meta.table;meta.file", "value": "744332p.cat" }, "Image_Ident": { "name": "Image_Ident", "datatype": "string", "ucd": "meta.id;obs.field", "value": "13h" }, "NExtensions": { "name": "NExtensions", "datatype": "int", "ucd": "meta.record", "value": 36 }, "NAxis": { "name": "NAxis", "datatype": "int", "ucd": "pos.wcs.naxis", "value": 2 }, "Lng_Axis": { "name": "Lng_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.ra", "value": 0 }, "Lat_Axis": { "name": "Lat_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.dec", "value": 1 }, "Ext_Header": { "name": "Ext_Header", "datatype": "boolean", "ucd": "meta.code", "value": false }, "NDetect": { "name": "NDetect", "datatype": "int", "ucd": "meta.number;src", "value": 30483 }, "Group": { "name": "Group", "datatype": "int", "ucd": "meta.id.parent;meta.dataset", "value": 1 }, "Astr_Instrum": { "name": "Astr_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "A1" }, "Phot_Instrum": { "name": "Phot_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "P1" }, "Photom_Flag": { "name": "Photom_Flag", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Photom_Link": { "name": "Photom_Link", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Observation_Date": { "name": "Observation_Date", "datatype": "double", "ucd": "time.epoch;obs.field", "unit": "yr", "value": 2004.353028 }, "Exposure_Time": { "name": "Exposure_Time", "datatype": "float", "ucd": "time.duration;obs.exposure", "value": 512.221000 }, "Air_Mass": { "name": "Air_Mass", "datatype": "float", "ucd": "obs.airMass", "value": 1.149000 }, "Field_Coordinates": { "name": "Field_Coordinates", "datatype": "double array", "ucd": "pos.eq.ra;pos.eq.dec;obs.field", "unit": "%s", "value": [ 203.642916, 37.904651 ] }, "Pixel_Scale": { "name": "Pixel_Scale", "datatype": "float array", "ucd": "instr.scale;instr.pixel;stat.mean", "unit": "%s", "value": [ 0.185975, 0.185975 ] }, "Max_Radius": { "name": "Max_Radius", "datatype": "float", "ucd": "phys.size.radius", "unit": "%s", "value": 42.646130 }, "ZeroPoint_Corr": { "name": "ZeroPoint_Corr", "datatype": "float", "ucd": "phot.mag;phot.calib;arith.zp", "unit": "mag", "value": 0.001583 }, "DPixel_Scale": { "name": "DPixel_Scale", "datatype": "float", "ucd": "instr.scale;instr.pixel;arith.ratio", "value": 0.999904 }, "DPos_Angle": { "name": "DPos_Angle", "datatype": "float", "ucd": "pos.posAng;obs.image;arith.diff", "unit": "deg", "value": 0.003023 }, "AS_Contrast": { "name": "AS_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 26.756355 }, "DX": { "name": "DX", "datatype": "float", "ucd": "pos.eq.ra;arith.diff", "unit": "deg", "value": 0.000206 }, "DY": { "name": "DY", "datatype": "float", "ucd": "pos.eq.dec;arith.diff", "unit": "deg", "value": 0.000020 }, "XY_Contrast": { "name": "XY_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 30.711672 }, "Shear": { "name": "Shear", "datatype": "float", "ucd": "phys.size.axisRatio;obs.image", "value": 0.000028 }, "Shear_PosAngle": { "name": "Shear_PosAngle", "datatype": "float", "ucd": "pos.posAng;obs.image", "unit": "deg", "value": 56.534084 }, "Chi2_Internal": { "name": "Chi2_Internal", "datatype": "float", "ucd": "stat.fit.chi2", "value": 98.173945 }, "NDeg_Internal": { "name": "NDeg_Internal", "datatype": "int", "ucd": "stat.fit.dof", "value": 153904 }, "Chi2_Internal_HighSN": { "name": "Chi2_Internal_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 479.482129 }, "NDeg_Internal_HighSN": { "name": "NDeg_Internal_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 16914 }, "AstromOffset_Reference": { "name": "AstromOffset_Reference", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.005346, -0.009579 ] }, "Astrom_Reference": { "name": "AstromSigma_Reference", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.309135, 0.347105 ] }, "AstromCorr_Reference": { "name": "AstromCorr_Reference", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.043200 }, "Chi2_Reference": { "name": "Chi2_Reference", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.341890 }, "NDeg_Reference": { "name": "NDeg_Reference", "datatype": "int", "ucd": "stat.fit.dof", "value": 783 }, "AstromOffset_Reference_HighSN": { "name": "AstromOffset_Reference_HighSN", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.004271, -0.009282 ] }, "AstromSigma_Reference_HighSN": { "name": "AstromSigma_Reference_HighSN", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.309101, 0.345320 ] }, "AstromCorr_Reference_HighSN": { "name": "AstromCorr_Reference_HighSN", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.042490 }, "Chi2_Reference_HighSN": { "name": "Chi2_Reference_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.315646 }, "NDeg_Reference_HighSN": { "name": "NDeg_Reference_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 780 }, "Set_Polygon": { "name": "Set_Polygon", "datatype": "float array", "ucd": "qsldfjklqksj", "unit": "%s", "value": [ [ [ 204.123896, 38.407071 ], [ 204.261472, 38.406620 ], [ 204.265164, 38.168727 ], [ 204.127623, 38.168895 ] ], [ [ 203.985099, 38.407222 ], [ 204.123753, 38.407041 ], [ 204.127166, 38.169071 ], [ 203.989100, 38.168708 ] ], [ [ 203.846000, 38.407225 ], [ 203.984969, 38.407234 ], [ 203.988491, 38.168581 ], [ 203.849189, 38.168150 ] ], [ [ 203.706425, 38.406621 ], [ 203.845601, 38.407153 ], [ 203.849219, 38.167951 ], [ 203.709728, 38.166905 ] ], [ [ 203.566830, 38.405523 ], [ 203.706290, 38.406660 ], [ 203.709420, 38.167041 ], [ 203.569778, 38.166308 ] ], [ [ 203.427093, 38.403987 ], [ 203.566605, 38.405312 ], [ 203.569536, 38.166269 ], [ 203.430246, 38.164635 ] ], [ [ 203.288252, 38.401452 ], [ 203.426933, 38.403958 ], [ 203.429573, 38.164931 ], [ 203.290770, 38.163219 ] ], [ [ 203.149823, 38.399167 ], [ 203.287477, 38.401609 ], [ 203.290718, 38.162753 ], [ 203.151828, 38.161076 ] ], [ [ 203.011924, 38.396266 ], [ 203.149483, 38.399152 ], [ 203.151896, 38.161054 ], [ 203.014016, 38.158724 ] ], [ [ 204.127813, 38.148719 ], [ 204.265228, 38.148384 ], [ 204.267177, 37.908352 ], [ 204.129904, 37.908278 ] ], [ [ 203.988997, 38.148418 ], [ 204.127526, 38.148695 ], [ 204.129545, 37.908157 ], [ 203.991347, 37.907577 ] ], [ [ 203.849533, 38.147841 ], [ 203.988683, 38.148291 ], [ 203.991150, 37.907533 ], [ 203.852277, 37.906710 ] ], [ [ 203.709847, 38.147047 ], [ 203.849359, 38.147752 ], [ 203.852007, 37.906755 ], [ 203.712808, 37.905901 ] ], [ [ 203.570104, 38.146186 ], [ 203.709590, 38.147099 ], [ 203.712609, 37.905926 ], [ 203.573244, 37.904787 ] ], [ [ 203.430348, 38.144469 ], [ 203.569836, 38.146004 ], [ 203.572913, 37.904676 ], [ 203.433893, 37.903421 ] ], [ [ 203.291145, 38.142999 ], [ 203.430047, 38.144650 ], [ 203.433563, 37.903542 ], [ 203.294847, 37.902200 ] ], [ [ 203.152203, 38.141069 ], [ 203.290920, 38.142948 ], [ 203.294611, 37.902172 ], [ 203.156431, 37.900581 ] ], [ [ 203.014544, 38.138591 ], [ 203.151991, 38.140906 ], [ 203.155925, 37.900390 ], [ 203.018464, 37.898653 ] ], [ [ 204.267428, 37.669860 ], [ 204.130278, 37.668066 ], [ 204.129967, 37.908965 ], [ 204.267154, 37.909018 ] ], [ [ 204.130290, 37.668364 ], [ 203.992767, 37.667519 ], [ 203.991548, 37.908410 ], [ 204.129653, 37.909057 ] ], [ [ 203.992359, 37.667549 ], [ 203.854374, 37.666509 ], [ 203.852351, 37.907484 ], [ 203.991165, 37.908382 ] ], [ [ 203.853842, 37.666330 ], [ 203.715500, 37.665543 ], [ 203.712730, 37.906577 ], [ 203.851954, 37.907482 ] ], [ [ 203.715282, 37.665475 ], [ 203.576653, 37.664264 ], [ 203.573245, 37.905542 ], [ 203.712570, 37.906671 ] ], [ [ 203.576377, 37.664402 ], [ 203.438048, 37.663366 ], [ 203.433598, 37.904334 ], [ 203.572522, 37.905614 ] ], [ [ 203.437695, 37.663399 ], [ 203.299501, 37.662095 ], [ 203.294833, 37.902942 ], [ 203.433565, 37.904355 ] ], [ [ 203.299254, 37.662106 ], [ 203.161833, 37.660914 ], [ 203.156133, 37.901239 ], [ 203.294413, 37.902919 ] ], [ [ 203.161515, 37.660953 ], [ 203.024615, 37.659877 ], [ 203.018677, 37.899507 ], [ 203.156127, 37.901215 ] ], [ [ 204.266156, 37.411390 ], [ 204.130070, 37.410062 ], [ 204.130478, 37.648086 ], [ 204.267087, 37.649708 ] ], [ [ 204.129736, 37.410620 ], [ 203.993415, 37.408829 ], [ 203.992791, 37.647470 ], [ 204.130438, 37.648380 ] ], [ [ 203.993157, 37.408696 ], [ 203.855890, 37.407019 ], [ 203.854499, 37.646563 ], [ 203.992406, 37.647518 ] ], [ [ 203.856056, 37.407237 ], [ 203.718070, 37.405860 ], [ 203.715726, 37.645154 ], [ 203.854178, 37.646169 ] ], [ [ 203.717774, 37.406195 ], [ 203.580608, 37.405067 ], [ 203.576737, 37.644214 ], [ 203.715458, 37.645303 ] ], [ [ 203.580192, 37.404910 ], [ 203.442456, 37.404116 ], [ 203.438005, 37.642784 ], [ 203.576633, 37.644008 ] ], [ [ 203.442532, 37.404210 ], [ 203.305450, 37.403412 ], [ 203.299983, 37.641943 ], [ 203.437962, 37.643250 ] ], [ [ 203.304994, 37.403610 ], [ 203.168347, 37.402253 ], [ 203.162123, 37.640858 ], [ 203.299516, 37.641951 ] ], [ [ 203.168018, 37.402375 ], [ 203.032020, 37.401596 ], [ 203.024905, 37.639746 ], [ 203.161894, 37.640870 ] ] ] } }, { "Catalog_Number": { "name": "Catalog_Number", "datatype": "int", "ucd": "meta.record;meta.table;meta.file", "value": 3 }, "Catalog_Name": { "name": "Catalog_Name", "datatype": "string", "ucd": "meta.id;meta.table;meta.file", "value": "744333p.cat" }, "Image_Ident": { "name": "Image_Ident", "datatype": "string", "ucd": "meta.id;obs.field", "value": "13h" }, "NExtensions": { "name": "NExtensions", "datatype": "int", "ucd": "meta.record", "value": 36 }, "NAxis": { "name": "NAxis", "datatype": "int", "ucd": "pos.wcs.naxis", "value": 2 }, "Lng_Axis": { "name": "Lng_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.ra", "value": 0 }, "Lat_Axis": { "name": "Lat_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.dec", "value": 1 }, "Ext_Header": { "name": "Ext_Header", "datatype": "boolean", "ucd": "meta.code", "value": false }, "NDetect": { "name": "NDetect", "datatype": "int", "ucd": "meta.number;src", "value": 30592 }, "Group": { "name": "Group", "datatype": "int", "ucd": "meta.id.parent;meta.dataset", "value": 1 }, "Astr_Instrum": { "name": "Astr_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "A1" }, "Phot_Instrum": { "name": "Phot_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "P1" }, "Photom_Flag": { "name": "Photom_Flag", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Photom_Link": { "name": "Photom_Link", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Observation_Date": { "name": "Observation_Date", "datatype": "double", "ucd": "time.epoch;obs.field", "unit": "yr", "value": 2004.353046 }, "Exposure_Time": { "name": "Exposure_Time", "datatype": "float", "ucd": "time.duration;obs.exposure", "value": 512.224000 }, "Air_Mass": { "name": "Air_Mass", "datatype": "float", "ucd": "obs.airMass", "value": 1.168000 }, "Field_Coordinates": { "name": "Field_Coordinates", "datatype": "double array", "ucd": "pos.eq.ra;pos.eq.dec;obs.field", "unit": "%s", "value": [ 203.653761, 37.887970 ] }, "Pixel_Scale": { "name": "Pixel_Scale", "datatype": "float array", "ucd": "instr.scale;instr.pixel;stat.mean", "unit": "%s", "value": [ 0.185948, 0.185948 ] }, "Max_Radius": { "name": "Max_Radius", "datatype": "float", "ucd": "phys.size.radius", "unit": "%s", "value": 42.658742 }, "ZeroPoint_Corr": { "name": "ZeroPoint_Corr", "datatype": "float", "ucd": "phot.mag;phot.calib;arith.zp", "unit": "mag", "value": 0.001380 }, "DPixel_Scale": { "name": "DPixel_Scale", "datatype": "float", "ucd": "instr.scale;instr.pixel;arith.ratio", "value": 0.999923 }, "DPos_Angle": { "name": "DPos_Angle", "datatype": "float", "ucd": "pos.posAng;obs.image;arith.diff", "unit": "deg", "value": 0.003499 }, "AS_Contrast": { "name": "AS_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 25.249022 }, "DX": { "name": "DX", "datatype": "float", "ucd": "pos.eq.ra;arith.diff", "unit": "deg", "value": 0.000200 }, "DY": { "name": "DY", "datatype": "float", "ucd": "pos.eq.dec;arith.diff", "unit": "deg", "value": 0.000024 }, "XY_Contrast": { "name": "XY_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 31.937853 }, "Shear": { "name": "Shear", "datatype": "float", "ucd": "phys.size.axisRatio;obs.image", "value": 0.000027 }, "Shear_PosAngle": { "name": "Shear_PosAngle", "datatype": "float", "ucd": "pos.posAng;obs.image", "unit": "deg", "value": 39.448219 }, "Chi2_Internal": { "name": "Chi2_Internal", "datatype": "float", "ucd": "stat.fit.chi2", "value": 89.507452 }, "NDeg_Internal": { "name": "NDeg_Internal", "datatype": "int", "ucd": "stat.fit.dof", "value": 152593 }, "Chi2_Internal_HighSN": { "name": "Chi2_Internal_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 432.488930 }, "NDeg_Internal_HighSN": { "name": "NDeg_Internal_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 16581 }, "AstromOffset_Reference": { "name": "AstromOffset_Reference", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.013583, -0.013115 ] }, "Astrom_Reference": { "name": "AstromSigma_Reference", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.295628, 0.343247 ] }, "AstromCorr_Reference": { "name": "AstromCorr_Reference", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.068280 }, "Chi2_Reference": { "name": "Chi2_Reference", "datatype": "float", "ucd": "stat.fit.chi2", "value": 2.857869 }, "NDeg_Reference": { "name": "NDeg_Reference", "datatype": "int", "ucd": "stat.fit.dof", "value": 769 }, "AstromOffset_Reference_HighSN": { "name": "AstromOffset_Reference_HighSN", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.013372, -0.012726 ] }, "AstromSigma_Reference_HighSN": { "name": "AstromSigma_Reference_HighSN", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.295352, 0.341434 ] }, "AstromCorr_Reference_HighSN": { "name": "AstromCorr_Reference_HighSN", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.070276 }, "Chi2_Reference_HighSN": { "name": "Chi2_Reference_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 2.831429 }, "NDeg_Reference_HighSN": { "name": "NDeg_Reference_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 766 }, "Set_Polygon": { "name": "Set_Polygon", "datatype": "float array", "ucd": "qsldfjklqksj", "unit": "%s", "value": [ [ [ 204.134751, 38.390634 ], [ 204.272284, 38.389984 ], [ 204.275803, 38.151911 ], [ 204.138309, 38.152281 ] ], [ [ 203.995927, 38.390523 ], [ 204.134477, 38.390392 ], [ 204.137860, 38.152448 ], [ 203.999898, 38.152036 ] ], [ [ 203.856860, 38.390518 ], [ 203.995795, 38.390467 ], [ 203.999266, 38.151914 ], [ 203.859998, 38.151544 ] ], [ [ 203.717299, 38.389981 ], [ 203.856547, 38.390487 ], [ 203.860057, 38.151240 ], [ 203.720497, 38.150218 ] ], [ [ 203.577763, 38.388679 ], [ 203.717016, 38.389741 ], [ 203.720158, 38.150383 ], [ 203.580724, 38.149726 ] ], [ [ 203.437975, 38.387205 ], [ 203.577394, 38.388573 ], [ 203.580336, 38.149696 ], [ 203.441138, 38.148019 ] ], [ [ 203.299087, 38.384920 ], [ 203.437768, 38.387390 ], [ 203.440482, 38.148165 ], [ 203.301677, 38.146490 ] ], [ [ 203.160846, 38.382468 ], [ 203.298427, 38.384830 ], [ 203.301574, 38.146096 ], [ 203.162760, 38.144494 ] ], [ [ 203.022696, 38.379687 ], [ 203.160507, 38.382481 ], [ 203.163076, 38.144271 ], [ 203.024935, 38.142032 ] ], [ [ 204.138535, 38.132026 ], [ 204.275835, 38.131604 ], [ 204.277876, 37.891612 ], [ 204.140722, 37.891628 ] ], [ [ 203.999710, 38.131694 ], [ 204.138148, 38.131985 ], [ 204.140284, 37.891487 ], [ 204.002178, 37.890891 ] ], [ [ 203.860354, 38.131134 ], [ 203.999421, 38.131682 ], [ 204.001886, 37.890869 ], [ 203.863095, 37.889946 ] ], [ [ 203.720650, 38.130418 ], [ 203.860157, 38.131184 ], [ 203.862870, 37.889982 ], [ 203.723675, 37.889066 ] ], [ [ 203.580950, 38.129494 ], [ 203.720419, 38.130382 ], [ 203.723444, 37.889247 ], [ 203.584095, 37.888132 ] ], [ [ 203.441259, 38.127847 ], [ 203.580692, 38.129290 ], [ 203.583761, 37.887930 ], [ 203.444795, 37.886766 ] ], [ [ 203.302071, 38.126240 ], [ 203.440944, 38.127943 ], [ 203.444454, 37.886978 ], [ 203.305767, 37.885583 ] ], [ [ 203.163133, 38.124409 ], [ 203.301897, 38.126305 ], [ 203.305601, 37.885507 ], [ 203.167375, 37.883898 ] ], [ [ 203.025506, 38.121902 ], [ 203.163129, 38.124121 ], [ 203.167029, 37.883669 ], [ 203.029392, 37.882034 ] ], [ [ 204.278190, 37.653129 ], [ 204.141157, 37.651406 ], [ 204.140619, 37.892334 ], [ 204.277694, 37.892319 ] ], [ [ 204.141193, 37.651639 ], [ 204.003601, 37.650803 ], [ 204.002095, 37.891730 ], [ 204.140271, 37.892370 ] ], [ [ 204.003228, 37.650852 ], [ 203.865344, 37.649775 ], [ 203.863084, 37.890792 ], [ 204.001798, 37.891729 ] ], [ [ 203.864649, 37.649578 ], [ 203.726363, 37.648813 ], [ 203.723565, 37.889968 ], [ 203.862732, 37.890850 ] ], [ [ 203.726211, 37.648811 ], [ 203.587527, 37.647596 ], [ 203.584034, 37.888848 ], [ 203.723413, 37.889979 ] ], [ [ 203.587207, 37.647685 ], [ 203.448735, 37.646586 ], [ 203.444449, 37.887657 ], [ 203.583519, 37.889003 ] ], [ [ 203.448622, 37.646701 ], [ 203.310448, 37.645376 ], [ 203.305751, 37.886291 ], [ 203.444462, 37.887723 ] ], [ [ 203.310153, 37.645449 ], [ 203.172800, 37.644225 ], [ 203.167086, 37.884582 ], [ 203.305299, 37.886294 ] ], [ [ 203.172397, 37.644340 ], [ 203.035396, 37.643205 ], [ 203.029606, 37.882859 ], [ 203.167158, 37.884629 ] ], [ [ 204.276831, 37.394701 ], [ 204.140614, 37.393603 ], [ 204.141237, 37.631464 ], [ 204.277979, 37.632851 ] ], [ [ 204.140537, 37.394056 ], [ 204.004165, 37.392192 ], [ 204.003523, 37.630742 ], [ 204.141220, 37.631727 ] ], [ [ 204.003975, 37.392120 ], [ 203.866662, 37.390293 ], [ 203.865242, 37.629713 ], [ 204.003192, 37.630821 ] ], [ [ 203.866887, 37.390430 ], [ 203.728997, 37.389099 ], [ 203.726561, 37.628538 ], [ 203.864919, 37.629508 ] ], [ [ 203.728436, 37.389677 ], [ 203.591426, 37.388507 ], [ 203.587749, 37.627440 ], [ 203.726306, 37.628571 ] ], [ [ 203.591082, 37.387963 ], [ 203.453362, 37.387171 ], [ 203.448975, 37.626208 ], [ 203.587592, 37.627430 ] ], [ [ 203.453424, 37.387666 ], [ 203.316493, 37.386969 ], [ 203.310802, 37.625213 ], [ 203.448627, 37.626408 ] ], [ [ 203.315969, 37.387201 ], [ 203.179260, 37.385504 ], [ 203.172934, 37.623854 ], [ 203.310389, 37.625287 ] ], [ [ 203.178883, 37.385643 ], [ 203.043045, 37.385013 ], [ 203.035972, 37.623148 ], [ 203.172802, 37.624128 ] ] ] } }, { "Catalog_Number": { "name": "Catalog_Number", "datatype": "int", "ucd": "meta.record;meta.table;meta.file", "value": 4 }, "Catalog_Name": { "name": "Catalog_Name", "datatype": "string", "ucd": "meta.id;meta.table;meta.file", "value": "744334p.cat" }, "Image_Ident": { "name": "Image_Ident", "datatype": "string", "ucd": "meta.id;obs.field", "value": "13h" }, "NExtensions": { "name": "NExtensions", "datatype": "int", "ucd": "meta.record", "value": 36 }, "NAxis": { "name": "NAxis", "datatype": "int", "ucd": "pos.wcs.naxis", "value": 2 }, "Lng_Axis": { "name": "Lng_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.ra", "value": 0 }, "Lat_Axis": { "name": "Lat_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.dec", "value": 1 }, "Ext_Header": { "name": "Ext_Header", "datatype": "boolean", "ucd": "meta.code", "value": false }, "NDetect": { "name": "NDetect", "datatype": "int", "ucd": "meta.number;src", "value": 29519 }, "Group": { "name": "Group", "datatype": "int", "ucd": "meta.id.parent;meta.dataset", "value": 1 }, "Astr_Instrum": { "name": "Astr_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "A1" }, "Phot_Instrum": { "name": "Phot_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "P1" }, "Photom_Flag": { "name": "Photom_Flag", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Photom_Link": { "name": "Photom_Link", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Observation_Date": { "name": "Observation_Date", "datatype": "double", "ucd": "time.epoch;obs.field", "unit": "yr", "value": 2004.353064 }, "Exposure_Time": { "name": "Exposure_Time", "datatype": "float", "ucd": "time.duration;obs.exposure", "value": 512.234000 }, "Air_Mass": { "name": "Air_Mass", "datatype": "float", "ucd": "obs.airMass", "value": 1.190000 }, "Field_Coordinates": { "name": "Field_Coordinates", "datatype": "double array", "ucd": "pos.eq.ra;pos.eq.dec;obs.field", "unit": "%s", "value": [ 203.658406, 37.921448 ] }, "Pixel_Scale": { "name": "Pixel_Scale", "datatype": "float array", "ucd": "instr.scale;instr.pixel;stat.mean", "unit": "%s", "value": [ 0.186017, 0.186017 ] }, "Max_Radius": { "name": "Max_Radius", "datatype": "float", "ucd": "phys.size.radius", "unit": "%s", "value": 42.654234 }, "ZeroPoint_Corr": { "name": "ZeroPoint_Corr", "datatype": "float", "ucd": "phot.mag;phot.calib;arith.zp", "unit": "mag", "value": 0.002140 }, "DPixel_Scale": { "name": "DPixel_Scale", "datatype": "float", "ucd": "instr.scale;instr.pixel;arith.ratio", "value": 0.999902 }, "DPos_Angle": { "name": "DPos_Angle", "datatype": "float", "ucd": "pos.posAng;obs.image;arith.diff", "unit": "deg", "value": 0.003476 }, "AS_Contrast": { "name": "AS_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 30.621208 }, "DX": { "name": "DX", "datatype": "float", "ucd": "pos.eq.ra;arith.diff", "unit": "deg", "value": 0.000198 }, "DY": { "name": "DY", "datatype": "float", "ucd": "pos.eq.dec;arith.diff", "unit": "deg", "value": 0.000014 }, "XY_Contrast": { "name": "XY_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 31.006329 }, "Shear": { "name": "Shear", "datatype": "float", "ucd": "phys.size.axisRatio;obs.image", "value": 0.000038 }, "Shear_PosAngle": { "name": "Shear_PosAngle", "datatype": "float", "ucd": "pos.posAng;obs.image", "unit": "deg", "value": 74.910965 }, "Chi2_Internal": { "name": "Chi2_Internal", "datatype": "float", "ucd": "stat.fit.chi2", "value": 86.744039 }, "NDeg_Internal": { "name": "NDeg_Internal", "datatype": "int", "ucd": "stat.fit.dof", "value": 153347 }, "Chi2_Internal_HighSN": { "name": "Chi2_Internal_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 425.346816 }, "NDeg_Internal_HighSN": { "name": "NDeg_Internal_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 16591 }, "AstromOffset_Reference": { "name": "AstromOffset_Reference", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.003200, -0.004815 ] }, "Astrom_Reference": { "name": "AstromSigma_Reference", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.304469, 0.334024 ] }, "AstromCorr_Reference": { "name": "AstromCorr_Reference", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.002981 }, "Chi2_Reference": { "name": "Chi2_Reference", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.119264 }, "NDeg_Reference": { "name": "NDeg_Reference", "datatype": "int", "ucd": "stat.fit.dof", "value": 800 }, "AstromOffset_Reference_HighSN": { "name": "AstromOffset_Reference_HighSN", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.002486, -0.005183 ] }, "AstromSigma_Reference_HighSN": { "name": "AstromSigma_Reference_HighSN", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.304723, 0.333221 ] }, "AstromCorr_Reference_HighSN": { "name": "AstromCorr_Reference_HighSN", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.004477 }, "Chi2_Reference_HighSN": { "name": "Chi2_Reference_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.108740 }, "NDeg_Reference_HighSN": { "name": "NDeg_Reference_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 797 }, "Set_Polygon": { "name": "Set_Polygon", "datatype": "float array", "ucd": "qsldfjklqksj", "unit": "%s", "value": [ [ [ 204.139522, 38.423996 ], [ 204.277265, 38.423531 ], [ 204.280800, 38.185438 ], [ 204.143091, 38.185621 ] ], [ [ 204.000800, 38.423954 ], [ 204.139423, 38.424031 ], [ 204.142624, 38.185964 ], [ 204.004585, 38.185342 ] ], [ [ 203.861871, 38.424110 ], [ 204.000753, 38.424087 ], [ 204.003991, 38.185303 ], [ 203.864771, 38.184903 ] ], [ [ 203.721844, 38.423435 ], [ 203.861236, 38.423794 ], [ 203.864824, 38.184819 ], [ 203.725120, 38.183951 ] ], [ [ 203.582300, 38.422232 ], [ 203.721743, 38.423313 ], [ 203.724883, 38.183911 ], [ 203.585259, 38.183236 ] ], [ [ 203.442546, 38.420780 ], [ 203.582106, 38.422092 ], [ 203.584982, 38.183068 ], [ 203.445644, 38.181449 ] ], [ [ 203.303539, 38.418658 ], [ 203.442350, 38.421100 ], [ 203.445100, 38.181479 ], [ 203.306165, 38.179835 ] ], [ [ 203.165096, 38.415856 ], [ 203.302723, 38.418387 ], [ 203.306116, 38.179582 ], [ 203.167253, 38.177813 ] ], [ [ 203.027418, 38.412810 ], [ 203.164998, 38.415800 ], [ 203.167122, 38.178035 ], [ 203.029232, 38.175600 ] ], [ [ 204.143355, 38.165416 ], [ 204.280715, 38.165071 ], [ 204.282775, 37.925132 ], [ 204.145558, 37.925070 ] ], [ [ 204.004574, 38.165135 ], [ 204.143101, 38.165483 ], [ 204.145096, 37.924958 ], [ 204.006900, 37.924307 ] ], [ [ 203.865089, 38.164628 ], [ 204.004229, 38.165094 ], [ 204.006718, 37.924287 ], [ 203.867854, 37.923446 ] ], [ [ 203.725364, 38.163784 ], [ 203.864938, 38.164565 ], [ 203.867572, 37.923490 ], [ 203.728312, 37.922559 ] ], [ [ 203.585584, 38.162963 ], [ 203.725113, 38.163876 ], [ 203.728133, 37.922735 ], [ 203.588725, 37.921596 ] ], [ [ 203.445735, 38.161287 ], [ 203.585290, 38.162792 ], [ 203.588412, 37.921455 ], [ 203.449325, 37.920231 ] ], [ [ 203.306472, 38.159835 ], [ 203.445473, 38.161465 ], [ 203.449021, 37.920304 ], [ 203.310206, 37.918983 ] ], [ [ 203.167627, 38.157855 ], [ 203.306487, 38.159784 ], [ 203.310024, 37.918983 ], [ 203.171704, 37.917341 ] ], [ [ 203.029902, 38.155401 ], [ 203.167381, 38.157656 ], [ 203.171285, 37.917204 ], [ 203.033793, 37.915526 ] ], [ [ 204.282985, 37.686535 ], [ 204.146068, 37.684919 ], [ 204.145643, 37.925796 ], [ 204.282607, 37.925680 ] ], [ [ 204.145935, 37.685203 ], [ 204.008305, 37.684304 ], [ 204.007070, 37.925141 ], [ 204.145282, 37.925843 ] ], [ [ 204.008018, 37.684225 ], [ 203.870014, 37.683257 ], [ 203.867740, 37.924382 ], [ 204.006578, 37.925210 ] ], [ [ 203.869435, 37.683134 ], [ 203.730977, 37.682255 ], [ 203.728170, 37.923309 ], [ 203.867511, 37.924305 ] ], [ [ 203.730820, 37.682312 ], [ 203.592143, 37.681105 ], [ 203.588706, 37.922267 ], [ 203.728078, 37.923392 ] ], [ [ 203.591844, 37.681231 ], [ 203.453335, 37.680047 ], [ 203.448987, 37.921055 ], [ 203.588094, 37.922486 ] ], [ [ 203.453207, 37.680248 ], [ 203.314956, 37.678883 ], [ 203.310157, 37.919670 ], [ 203.448946, 37.921144 ] ], [ [ 203.314659, 37.678955 ], [ 203.177140, 37.677725 ], [ 203.171514, 37.918051 ], [ 203.309896, 37.919771 ] ], [ [ 203.176843, 37.677702 ], [ 203.039924, 37.676488 ], [ 203.034015, 37.916283 ], [ 203.171486, 37.918131 ] ], [ [ 204.281677, 37.428332 ], [ 204.145553, 37.427039 ], [ 204.146150, 37.664796 ], [ 204.282795, 37.666384 ] ], [ [ 204.145325, 37.427431 ], [ 204.008914, 37.425607 ], [ 204.008419, 37.664243 ], [ 204.146157, 37.665188 ] ], [ [ 204.008723, 37.425737 ], [ 203.871381, 37.423987 ], [ 203.869995, 37.663205 ], [ 204.007975, 37.664238 ] ], [ [ 203.871490, 37.423995 ], [ 203.733626, 37.422813 ], [ 203.731324, 37.661990 ], [ 203.869651, 37.662810 ] ], [ [ 203.733113, 37.423109 ], [ 203.596072, 37.421919 ], [ 203.592376, 37.660947 ], [ 203.730967, 37.662098 ] ], [ [ 203.595776, 37.421642 ], [ 203.458050, 37.420777 ], [ 203.453477, 37.659564 ], [ 203.592096, 37.660860 ] ], [ [ 203.457846, 37.421054 ], [ 203.320904, 37.420189 ], [ 203.315459, 37.658795 ], [ 203.453298, 37.660167 ] ], [ [ 203.320482, 37.420420 ], [ 203.183717, 37.419155 ], [ 203.177411, 37.657656 ], [ 203.314921, 37.658661 ] ], [ [ 203.183357, 37.419234 ], [ 203.047511, 37.418612 ], [ 203.040316, 37.656483 ], [ 203.177157, 37.657461 ] ] ] } }, { "Catalog_Number": { "name": "Catalog_Number", "datatype": "int", "ucd": "meta.record;meta.table;meta.file", "value": 5 }, "Catalog_Name": { "name": "Catalog_Name", "datatype": "string", "ucd": "meta.id;meta.table;meta.file", "value": "744335p.cat" }, "Image_Ident": { "name": "Image_Ident", "datatype": "string", "ucd": "meta.id;obs.field", "value": "13h" }, "NExtensions": { "name": "NExtensions", "datatype": "int", "ucd": "meta.record", "value": 36 }, "NAxis": { "name": "NAxis", "datatype": "int", "ucd": "pos.wcs.naxis", "value": 2 }, "Lng_Axis": { "name": "Lng_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.ra", "value": 0 }, "Lat_Axis": { "name": "Lat_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.dec", "value": 1 }, "Ext_Header": { "name": "Ext_Header", "datatype": "boolean", "ucd": "meta.code", "value": false }, "NDetect": { "name": "NDetect", "datatype": "int", "ucd": "meta.number;src", "value": 29174 }, "Group": { "name": "Group", "datatype": "int", "ucd": "meta.id.parent;meta.dataset", "value": 1 }, "Astr_Instrum": { "name": "Astr_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "A1" }, "Phot_Instrum": { "name": "Phot_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "P1" }, "Photom_Flag": { "name": "Photom_Flag", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Photom_Link": { "name": "Photom_Link", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Observation_Date": { "name": "Observation_Date", "datatype": "double", "ucd": "time.epoch;obs.field", "unit": "yr", "value": 2004.353082 }, "Exposure_Time": { "name": "Exposure_Time", "datatype": "float", "ucd": "time.duration;obs.exposure", "value": 512.218000 }, "Air_Mass": { "name": "Air_Mass", "datatype": "float", "ucd": "obs.airMass", "value": 1.214000 }, "Field_Coordinates": { "name": "Field_Coordinates", "datatype": "double array", "ucd": "pos.eq.ra;pos.eq.dec;obs.field", "unit": "%s", "value": [ 203.647564, 37.938111 ] }, "Pixel_Scale": { "name": "Pixel_Scale", "datatype": "float array", "ucd": "instr.scale;instr.pixel;stat.mean", "unit": "%s", "value": [ 0.185978, 0.185978 ] }, "Max_Radius": { "name": "Max_Radius", "datatype": "float", "ucd": "phys.size.radius", "unit": "%s", "value": 42.654239 }, "ZeroPoint_Corr": { "name": "ZeroPoint_Corr", "datatype": "float", "ucd": "phot.mag;phot.calib;arith.zp", "unit": "mag", "value": 0.002283 }, "DPixel_Scale": { "name": "DPixel_Scale", "datatype": "float", "ucd": "instr.scale;instr.pixel;arith.ratio", "value": 0.999883 }, "DPos_Angle": { "name": "DPos_Angle", "datatype": "float", "ucd": "pos.posAng;obs.image;arith.diff", "unit": "deg", "value": 0.003545 }, "AS_Contrast": { "name": "AS_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 27.937881 }, "DX": { "name": "DX", "datatype": "float", "ucd": "pos.eq.ra;arith.diff", "unit": "deg", "value": 0.000200 }, "DY": { "name": "DY", "datatype": "float", "ucd": "pos.eq.dec;arith.diff", "unit": "deg", "value": 0.000011 }, "XY_Contrast": { "name": "XY_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 31.225740 }, "Shear": { "name": "Shear", "datatype": "float", "ucd": "phys.size.axisRatio;obs.image", "value": 0.000043 }, "Shear_PosAngle": { "name": "Shear_PosAngle", "datatype": "float", "ucd": "pos.posAng;obs.image", "unit": "deg", "value": 56.434708 }, "Chi2_Internal": { "name": "Chi2_Internal", "datatype": "float", "ucd": "stat.fit.chi2", "value": 98.651058 }, "NDeg_Internal": { "name": "NDeg_Internal", "datatype": "int", "ucd": "stat.fit.dof", "value": 150401 }, "Chi2_Internal_HighSN": { "name": "Chi2_Internal_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 493.302022 }, "NDeg_Internal_HighSN": { "name": "NDeg_Internal_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 16308 }, "AstromOffset_Reference": { "name": "AstromOffset_Reference", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ 0.000113, -0.002746 ] }, "Astrom_Reference": { "name": "AstromSigma_Reference", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.303426, 0.351304 ] }, "AstromCorr_Reference": { "name": "AstromCorr_Reference", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.060508 }, "Chi2_Reference": { "name": "Chi2_Reference", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.140107 }, "NDeg_Reference": { "name": "NDeg_Reference", "datatype": "int", "ucd": "stat.fit.dof", "value": 749 }, "AstromOffset_Reference_HighSN": { "name": "AstromOffset_Reference_HighSN", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ 0.001570, -0.003465 ] }, "AstromSigma_Reference_HighSN": { "name": "AstromSigma_Reference_HighSN", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.303414, 0.349615 ] }, "AstromCorr_Reference_HighSN": { "name": "AstromCorr_Reference_HighSN", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.061476 }, "Chi2_Reference_HighSN": { "name": "Chi2_Reference_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.120611 }, "NDeg_Reference_HighSN": { "name": "NDeg_Reference_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 744 }, "Set_Polygon": { "name": "Set_Polygon", "datatype": "float array", "ucd": "qsldfjklqksj", "unit": "%s", "value": [ [ [ 204.128820, 38.440667 ], [ 204.266458, 38.440177 ], [ 204.270102, 38.202112 ], [ 204.132498, 38.202319 ] ], [ [ 203.989902, 38.440801 ], [ 204.128761, 38.440517 ], [ 204.132035, 38.202300 ], [ 203.993764, 38.202037 ] ], [ [ 203.850723, 38.440569 ], [ 203.989617, 38.440382 ], [ 203.993370, 38.202074 ], [ 203.854158, 38.201844 ] ], [ [ 203.711051, 38.439990 ], [ 203.850244, 38.440380 ], [ 203.853966, 38.201643 ], [ 203.714460, 38.200749 ] ], [ [ 203.571470, 38.438873 ], [ 203.710936, 38.440086 ], [ 203.714014, 38.200543 ], [ 203.574365, 38.199734 ] ], [ [ 203.431593, 38.437332 ], [ 203.571279, 38.438574 ], [ 203.574170, 38.199797 ], [ 203.434708, 38.198248 ] ], [ [ 203.292622, 38.435355 ], [ 203.431459, 38.437838 ], [ 203.434242, 38.198153 ], [ 203.295281, 38.196467 ] ], [ [ 203.154347, 38.432401 ], [ 203.291708, 38.435049 ], [ 203.295078, 38.196370 ], [ 203.156494, 38.194466 ] ], [ [ 203.016346, 38.429689 ], [ 203.154046, 38.432483 ], [ 203.156269, 38.194465 ], [ 203.018249, 38.192227 ] ], [ [ 204.132636, 38.182097 ], [ 204.269948, 38.181644 ], [ 204.272163, 37.941872 ], [ 204.135001, 37.941925 ] ], [ [ 203.993713, 38.181783 ], [ 204.132394, 38.182149 ], [ 204.134482, 37.941658 ], [ 203.996132, 37.940988 ] ], [ [ 203.854308, 38.181308 ], [ 203.993487, 38.181830 ], [ 203.995958, 37.940953 ], [ 203.857056, 37.940056 ] ], [ [ 203.714512, 38.180428 ], [ 203.854130, 38.181134 ], [ 203.856788, 37.940218 ], [ 203.717484, 37.939363 ] ], [ [ 203.574717, 38.179676 ], [ 203.714240, 38.180462 ], [ 203.717302, 37.939366 ], [ 203.577901, 37.938355 ] ], [ [ 203.434829, 38.177930 ], [ 203.574447, 38.179457 ], [ 203.577565, 37.938110 ], [ 203.438415, 37.936864 ] ], [ [ 203.295545, 38.176497 ], [ 203.434478, 38.178172 ], [ 203.438165, 37.936985 ], [ 203.299418, 37.935619 ] ], [ [ 203.156588, 38.174480 ], [ 203.295488, 38.176344 ], [ 203.299108, 37.935676 ], [ 203.160747, 37.934098 ] ], [ [ 203.018878, 38.171997 ], [ 203.156423, 38.174239 ], [ 203.160389, 37.933848 ], [ 203.022830, 37.932185 ] ], [ [ 204.272273, 37.703085 ], [ 204.135335, 37.701716 ], [ 204.135038, 37.942607 ], [ 204.272025, 37.942244 ] ], [ [ 204.135264, 37.701881 ], [ 203.997709, 37.700962 ], [ 203.996307, 37.941792 ], [ 204.134444, 37.942512 ] ], [ [ 203.997338, 37.700924 ], [ 203.859216, 37.699863 ], [ 203.856954, 37.940992 ], [ 203.995909, 37.941913 ] ], [ [ 203.858683, 37.699771 ], [ 203.720187, 37.699001 ], [ 203.717340, 37.940017 ], [ 203.856720, 37.940905 ] ], [ [ 203.720088, 37.698958 ], [ 203.581307, 37.697787 ], [ 203.577745, 37.938951 ], [ 203.717224, 37.940039 ] ], [ [ 203.580999, 37.697897 ], [ 203.442542, 37.696782 ], [ 203.438118, 37.937748 ], [ 203.577172, 37.939107 ] ], [ [ 203.442233, 37.696844 ], [ 203.303965, 37.695632 ], [ 203.299285, 37.936412 ], [ 203.438089, 37.937732 ] ], [ [ 203.303749, 37.695637 ], [ 203.166170, 37.694329 ], [ 203.160520, 37.934696 ], [ 203.298964, 37.936495 ] ], [ [ 203.165883, 37.694389 ], [ 203.028995, 37.693323 ], [ 203.023029, 37.932973 ], [ 203.160467, 37.934671 ] ], [ [ 204.270878, 37.445063 ], [ 204.134829, 37.443692 ], [ 204.135452, 37.681369 ], [ 204.272019, 37.683039 ] ], [ [ 204.134570, 37.444101 ], [ 203.998126, 37.442312 ], [ 203.997646, 37.680919 ], [ 204.135418, 37.681828 ] ], [ [ 203.997930, 37.442194 ], [ 203.860765, 37.440508 ], [ 203.859281, 37.679928 ], [ 203.997087, 37.680895 ] ], [ [ 203.860738, 37.440732 ], [ 203.722724, 37.439532 ], [ 203.720398, 37.678599 ], [ 203.858877, 37.679437 ] ], [ [ 203.722396, 37.439696 ], [ 203.585202, 37.438549 ], [ 203.581380, 37.677700 ], [ 203.720128, 37.678807 ] ], [ [ 203.585058, 37.438403 ], [ 203.447155, 37.437524 ], [ 203.442456, 37.676137 ], [ 203.581251, 37.677447 ] ], [ [ 203.447027, 37.437784 ], [ 203.310012, 37.436925 ], [ 203.304503, 37.675355 ], [ 203.442415, 37.676718 ] ], [ [ 203.309545, 37.437081 ], [ 203.172895, 37.435974 ], [ 203.166532, 37.674378 ], [ 203.303926, 37.675225 ] ], [ [ 203.172489, 37.435952 ], [ 203.036488, 37.435249 ], [ 203.029265, 37.673131 ], [ 203.166260, 37.674184 ] ] ] } }, { "Catalog_Number": { "name": "Catalog_Number", "datatype": "int", "ucd": "meta.record;meta.table;meta.file", "value": 6 }, "Catalog_Name": { "name": "Catalog_Name", "datatype": "string", "ucd": "meta.id;meta.table;meta.file", "value": "744336p.cat" }, "Image_Ident": { "name": "Image_Ident", "datatype": "string", "ucd": "meta.id;obs.field", "value": "13h" }, "NExtensions": { "name": "NExtensions", "datatype": "int", "ucd": "meta.record", "value": 36 }, "NAxis": { "name": "NAxis", "datatype": "int", "ucd": "pos.wcs.naxis", "value": 2 }, "Lng_Axis": { "name": "Lng_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.ra", "value": 0 }, "Lat_Axis": { "name": "Lat_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.dec", "value": 1 }, "Ext_Header": { "name": "Ext_Header", "datatype": "boolean", "ucd": "meta.code", "value": false }, "NDetect": { "name": "NDetect", "datatype": "int", "ucd": "meta.number;src", "value": 28830 }, "Group": { "name": "Group", "datatype": "int", "ucd": "meta.id.parent;meta.dataset", "value": 1 }, "Astr_Instrum": { "name": "Astr_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "A1" }, "Phot_Instrum": { "name": "Phot_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "P1" }, "Photom_Flag": { "name": "Photom_Flag", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Photom_Link": { "name": "Photom_Link", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Observation_Date": { "name": "Observation_Date", "datatype": "double", "ucd": "time.epoch;obs.field", "unit": "yr", "value": 2004.353099 }, "Exposure_Time": { "name": "Exposure_Time", "datatype": "float", "ucd": "time.duration;obs.exposure", "value": 512.236000 }, "Air_Mass": { "name": "Air_Mass", "datatype": "float", "ucd": "obs.airMass", "value": 1.241000 }, "Field_Coordinates": { "name": "Field_Coordinates", "datatype": "double array", "ucd": "pos.eq.ra;pos.eq.dec;obs.field", "unit": "%s", "value": [ 203.645346, 37.929610 ] }, "Pixel_Scale": { "name": "Pixel_Scale", "datatype": "float array", "ucd": "instr.scale;instr.pixel;stat.mean", "unit": "%s", "value": [ 0.185971, 0.185971 ] }, "Max_Radius": { "name": "Max_Radius", "datatype": "float", "ucd": "phys.size.radius", "unit": "%s", "value": 42.653297 }, "ZeroPoint_Corr": { "name": "ZeroPoint_Corr", "datatype": "float", "ucd": "phot.mag;phot.calib;arith.zp", "unit": "mag", "value": 0.001876 }, "DPixel_Scale": { "name": "DPixel_Scale", "datatype": "float", "ucd": "instr.scale;instr.pixel;arith.ratio", "value": 0.999907 }, "DPos_Angle": { "name": "DPos_Angle", "datatype": "float", "ucd": "pos.posAng;obs.image;arith.diff", "unit": "deg", "value": 0.003225 }, "AS_Contrast": { "name": "AS_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 28.232405 }, "DX": { "name": "DX", "datatype": "float", "ucd": "pos.eq.ra;arith.diff", "unit": "deg", "value": 0.000195 }, "DY": { "name": "DY", "datatype": "float", "ucd": "pos.eq.dec;arith.diff", "unit": "deg", "value": 0.000015 }, "XY_Contrast": { "name": "XY_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 31.553288 }, "Shear": { "name": "Shear", "datatype": "float", "ucd": "phys.size.axisRatio;obs.image", "value": 0.000027 }, "Shear_PosAngle": { "name": "Shear_PosAngle", "datatype": "float", "ucd": "pos.posAng;obs.image", "unit": "deg", "value": 85.449661 }, "Chi2_Internal": { "name": "Chi2_Internal", "datatype": "float", "ucd": "stat.fit.chi2", "value": 96.263789 }, "NDeg_Internal": { "name": "NDeg_Internal", "datatype": "int", "ucd": "stat.fit.dof", "value": 148349 }, "Chi2_Internal_HighSN": { "name": "Chi2_Internal_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 469.545961 }, "NDeg_Internal_HighSN": { "name": "NDeg_Internal_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 16109 }, "AstromOffset_Reference": { "name": "AstromOffset_Reference", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.001220, 0.002747 ] }, "Astrom_Reference": { "name": "AstromSigma_Reference", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.302599, 0.357099 ] }, "AstromCorr_Reference": { "name": "AstromCorr_Reference", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.030713 }, "Chi2_Reference": { "name": "Chi2_Reference", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.134238 }, "NDeg_Reference": { "name": "NDeg_Reference", "datatype": "int", "ucd": "stat.fit.dof", "value": 752 }, "AstromOffset_Reference_HighSN": { "name": "AstromOffset_Reference_HighSN", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ 0.001403, 0.002390 ] }, "AstromSigma_Reference_HighSN": { "name": "AstromSigma_Reference_HighSN", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.301801, 0.355210 ] }, "AstromCorr_Reference_HighSN": { "name": "AstromCorr_Reference_HighSN", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.027330 }, "Chi2_Reference_HighSN": { "name": "Chi2_Reference_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.115596 }, "NDeg_Reference_HighSN": { "name": "NDeg_Reference_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 744 }, "Set_Polygon": { "name": "Set_Polygon", "datatype": "float array", "ucd": "qsldfjklqksj", "unit": "%s", "value": [ [ [ 204.126380, 38.432303 ], [ 204.264025, 38.431819 ], [ 204.267605, 38.193722 ], [ 204.129995, 38.193923 ] ], [ [ 203.987496, 38.432369 ], [ 204.126302, 38.432138 ], [ 204.129616, 38.194040 ], [ 203.991398, 38.193724 ] ], [ [ 203.855643, 38.429018 ], [ 203.995123, 38.428244 ], [ 203.994881, 38.187684 ], [ 203.854954, 38.188010 ] ], [ [ 203.708659, 38.431626 ], [ 203.847919, 38.432117 ], [ 203.851559, 38.193205 ], [ 203.711986, 38.192204 ] ], [ [ 203.569088, 38.430493 ], [ 203.708504, 38.431666 ], [ 203.711697, 38.192277 ], [ 203.572100, 38.191509 ] ], [ [ 203.429272, 38.429137 ], [ 203.568904, 38.430334 ], [ 203.571751, 38.191305 ], [ 203.432341, 38.189801 ] ], [ [ 203.290321, 38.426658 ], [ 203.429115, 38.429117 ], [ 203.431837, 38.189897 ], [ 203.292920, 38.188233 ] ], [ [ 203.151979, 38.424167 ], [ 203.289467, 38.426643 ], [ 203.292794, 38.187956 ], [ 203.154077, 38.186234 ] ], [ [ 203.013974, 38.421199 ], [ 203.151596, 38.424074 ], [ 203.153879, 38.186204 ], [ 203.015939, 38.183885 ] ], [ [ 204.130249, 38.173776 ], [ 204.267557, 38.173371 ], [ 204.269617, 37.933414 ], [ 204.132455, 37.933415 ] ], [ [ 203.991377, 38.173515 ], [ 204.129967, 38.173759 ], [ 204.132003, 37.933245 ], [ 203.993745, 37.932697 ] ], [ [ 203.851884, 38.172970 ], [ 203.991077, 38.173427 ], [ 203.993542, 37.932579 ], [ 203.854626, 37.931747 ] ], [ [ 203.712097, 38.172067 ], [ 203.851693, 38.172759 ], [ 203.854404, 37.931882 ], [ 203.715120, 37.931040 ] ], [ [ 203.572400, 38.171362 ], [ 203.711906, 38.172098 ], [ 203.714866, 37.930936 ], [ 203.575481, 37.929975 ] ], [ [ 203.432504, 38.169598 ], [ 203.572062, 38.171126 ], [ 203.575175, 37.929758 ], [ 203.436085, 37.928510 ] ], [ [ 203.293225, 38.168067 ], [ 203.432182, 38.169844 ], [ 203.435775, 37.928706 ], [ 203.297005, 37.927238 ] ], [ [ 203.154303, 38.166111 ], [ 203.293255, 38.167945 ], [ 203.296737, 37.927316 ], [ 203.158326, 37.925768 ] ], [ [ 203.016538, 38.163729 ], [ 203.154062, 38.165935 ], [ 203.158057, 37.925425 ], [ 203.020518, 37.923796 ] ], [ [ 204.269806, 37.694704 ], [ 204.133013, 37.693415 ], [ 204.132650, 37.934314 ], [ 204.269497, 37.933875 ] ], [ [ 204.132787, 37.693521 ], [ 203.995197, 37.692621 ], [ 203.993917, 37.933437 ], [ 204.132090, 37.934138 ] ], [ [ 203.994900, 37.692575 ], [ 203.856823, 37.691495 ], [ 203.854515, 37.932586 ], [ 203.993425, 37.933526 ] ], [ [ 203.856237, 37.691460 ], [ 203.717847, 37.690669 ], [ 203.714952, 37.931617 ], [ 203.854225, 37.932525 ] ], [ [ 203.717675, 37.690579 ], [ 203.578896, 37.689400 ], [ 203.575407, 37.930603 ], [ 203.714883, 37.931700 ] ], [ [ 203.578623, 37.689514 ], [ 203.440193, 37.688416 ], [ 203.435769, 37.929401 ], [ 203.574795, 37.930744 ] ], [ [ 203.439850, 37.688497 ], [ 203.301624, 37.687219 ], [ 203.296961, 37.928028 ], [ 203.435724, 37.929414 ] ], [ [ 203.301399, 37.687239 ], [ 203.163860, 37.686060 ], [ 203.158183, 37.926298 ], [ 203.296583, 37.927966 ] ], [ [ 203.163502, 37.685956 ], [ 203.026641, 37.684891 ], [ 203.020760, 37.924636 ], [ 203.158172, 37.926334 ] ], [ [ 204.268441, 37.436702 ], [ 204.132396, 37.435334 ], [ 204.133034, 37.673019 ], [ 204.269596, 37.674686 ] ], [ [ 204.132251, 37.435528 ], [ 203.995817, 37.433762 ], [ 203.995143, 37.672616 ], [ 204.132906, 37.673498 ] ], [ [ 203.995396, 37.434095 ], [ 203.858167, 37.432335 ], [ 203.856925, 37.671409 ], [ 203.994793, 37.672455 ] ], [ [ 203.858339, 37.432312 ], [ 203.720279, 37.431014 ], [ 203.717989, 37.670225 ], [ 203.856515, 37.671159 ] ], [ [ 203.720051, 37.431331 ], [ 203.582853, 37.430198 ], [ 203.579032, 37.669355 ], [ 203.717785, 37.670447 ] ], [ [ 203.582572, 37.430010 ], [ 203.444783, 37.429117 ], [ 203.440134, 37.667783 ], [ 203.578814, 37.669107 ] ], [ [ 203.444642, 37.429339 ], [ 203.307706, 37.428687 ], [ 203.302203, 37.667217 ], [ 203.440034, 37.668370 ] ], [ [ 203.307115, 37.428594 ], [ 203.170494, 37.427450 ], [ 203.164219, 37.666039 ], [ 203.301585, 37.666919 ] ], [ [ 203.170138, 37.427542 ], [ 203.034211, 37.426961 ], [ 203.026971, 37.664792 ], [ 203.163892, 37.665727 ] ] ] } }, { "Catalog_Number": { "name": "Catalog_Number", "datatype": "int", "ucd": "meta.record;meta.table;meta.file", "value": 7 }, "Catalog_Name": { "name": "Catalog_Name", "datatype": "string", "ucd": "meta.id;meta.table;meta.file", "value": "744337p.cat" }, "Image_Ident": { "name": "Image_Ident", "datatype": "string", "ucd": "meta.id;obs.field", "value": "13h" }, "NExtensions": { "name": "NExtensions", "datatype": "int", "ucd": "meta.record", "value": 36 }, "NAxis": { "name": "NAxis", "datatype": "int", "ucd": "pos.wcs.naxis", "value": 2 }, "Lng_Axis": { "name": "Lng_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.ra", "value": 0 }, "Lat_Axis": { "name": "Lat_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.dec", "value": 1 }, "Ext_Header": { "name": "Ext_Header", "datatype": "boolean", "ucd": "meta.code", "value": false }, "NDetect": { "name": "NDetect", "datatype": "int", "ucd": "meta.number;src", "value": 28234 }, "Group": { "name": "Group", "datatype": "int", "ucd": "meta.id.parent;meta.dataset", "value": 1 }, "Astr_Instrum": { "name": "Astr_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "A1" }, "Phot_Instrum": { "name": "Phot_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "P1" }, "Photom_Flag": { "name": "Photom_Flag", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Photom_Link": { "name": "Photom_Link", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Observation_Date": { "name": "Observation_Date", "datatype": "double", "ucd": "time.epoch;obs.field", "unit": "yr", "value": 2004.353117 }, "Exposure_Time": { "name": "Exposure_Time", "datatype": "float", "ucd": "time.duration;obs.exposure", "value": 512.224000 }, "Air_Mass": { "name": "Air_Mass", "datatype": "float", "ucd": "obs.airMass", "value": 1.270000 }, "Field_Coordinates": { "name": "Field_Coordinates", "datatype": "double array", "ucd": "pos.eq.ra;pos.eq.dec;obs.field", "unit": "%s", "value": [ 203.656265, 37.896258 ] }, "Pixel_Scale": { "name": "Pixel_Scale", "datatype": "float array", "ucd": "instr.scale;instr.pixel;stat.mean", "unit": "%s", "value": [ 0.186024, 0.186024 ] }, "Max_Radius": { "name": "Max_Radius", "datatype": "float", "ucd": "phys.size.radius", "unit": "%s", "value": 42.649248 }, "ZeroPoint_Corr": { "name": "ZeroPoint_Corr", "datatype": "float", "ucd": "phot.mag;phot.calib;arith.zp", "unit": "mag", "value": 0.004953 }, "DPixel_Scale": { "name": "DPixel_Scale", "datatype": "float", "ucd": "instr.scale;instr.pixel;arith.ratio", "value": 0.999917 }, "DPos_Angle": { "name": "DPos_Angle", "datatype": "float", "ucd": "pos.posAng;obs.image;arith.diff", "unit": "deg", "value": 0.002670 }, "AS_Contrast": { "name": "AS_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 25.708792 }, "DX": { "name": "DX", "datatype": "float", "ucd": "pos.eq.ra;arith.diff", "unit": "deg", "value": 0.000205 }, "DY": { "name": "DY", "datatype": "float", "ucd": "pos.eq.dec;arith.diff", "unit": "deg", "value": 0.000015 }, "XY_Contrast": { "name": "XY_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 30.607275 }, "Shear": { "name": "Shear", "datatype": "float", "ucd": "phys.size.axisRatio;obs.image", "value": 0.000039 }, "Shear_PosAngle": { "name": "Shear_PosAngle", "datatype": "float", "ucd": "pos.posAng;obs.image", "unit": "deg", "value": 66.127304 }, "Chi2_Internal": { "name": "Chi2_Internal", "datatype": "float", "ucd": "stat.fit.chi2", "value": 86.886943 }, "NDeg_Internal": { "name": "NDeg_Internal", "datatype": "int", "ucd": "stat.fit.dof", "value": 146147 }, "Chi2_Internal_HighSN": { "name": "Chi2_Internal_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 416.548141 }, "NDeg_Internal_HighSN": { "name": "NDeg_Internal_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 15788 }, "AstromOffset_Reference": { "name": "AstromOffset_Reference", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.002305, -0.006675 ] }, "Astrom_Reference": { "name": "AstromSigma_Reference", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.304852, 0.333695 ] }, "AstromCorr_Reference": { "name": "AstromCorr_Reference", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.005730 }, "Chi2_Reference": { "name": "Chi2_Reference", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.265490 }, "NDeg_Reference": { "name": "NDeg_Reference", "datatype": "int", "ucd": "stat.fit.dof", "value": 811 }, "AstromOffset_Reference_HighSN": { "name": "AstromOffset_Reference_HighSN", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.002168, -0.007558 ] }, "AstromSigma_Reference_HighSN": { "name": "AstromSigma_Reference_HighSN", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.303546, 0.330101 ] }, "AstromCorr_Reference_HighSN": { "name": "AstromCorr_Reference_HighSN", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.009887 }, "Chi2_Reference_HighSN": { "name": "Chi2_Reference_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.242532 }, "NDeg_Reference_HighSN": { "name": "NDeg_Reference_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 802 }, "Set_Polygon": { "name": "Set_Polygon", "datatype": "float array", "ucd": "qsldfjklqksj", "unit": "%s", "value": [ [ [ 204.137028, 38.398847 ], [ 204.274623, 38.398393 ], [ 204.278313, 38.160449 ], [ 204.140753, 38.160620 ] ], [ [ 204.005426, 38.395609 ], [ 204.145241, 38.395247 ], [ 204.143624, 38.154977 ], [ 204.004317, 38.154751 ] ], [ [ 203.859225, 38.398855 ], [ 203.998116, 38.398902 ], [ 204.001648, 38.160379 ], [ 203.862425, 38.159910 ] ], [ [ 203.719557, 38.398305 ], [ 203.858841, 38.398851 ], [ 203.862454, 38.159705 ], [ 203.722856, 38.158645 ] ], [ [ 203.580115, 38.397205 ], [ 203.719492, 38.398273 ], [ 203.722583, 38.158763 ], [ 203.583023, 38.158100 ] ], [ [ 203.440358, 38.395662 ], [ 203.579789, 38.397048 ], [ 203.582674, 38.158084 ], [ 203.443465, 38.156391 ] ], [ [ 203.301377, 38.393368 ], [ 203.440118, 38.395876 ], [ 203.442820, 38.156564 ], [ 203.303956, 38.154852 ] ], [ [ 203.163114, 38.390659 ], [ 203.300713, 38.393274 ], [ 203.303870, 38.154607 ], [ 203.165041, 38.152746 ] ], [ [ 203.025143, 38.388146 ], [ 203.162725, 38.390880 ], [ 203.165078, 38.152724 ], [ 203.027172, 38.150546 ] ], [ [ 204.141065, 38.140496 ], [ 204.278411, 38.140106 ], [ 204.280294, 37.900008 ], [ 204.143091, 37.899990 ] ], [ [ 204.002138, 38.140117 ], [ 204.140639, 38.140435 ], [ 204.142705, 37.899872 ], [ 204.004535, 37.899250 ] ], [ [ 203.862745, 38.139544 ], [ 204.001868, 38.140103 ], [ 204.004323, 37.899275 ], [ 203.865476, 37.898341 ] ], [ [ 203.723077, 38.138687 ], [ 203.862546, 38.139510 ], [ 203.865223, 37.898468 ], [ 203.726066, 37.897495 ] ], [ [ 203.583296, 38.137925 ], [ 203.722806, 38.138817 ], [ 203.725838, 37.897639 ], [ 203.586449, 37.896521 ] ], [ [ 203.443580, 38.136225 ], [ 203.583021, 38.137727 ], [ 203.586100, 37.896376 ], [ 203.447126, 37.895154 ] ], [ [ 203.304369, 38.134679 ], [ 203.443259, 38.136402 ], [ 203.446741, 37.895348 ], [ 203.308038, 37.893932 ] ], [ [ 203.165424, 38.132798 ], [ 203.304153, 38.134719 ], [ 203.307859, 37.893926 ], [ 203.169668, 37.892292 ] ], [ [ 203.027724, 38.130338 ], [ 203.165297, 38.132604 ], [ 203.169280, 37.892082 ], [ 203.031693, 37.890398 ] ], [ [ 204.280540, 37.661444 ], [ 204.143609, 37.659835 ], [ 204.143204, 37.900832 ], [ 204.280179, 37.900707 ] ], [ [ 204.143566, 37.660071 ], [ 204.005972, 37.659187 ], [ 204.004605, 37.900118 ], [ 204.142782, 37.900805 ] ], [ [ 204.005568, 37.659205 ], [ 203.867664, 37.658169 ], [ 203.865514, 37.899239 ], [ 204.004248, 37.900134 ] ], [ [ 203.867048, 37.657996 ], [ 203.728724, 37.657224 ], [ 203.725955, 37.898347 ], [ 203.865160, 37.899237 ] ], [ [ 203.728509, 37.657219 ], [ 203.589885, 37.656010 ], [ 203.586443, 37.897219 ], [ 203.725763, 37.898345 ] ], [ [ 203.589594, 37.656138 ], [ 203.451177, 37.655002 ], [ 203.446790, 37.896023 ], [ 203.585803, 37.897405 ] ], [ [ 203.450953, 37.655084 ], [ 203.312748, 37.653807 ], [ 203.308046, 37.894729 ], [ 203.446788, 37.896115 ] ], [ [ 203.312463, 37.653843 ], [ 203.174969, 37.652658 ], [ 203.169398, 37.892973 ], [ 203.307755, 37.894647 ] ], [ [ 203.174713, 37.652694 ], [ 203.037959, 37.651586 ], [ 203.031986, 37.891239 ], [ 203.169289, 37.892978 ] ], [ [ 204.279231, 37.403177 ], [ 204.143001, 37.402044 ], [ 204.143701, 37.639841 ], [ 204.280455, 37.641264 ] ], [ [ 204.142930, 37.402528 ], [ 204.006556, 37.400617 ], [ 204.005966, 37.639068 ], [ 204.143665, 37.640101 ] ], [ [ 204.006383, 37.400606 ], [ 203.869047, 37.398814 ], [ 203.867633, 37.638088 ], [ 204.005606, 37.639162 ] ], [ [ 203.869300, 37.398886 ], [ 203.731304, 37.397503 ], [ 203.728904, 37.636890 ], [ 203.867368, 37.637911 ] ], [ [ 203.730626, 37.398085 ], [ 203.593769, 37.396823 ], [ 203.590235, 37.635822 ], [ 203.728635, 37.637046 ] ], [ [ 203.593453, 37.396681 ], [ 203.455637, 37.395750 ], [ 203.451240, 37.634514 ], [ 203.589951, 37.635876 ] ], [ [ 203.455771, 37.396080 ], [ 203.318807, 37.395341 ], [ 203.313100, 37.633588 ], [ 203.450958, 37.634826 ] ], [ [ 203.318222, 37.395257 ], [ 203.181554, 37.393979 ], [ 203.175319, 37.632637 ], [ 203.312734, 37.633653 ] ], [ [ 203.181261, 37.394095 ], [ 203.045314, 37.393329 ], [ 203.038120, 37.631491 ], [ 203.175059, 37.632605 ] ] ] } }, { "Catalog_Number": { "name": "Catalog_Number", "datatype": "int", "ucd": "meta.record;meta.table;meta.file", "value": 8 }, "Catalog_Name": { "name": "Catalog_Name", "datatype": "string", "ucd": "meta.id;meta.table;meta.file", "value": "744338p.cat" }, "Image_Ident": { "name": "Image_Ident", "datatype": "string", "ucd": "meta.id;obs.field", "value": "13h" }, "NExtensions": { "name": "NExtensions", "datatype": "int", "ucd": "meta.record", "value": 36 }, "NAxis": { "name": "NAxis", "datatype": "int", "ucd": "pos.wcs.naxis", "value": 2 }, "Lng_Axis": { "name": "Lng_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.ra", "value": 0 }, "Lat_Axis": { "name": "Lat_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.dec", "value": 1 }, "Ext_Header": { "name": "Ext_Header", "datatype": "boolean", "ucd": "meta.code", "value": false }, "NDetect": { "name": "NDetect", "datatype": "int", "ucd": "meta.number;src", "value": 27896 }, "Group": { "name": "Group", "datatype": "int", "ucd": "meta.id.parent;meta.dataset", "value": 1 }, "Astr_Instrum": { "name": "Astr_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "A1" }, "Phot_Instrum": { "name": "Phot_Instrum", "datatype": "string", "ucd": "meta.id.parent;meta.dataset", "value": "P1" }, "Photom_Flag": { "name": "Photom_Flag", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Photom_Link": { "name": "Photom_Link", "datatype": "boolean", "ucd": "meta.code;phot", "value": false }, "Observation_Date": { "name": "Observation_Date", "datatype": "double", "ucd": "time.epoch;obs.field", "unit": "yr", "value": 2004.353135 }, "Exposure_Time": { "name": "Exposure_Time", "datatype": "float", "ucd": "time.duration;obs.exposure", "value": 512.232000 }, "Air_Mass": { "name": "Air_Mass", "datatype": "float", "ucd": "obs.airMass", "value": 1.302000 }, "Field_Coordinates": { "name": "Field_Coordinates", "datatype": "double array", "ucd": "pos.eq.ra;pos.eq.dec;obs.field", "unit": "%s", "value": [ 203.654049, 37.933141 ] }, "Pixel_Scale": { "name": "Pixel_Scale", "datatype": "float array", "ucd": "instr.scale;instr.pixel;stat.mean", "unit": "%s", "value": [ 0.185981, 0.185981 ] }, "Max_Radius": { "name": "Max_Radius", "datatype": "float", "ucd": "phys.size.radius", "unit": "%s", "value": 42.657439 }, "ZeroPoint_Corr": { "name": "ZeroPoint_Corr", "datatype": "float", "ucd": "phot.mag;phot.calib;arith.zp", "unit": "mag", "value": 0.002178 }, "DPixel_Scale": { "name": "DPixel_Scale", "datatype": "float", "ucd": "instr.scale;instr.pixel;arith.ratio", "value": 0.999896 }, "DPos_Angle": { "name": "DPos_Angle", "datatype": "float", "ucd": "pos.posAng;obs.image;arith.diff", "unit": "deg", "value": 0.003631 }, "AS_Contrast": { "name": "AS_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 29.235916 }, "DX": { "name": "DX", "datatype": "float", "ucd": "pos.eq.ra;arith.diff", "unit": "deg", "value": 0.000199 }, "DY": { "name": "DY", "datatype": "float", "ucd": "pos.eq.dec;arith.diff", "unit": "deg", "value": 0.000011 }, "XY_Contrast": { "name": "XY_Contrast", "datatype": "float", "ucd": "stat.correlation;arith.ratio", "value": 31.136786 }, "Shear": { "name": "Shear", "datatype": "float", "ucd": "phys.size.axisRatio;obs.image", "value": 0.000049 }, "Shear_PosAngle": { "name": "Shear_PosAngle", "datatype": "float", "ucd": "pos.posAng;obs.image", "unit": "deg", "value": 68.145279 }, "Chi2_Internal": { "name": "Chi2_Internal", "datatype": "float", "ucd": "stat.fit.chi2", "value": 84.044002 }, "NDeg_Internal": { "name": "NDeg_Internal", "datatype": "int", "ucd": "stat.fit.dof", "value": 149257 }, "Chi2_Internal_HighSN": { "name": "Chi2_Internal_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 417.930472 }, "NDeg_Internal_HighSN": { "name": "NDeg_Internal_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 16038 }, "AstromOffset_Reference": { "name": "AstromOffset_Reference", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.002321, -0.002046 ] }, "Astrom_Reference": { "name": "AstromSigma_Reference", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.294253, 0.348964 ] }, "AstromCorr_Reference": { "name": "AstromCorr_Reference", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.068640 }, "Chi2_Reference": { "name": "Chi2_Reference", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.190759 }, "NDeg_Reference": { "name": "NDeg_Reference", "datatype": "int", "ucd": "stat.fit.dof", "value": 820 }, "AstromOffset_Reference_HighSN": { "name": "AstromOffset_Reference_HighSN", "datatype": "float array", "ucd": "pos.eq.ra;pos.eq.dec;arith.diff;obs.field", "unit": "%s", "value": [ -0.001247, -0.002904 ] }, "AstromSigma_Reference_HighSN": { "name": "AstromSigma_Reference_HighSN", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "unit": "%s", "value": [ 0.294294, 0.346732 ] }, "AstromCorr_Reference_HighSN": { "name": "AstromCorr_Reference_HighSN", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.068268 }, "Chi2_Reference_HighSN": { "name": "Chi2_Reference_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.171185 }, "NDeg_Reference_HighSN": { "name": "NDeg_Reference_HighSN", "datatype": "int", "ucd": "stat.fit.dof", "value": 815 }, "Set_Polygon": { "name": "Set_Polygon", "datatype": "float array", "ucd": "qsldfjklqksj", "unit": "%s", "value": [ [ [ 204.135221, 38.435712 ], [ 204.272907, 38.435167 ], [ 204.276621, 38.197132 ], [ 204.138971, 38.197395 ] ], [ [ 203.996514, 38.435681 ], [ 204.135208, 38.435711 ], [ 204.138356, 38.197623 ], [ 204.000246, 38.197048 ] ], [ [ 203.857240, 38.435611 ], [ 203.996238, 38.435680 ], [ 203.999805, 38.197094 ], [ 203.860475, 38.196602 ] ], [ [ 203.717560, 38.435190 ], [ 203.856860, 38.435675 ], [ 203.860419, 38.196510 ], [ 203.720807, 38.195511 ] ], [ [ 203.578047, 38.433906 ], [ 203.717420, 38.435096 ], [ 203.720456, 38.195600 ], [ 203.580901, 38.194812 ] ], [ [ 203.438037, 38.432419 ], [ 203.577720, 38.433724 ], [ 203.580693, 38.194827 ], [ 203.441233, 38.193213 ] ], [ [ 203.299188, 38.430094 ], [ 203.437913, 38.432643 ], [ 203.440621, 38.193282 ], [ 203.301773, 38.191527 ] ], [ [ 203.160653, 38.427599 ], [ 203.298253, 38.430071 ], [ 203.301653, 38.191301 ], [ 203.162819, 38.189591 ] ], [ [ 203.022934, 38.424646 ], [ 203.160506, 38.427418 ], [ 203.162688, 38.189501 ], [ 203.024799, 38.187286 ] ], [ [ 204.139136, 38.177183 ], [ 204.276447, 38.176787 ], [ 204.278542, 37.936797 ], [ 204.141376, 37.936790 ] ], [ [ 204.000335, 38.176874 ], [ 204.138844, 38.177184 ], [ 204.140813, 37.936625 ], [ 204.002638, 37.936013 ] ], [ [ 203.860752, 38.176320 ], [ 203.999920, 38.176802 ], [ 204.002456, 37.935985 ], [ 203.863565, 37.935128 ] ], [ [ 203.720988, 38.175467 ], [ 203.860623, 38.176184 ], [ 203.863288, 37.935249 ], [ 203.723966, 37.934382 ] ], [ [ 203.581203, 38.174682 ], [ 203.720760, 38.175575 ], [ 203.723790, 37.934408 ], [ 203.584355, 37.933289 ] ], [ [ 203.441319, 38.172968 ], [ 203.580908, 38.174471 ], [ 203.584056, 37.933140 ], [ 203.444935, 37.931917 ] ], [ [ 203.302018, 38.171527 ], [ 203.441048, 38.173201 ], [ 203.444659, 37.932008 ], [ 203.305814, 37.930644 ] ], [ [ 203.163118, 38.169523 ], [ 203.302060, 38.171415 ], [ 203.305640, 37.930689 ], [ 203.167240, 37.929083 ] ], [ [ 203.025371, 38.167138 ], [ 203.162821, 38.169402 ], [ 203.166834, 37.928863 ], [ 203.029370, 37.927175 ] ], [ [ 204.278700, 37.698084 ], [ 204.141833, 37.696779 ], [ 204.141529, 37.937691 ], [ 204.278447, 37.937266 ] ], [ [ 204.141686, 37.696914 ], [ 204.004069, 37.695995 ], [ 204.002861, 37.936814 ], [ 204.141059, 37.937535 ] ], [ [ 204.003749, 37.695942 ], [ 203.865693, 37.694917 ], [ 203.863451, 37.936044 ], [ 204.002341, 37.936928 ] ], [ [ 203.865186, 37.694817 ], [ 203.726644, 37.693957 ], [ 203.723804, 37.935012 ], [ 203.863229, 37.935990 ] ], [ [ 203.726614, 37.693957 ], [ 203.587788, 37.692813 ], [ 203.584224, 37.934001 ], [ 203.723746, 37.935062 ] ], [ [ 203.587476, 37.692940 ], [ 203.449005, 37.691818 ], [ 203.444593, 37.932776 ], [ 203.583661, 37.934144 ] ], [ [ 203.448738, 37.691853 ], [ 203.310503, 37.690656 ], [ 203.305795, 37.931492 ], [ 203.444566, 37.932797 ] ], [ [ 203.310198, 37.690516 ], [ 203.172649, 37.689470 ], [ 203.167115, 37.929800 ], [ 203.305525, 37.931334 ] ], [ [ 203.172353, 37.689430 ], [ 203.035434, 37.688279 ], [ 203.029517, 37.927979 ], [ 203.166986, 37.929762 ] ], [ [ 204.277381, 37.440086 ], [ 204.141268, 37.438738 ], [ 204.141911, 37.676402 ], [ 204.278543, 37.678047 ] ], [ [ 204.141040, 37.439143 ], [ 204.004591, 37.437304 ], [ 204.004167, 37.675922 ], [ 204.141944, 37.676882 ] ], [ [ 204.004396, 37.437398 ], [ 203.867096, 37.435669 ], [ 203.865766, 37.674914 ], [ 204.003706, 37.675926 ] ], [ [ 203.867086, 37.435701 ], [ 203.729363, 37.434456 ], [ 203.727117, 37.673714 ], [ 203.865299, 37.674598 ] ], [ [ 203.728854, 37.434652 ], [ 203.591691, 37.433480 ], [ 203.587925, 37.672731 ], [ 203.726643, 37.673864 ] ], [ [ 203.591466, 37.433418 ], [ 203.453670, 37.432517 ], [ 203.449024, 37.671193 ], [ 203.587713, 37.672525 ] ], [ [ 203.453485, 37.432768 ], [ 203.316515, 37.431951 ], [ 203.311027, 37.670452 ], [ 203.448894, 37.671774 ] ], [ [ 203.315979, 37.432110 ], [ 203.179317, 37.430903 ], [ 203.173028, 37.669409 ], [ 203.310436, 37.670354 ] ], [ [ 203.178990, 37.430980 ], [ 203.043038, 37.430347 ], [ 203.035800, 37.668169 ], [ 203.172745, 37.669155 ] ] ] } } ], "Fgroups": [ { "Name": { "name": "Name", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "G1" }, "Index": { "name": "Index", "datatype": "int", "ucd": "meta.record;meta.dataset", "value": 1 }, "NFields": { "name": "NFields", "datatype": "int", "ucd": "meta.number;meta.dataset", "value": 8 }, "NAxis": { "name": "NAxis", "datatype": "int", "ucd": "pos.wcs.naxis", "value": 2 }, "Lng_Axis": { "name": "Lng_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.ra", "value": 0 }, "Lat_Axis": { "name": "Lat_Axis", "datatype": "int", "ucd": "meta.id;pos.eq.de", "value": 1 }, "Field_Coordinates": { "name": "Field_Coordinates", "datatype": "double array", "ucd": "pos.eq.ra;pos.eq.dec;obs.field", "value": [ 203.650923, 37.915512 ] }, "Pixel_Scale": { "name": "Pixel_Scale", "datatype": "float array", "ucd": "instr.pixel;obs.field;stat.mean", "value": [ 0.185996, 0.185996 ] }, "Max_Radius": { "name": "Max_Radius", "datatype": "float", "ucd": "phys.size.radius;obs.field", "value": 44.320241 }, "AstRef_Catalog": { "name": "AstRef_Catalog", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "2MASS" }, "AstRef_Band": { "name": "AstRef_Band", "datatype": "string", "ucd": "instr.bandpass", "value": "J" }, "AstromSigma_Internal": { "name": "AstromSigma_Internal", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "value": [ 0.094020, 0.129691 ] }, "AstromCorr_Internal": { "name": "AstromCorr_Internal", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": -0.008094 }, "AstromChi2_Internal": { "name": "AstromChi2_Internal", "datatype": "float", "ucd": "stat.fit.chi2", "value": 91.586733 }, "AstromNDets_Internal": { "name": "AstromNDets_Internal", "datatype": "int", "ucd": "meta.number;src", "value": 206040.000000 }, "AstromSigma_Internal_HighSN": { "name": "AstromSigma_Internal_HighSN", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "value": [ 0.080578, 0.116693 ] }, "AstromCorr_Internal_HighSN": { "name": "AstromCorr_Internal_HighSN", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": -0.080031 }, "AstromChi2_Internal_HighSN": { "name": "AstromChi2_Internal_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 448.308571 }, "AstromNDets_Internal_HighSN": { "name": "AstromNDets_Internal_HighSN", "datatype": "int", "ucd": "meta.number;src", "value": 21893.000000 }, "AstromOffset_Reference": { "name": "AstromOffset_Reference", "datatype": "float array", "ucd": "arith.diff;pos.eq;obs.field", "value": [ 0.001637, -0.010915 ] }, "AstromSigma_Reference": { "name": "AstromSigma_Reference", "datatype": "float array", "ucd": "stat.stdev;pos.eq;obs.field", "value": [ 0.305161, 0.337931 ] }, "AstromCorr_Reference": { "name": "AstromCorr_Reference", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.044493 }, "AstromChi2_Reference": { "name": "AstromChi2_Reference", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.197033 }, "AstromNDets_Reference": { "name": "AstromNDets_Reference", "datatype": "int", "ucd": "meta.number;src", "value": 954 }, "AstromOffset_Reference_HighSN": { "name": "AstromOffset_Reference_HighSN", "datatype": "float array", "ucd": "arith.diff;pos.eq;obs.field", "value": [ 0.002211, -0.010815 ] }, "AstromSigma_Reference_HighSN": { "name": "AstromSigma_Reference_HighSN", "datatype": "float array", "ucd": "stat.stDev;pos.eq;obs.field", "value": [ 0.304169, 0.335807 ] }, "AstromCorr_Reference_HighSN": { "name": "AstromCorr_Reference_HighSN", "datatype": "float", "ucd": "stat.correlation;pos.eq;obs.field", "value": 0.044727 }, "AstromChi2_Reference_HighSN": { "name": "AstromChi2_Reference_HighSN", "datatype": "float", "ucd": "stat.fit.chi2", "value": 3.169405 }, "AstromNDets_Reference_HighSN": { "name": "AstromNDets_Reference_HighSN", "datatype": "int", "ucd": "meta.number;src", "value": 949 }, "NPhotInstru": { "name": "NPhotInstru", "datatype": "int", "ucd": "meta.number;meta.em", "value": 1 }, "PhotInstru_Name": { "name": "PhotInstru_Name", "datatype": "string array", "ucd": "meta.id;instr.bandpass", "value": [ "P1" ] }, "PhotSigma_Internal": { "name": "PhotSigma_Internal", "datatype": "float array", "ucd": "stat.error;phot.mag", "value": [ 0.128496 ] }, "PhotChi2_Internal": { "name": "PhotChi2_Internal", "datatype": "float array", "ucd": "stat.chi2;phot.mag", "value": [ 8.600261 ] }, "PhotNDets_Internal": { "name": "PhotNDets_Internal", "datatype": "int array", "ucd": "meta.number;src", "value": [ 204250 ] }, "PhotSigma_Internal_HighSN": { "name": "PhotSigma_Internal_HighSN", "datatype": "float array", "ucd": "stat.error;phot.mag", "value": [ 0.035560 ] }, "PhotChi2_Internal_HighSN": { "name": "PhotChi2_Internal_HighSN", "datatype": "float array", "ucd": "stat.chi2;phot.mag", "value": [ 38.447707 ] }, "PhotNDets_Internal_HighSN": { "name": "PhotNDets_Internal_HighSN", "datatype": "int array", "ucd": "meta.number;src", "value": [ 21862 ] }, "PhotSigma_Reference": { "name": "PhotSigma_Reference", "datatype": "float array", "ucd": "stat.error;phot.mag", "value": [ 0.000000 ] }, "PhotChi2_Reference": { "name": "PhotChi2_Reference", "datatype": "float array", "ucd": "stat.chi2;phot.mag", "value": [ 0.000000 ] }, "PhotNDets_Reference": { "name": "PhotNDets_Reference", "datatype": "int array", "ucd": "meta.number;src", "value": [ 0.000000 ] }, "PhotSigma_Reference_HighSN": { "name": "PhotSigma_Reference_HighSN", "datatype": "float array", "ucd": "stat.error;phot.mag", "value": [ 0.000000 ] }, "PhotChi2_Reference_HighSN": { "name": "PhotChi2_Reference_HighSN", "datatype": "float array", "ucd": "stat.chi2;phot.mag", "value": [ 0.000000 ] }, "PhotNDets_Reference_HighSN": { "name": "PhotNDets_Reference_HighSN", "datatype": "int array", "ucd": "meta.number;src", "value": [ 0 ] }, "FgroupsPlot": { "name": "FgroupsPlot", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "fgroups_1.png" }, "Chi2Plot": { "name": "Chi2Plot", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "astr_chi2_1.png" }, "IntErr1DimPlot": { "name": "IntErr1DimPlot", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "astr_interror1d_1.png" }, "IntErr2DimPlot": { "name": "IntErr2DimPlot", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "astr_interror2d_1.png" }, "RefErr1DimPlot": { "name": "RefErr1DimPlot", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "astr_referror1d_1.png" }, "RefErr2DimPlot": { "name": "RefErr2DimPlot", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "astr_referror2d_1.png" }, "PhotErrPlot": { "name": "PhotErrPlot", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "psphot_error_1.png" } } ], "AstroInstruments": [ { "Name": { "name": "Name", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "A1" }, "Index": { "name": "Index", "datatype": "int", "ucd": "meta.record;meta.dataset", "value": 1 }, "NFields": { "name": "NFields", "datatype": "int", "ucd": "meta.number;meta.dataset", "value": 8 }, "MagZeroPoint_Output": { "name": "MagZeroPoint_Output", "datatype": "float", "ucd": "astr.mag;astr.calib;arith.zp", "value": 36.000000 }, "NKeys": { "name": "NKeys", "datatype": "int", "ucd": "meta.number", "value": 2 }, "Keys": { "name": "Keys", "datatype": "string array", "ucd": "meta.note", "value": [ "FILTER  = \'i.MP9701\'             ", "QRUNID  = \'04AQ04  \'             " ] }, "DistPlot": { "name": "DistPlot", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "distort_1.png" } } ], "PhotInstruments": [ { "Name": { "name": "Name", "datatype": "string", "ucd": "meta.id;meta.dataset", "value": "P1" }, "Index": { "name": "Index", "datatype": "int", "ucd": "meta.record;meta.dataset", "value": 1 }, "NFields": { "name": "NFields", "datatype": "int", "ucd": "meta.number;meta.dataset", "value": 8 }, "MagZeroPoint_Output": { "name": "MagZeroPoint_Output", "datatype": "float", "ucd": "phot.mag;phot.calib;arith.zp", "value": 0.000000 }, "NKeys": { "name": "NKeys", "datatype": "int", "ucd": "meta.number", "value": 1 }, "Keys": { "name": "Keys", "datatype": "string array", "ucd": "meta.note", "value": [ "FILTER  = \'i.MP9701\'             " ] } } ], "Warnings": [ { "Date": { "name": "Date", "datatype": "string", "ucd": "meta;time.end", "value": "2018-10-08" }, "Time": { "name": "Time", "datatype": "string", "ucd": "meta;time.end", "value": "09:27:48" }, "Text": { "name": "Text", "datatype": "string", "ucd": "meta", "value": "This executable has been compiled using a version of the ATLAS library without sPerformance will be degraded." } } ], "CommandLine": ".\/src\/scamp tests\/extra\/744331p.cat tests\/extra\/744332p.cat tests\/extra\/744333p.cat tests\/extra\/744334p.cat tests\/extra\/744335p.cat tests\/extra\/744336p.cat tests\/extra\/744337p.cat tests\/extra\/744338p.cat ", "Configuration": [ { "datatype": "string", "name": "AHEADER_GLOBAL", "value": "scamp.ahead" }, { "datatype": "string array", "name": "AHEADER_NAME", "value": [ ] }, { "datatype": "string", "name": "AHEADER_SUFFIX", "value": ".ahead" }, { "datatype": "string", "name": "AIRMASS_KEY", "value": "AIRMASS" }, { "datatype": "float", "name": "ASTR_ACCURACY", "value": 89128.960938 }, { "datatype": "int", "name": "ASTR_FLAGSMASK", "value": 252 }, { "datatype": "int", "name": "ASTR_IMAFLAGSMASK", "value": 0 }, { "datatype": "string", "name": "ASTRACCURACY_KEY", "value": "ASTRACCU" }, { "datatype": "string", "name": "ASTRACCURACY_TYPE", "value": "SIGMA-PIXEL" }, { "datatype": "float", "name": "ASTRCLIP_NSIGMA", "value": 0.000000 }, { "datatype": "string", "name": "ASTREF_BAND", "value": "DEFAULT" }, { "datatype": "string", "name": "ASTREF_CATALOG", "value": "2MASS" }, { "datatype": "float", "name": "ASTREF_WEIGHT", "value": 0.000000 }, { "datatype": "string array", "name": "ASTREFCAT_NAME", "value": [ "astrefcat.cat" ] }, { "datatype": "string array", "name": "ASTREFCENT_KEYS", "value": [ "X_WORLD", "Y_WORLD" ] }, { "datatype": "string array", "name": "ASTREFERR_KEYS", "value": [ "ERRA_WORLD", "ERRB_WORLD", "ERRTHETA_WORLD" ] }, { "datatype": "string", "name": "ASTREFMAG_KEY", "value": "MAG" }, { "datatype": "string", "name": "ASTREFMAGERR_KEY", "value": "MAGERR" }, { "datatype": "string", "name": "ASTREFOBSDATE_KEY", "value": "OBSDATE" }, { "datatype": "float array", "name": "ASTREFMAG_LIMITS", "value": [ -99.000000, 99.000000 ] }, { "datatype": "string array", "name": "ASTRINSTRU_KEY", "value": [ "FILTER", "QRUNID" ] }, { "datatype": "string array", "name": "CENTROID_KEYS", "value": [ "XWIN_IMAGE", "YWIN_IMAGE" ] }, { "datatype": "string array", "name": "CENTROIDERR_KEYS", "value": [ "ERRAWIN_IMAGE", "ERRBWIN_IMAGE", "ERRTHETAWIN_IMAGE" ] }, { "datatype": "string array", "name": "CHECKIMAGE_NAME", "value": [ "check.fits" ] }, { "datatype": "string array", "name": "CHECKIMAGE_TYPE", "value": [ "NONE" ] }, { "datatype": "boolean", "name": "CHECKPLOT_ANTIALIAS", "value": true }, { "datatype": "string", "name": "CHECKPLOT_CKEY", "value": "SCAMPCOL" }, { "datatype": "string array", "name": "CHECKPLOT_DEV", "value": [ "PNG" ] }, { "datatype": "string array", "name": "CHECKPLOT_NAME", "value": [ "fgroups", "distort", "astr_interror2d", "astr_interror1d", "astr_referror2d", "astr_referror1d", "astr_chi2", "psphot_error" ] }, { "datatype": "int array", "name": "CHECKPLOT_RES", "value": [ 0 ] }, { "datatype": "string array", "name": "CHECKPLOT_TYPE", "value": [ "FGROUPS", "DISTORTION", "ASTR_INTERROR2D", "ASTR_INTERROR1D", "ASTR_REFERROR2D", "ASTR_REFERROR1D", "ASTR_CHI2", "PHOT_ERROR" ] }, { "datatype": "boolean", "name": "COMPUTE_PARALLAXES", "value": false }, { "datatype": "boolean", "name": "COMPUTE_PROPERMOTIONS", "value": false }, { "datatype": "boolean", "name": "CORRECT_COLOURSHIFTS", "value": false }, { "datatype": "float", "name": "CROSSID_RADIUS", "value": 0.000000 }, { "datatype": "string", "name": "DGEOMAP_NAME", "value": "dgeo.fits" }, { "datatype": "int", "name": "DGEOMAP_NNEAREST", "value": 21 }, { "datatype": "int", "name": "DGEOMAP_STEP", "value": 2 }, { "datatype": "string array", "name": "DISTORT_KEYS", "value": [ "XWIN_IMAGE", "YWIN_IMAGE" ] }, { "datatype": "int array", "name": "DISTORT_GROUPS", "value": [ 1, 1 ] }, { "datatype": "int array", "name": "DISTORT_DEGREES", "value": [ 3 ] }, { "datatype": "float", "name": "ELLIPTICITY_MAX", "value": 0.000000 }, { "datatype": "string", "name": "EXPOTIME_KEY", "value": "EXPTIME" }, { "datatype": "string", "name": "EXTINCT_KEY", "value": "PHOT_K" }, { "datatype": "int", "name": "FIXFOCALPLANE_NMIN", "value": 1 }, { "datatype": "int", "name": "FLAGS_MASK", "value": 240 }, { "datatype": "int", "name": "FOCDISTORT_DEGREE", "value": 1 }, { "datatype": "string", "name": "FULLOUTCAT_NAME", "value": "full.cat" }, { "datatype": "string", "name": "FULLOUTCAT_TYPE", "value": "NONE" }, { "datatype": "float array", "name": "FWHM_THRESHOLDS", "value": [ 0.000000, 100.000000 ] }, { "datatype": "string", "name": "HEADER_SUFFIX", "value": ".head" }, { "datatype": "string", "name": "HEADER_TYPE", "value": "NORMAL" }, { "datatype": "string array", "name": "HEADER_NAME", "value": [ ] }, { "datatype": "int", "name": "IMAFLAGS_MASK", "value": 0 }, { "datatype": "boolean", "name": "INCLUDE_ASTREFCATALOG", "value": true }, { "datatype": "string", "name": "MAGZERO_KEY", "value": "PHOT_C" }, { "datatype": "float array", "name": "MAGZERO_OUT", "value": [ 0.000000 ] }, { "datatype": "float array", "name": "MAGZERO_INTERR", "value": [ 0.010000 ] }, { "datatype": "float array", "name": "MAGZERO_REFERR", "value": [ 0.030000 ] }, { "datatype": "boolean", "name": "MATCH", "value": true }, { "datatype": "boolean", "name": "MATCH_FLIPPED", "value": false }, { "datatype": "int", "name": "MATCH_NMAX", "value": 0 }, { "datatype": "float", "name": "MATCH_RESOL", "value": 0.000000 }, { "datatype": "string", "name": "MERGEDOUTCAT_NAME", "value": "merged.cat" }, { "datatype": "string", "name": "MERGEDOUTCAT_TYPE", "value": "NONE" }, { "datatype": "string array", "name": "MOSAIC_TYPE", "value": [ "UNCHANGED" ] }, { "datatype": "int", "name": "NTHREADS", "value": 4 }, { "datatype": "float", "name": "PHOT_ACCURACY", "value": -518969491456.000000 }, { "datatype": "int", "name": "PHOT_FLAGSMASK", "value": 252 }, { "datatype": "int", "name": "PHOT_IMAFLAGSMASK", "value": 0 }, { "datatype": "float", "name": "PHOTCLIP_NSIGMA", "value": 0.000000 }, { "datatype": "string", "name": "PHOTFLUX_KEY", "value": "FLUX_AUTO" }, { "datatype": "string", "name": "PHOTFLUXERR_KEY", "value": "FLUXERR_AUTO" }, { "datatype": "string array", "name": "PHOTINSTRU_KEY", "value": [ "FILTER" ] }, { "datatype": "string", "name": "PHOTOMFLAG_KEY", "value": "PHOTFLAG" }, { "datatype": "float", "name": "PIXSCALE_MAXERR", "value": 0.000000 }, { "datatype": "float", "name": "POSANGLE_MAXERR", "value": 0.000000 }, { "datatype": "float array", "name": "POSITION_MAXERR", "value": [ 0.016667, 0.016667 ] }, { "datatype": "string array", "name": "PROJECTION_TYPE", "value": [ "SAME" ] }, { "datatype": "string", "name": "REFOUT_CATPATH", "value": "." }, { "datatype": "string array", "name": "REF_SERVER", "value": [ "vizier.unistra.fr" ] }, { "datatype": "float array", "name": "REF_TIMEOUT", "value": [ 10.000000 ] }, { "datatype": "boolean", "name": "SAVE_REFCATALOG", "value": false }, { "datatype": "boolean", "name": "SAVE_DGEOMAP", "value": false }, { "datatype": "float array", "name": "SN_THRESHOLDS", "value": [ 10.000000, 100.000000 ] }, { "datatype": "boolean", "name": "SOLVE_ASTROM", "value": true }, { "datatype": "boolean", "name": "SOLVE_PHOTOM", "value": true }, { "datatype": "string array", "name": "STABILITY_TYPE", "value": [ "INSTRUMENT" ] }, { "datatype": "string", "name": "VERBOSE_TYPE", "value": "NORMAL" }, { "datatype": "int", "name": "WEIGHTFLAGS_MASK", "value": 255 }, { "datatype": "boolean", "name": "WRITE_XML", "value": true }, { "datatype": "boolean", "name": "WRITE_JSON", "value": true }, { "datatype": "boolean", "name": "WRITE_HTML", "value": true }, { "datatype": "string", "name": "XML_NAME", "value": "scamp.xml" }, { "datatype": "string", "name": "JSON_NAME", "value": "scamp.json" }, { "datatype": "string", "name": "HTML_NAME", "value": "scamp.html" }, { "datatype": "string", "name": "HTML_TPL", "value": ".\/html\/scamp.html.tpl" }, { "datatype": "string", "name": "XSL_URL", "value": "file:\/\/\/scamp\/scamp.xsl" } ] }');
</script>
